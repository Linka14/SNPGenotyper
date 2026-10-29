
# Shared functions of the pipeline, loaded at the start of a stage with source("Code/R_scripts/asp_functions.R", local = TRUE). Updating
# functions is now a lot easier because every script shares this script, its finally unified.

### Timing function for speed optimization analysis, needs timing_log <- list() in the stage before the first call
log_time <- function(label, t_start) { # takes a label and start time, computes elapsed seconds and appends the row 
  elapsed <- as.numeric(Sys.time() - t_start, units = "secs")
  timing_log[[length(timing_log) + 1]] <<- data.frame(step = label, seconds = elapsed) 
  cat(sprintf("[TIMING] %s: %.2f sec\n", label, elapsed))
}

### This helper function gives the actual name of the sample not just barcode01 or 01.
## needs barcode_map (sample sheet filtered to the dataset,with a barcode_num column) in every .Rmd stage that uses it
get_actual_sample_name <- function(bc) {
  bc_num <- sub("^barcode", "", bc)
  hit <- barcode_map$sample_id[barcode_map$barcode_num == bc_num]
  if (length(hit) == 1) return(hit)
  warning(sprintf("no unique sample name for %s, keeping the barcode name", bc))
  bc
}

### This function creates folder names for every run, for both ASP and Bayes runs. They can be build with different settings, this function put
## every results in its separate folder
build_run_names <- function(input_type, use_anchors, anchor_length, threshold_type, error_rate, min_score_fraction,
                            prior_type = "flat", prior_het_share = NA, prior_population = NA, prior_floor = NA,
                            probe_type = "asp_primers", ancestry_prior_shape = 1, phi_hom_structure = "shared") {
  asp_settings_label <- paste0(input_type,
                               if (probe_type == "reference") "_reference" # reference probes get their own folders, and have no separate anchor (centred window)
                               else if (use_anchors) paste0("_anchor", anchor_length) else "_noanchor",
                               if (threshold_type == "errors") paste0("_errors", error_rate) else paste0("_fraction", min_score_fraction))
  prior_label <- switch(prior_type,
                        flat = "priorflat",
                        het_share = paste0("priorhet", prior_het_share),
                        hwe = paste0("priorhwe_", prior_population, "_floor", prior_floor),
                        ancestry = paste0("priorancestry_", prior_population, "_FreqDOG_floor", prior_floor, "_beta", ancestry_prior_shape),
                        stop("prior_type must be \"flat\", \"het_share\", \"hwe\" or \"ancestry\"")) # used in folder names, so runs with different priors never overwrite each other
  if (phi_hom_structure != "shared") prior_label <- paste0(prior_label, "_phihom", phi_hom_structure) # homozygote overdispersion structure, only added when not the default
  list(asp_run_name = paste0("Identified_ASP_genotypes_", asp_settings_label),
       prior_label = prior_label,
       bayes_run_name = paste0("Base_Mixture_STAN_Hierarchical_ASP_genotypes_", asp_settings_label, "_", prior_label))
}

### Helper function for making .fastq files names, input_type picks the padded or the unpadded ObiTools4 identified .fastq files.
identified_file_names_f <- function(bc) {
  if (input_type == "padded") {
    sprintf("Data/Processed_data/%s/Demultiplexed_sequence_data/Identified_padded_sample_sequences/barcode%s_padded.fastq", dataset_name, bc)
  } else {
    sprintf("Data/Processed_data/%s/Demultiplexed_sequence_data/Identified_sample_sequences_Obi4/barcode%s_identified.fastq", dataset_name, bc)
  }
}

### Helper function for primer trimming, used wih ASP primer
trim_to_common <- function(a, b, input_type, cut = 0L) {
  n <- min(nchar(a), nchar(b)) # Grabs the length of the smaller shorter primer, both primers will be trimmed with n
  if (input_type == "unpadded") n <- n - cut # unpadded reads: ObiTools trims the STA primer, which also removes the ASP bases that overlap the STA's 3' end
  c(substr(a, nchar(a) - n + 1, nchar(a)), substr(b, nchar(b) - n + 1, nchar(b))) # cuts out first part of the string
}

### Helper function for computing how many bases at the 3' end of the STA primer are the same as the 5' start of the ASP primer
sta_overlap_f <- function(sta, asp) {
  if (is.na(sta) || is.na(asp)) return(0L) # no STA primer in the panel: nothing to overlap
  k <- min(nchar(sta), nchar(asp))
  while (k > 0 && substr(sta, nchar(sta) - k + 1, nchar(sta)) != substr(asp, 1, k)) k <- k - 1 # shortens the tested overlap until both ends match
  k
}

### Primer panel trimming and anchoring preparation: reads the panel, trims ASP1/ASP2 to the same length (unpadded: also cuts the bases ObiTools 
## removes with the STA primer) and adds the anchor bases after the SNP from the reference.
# Returns the finished panel, one row per locus. Used by the ASP genotyping and the visualization Rmd.
prepare_primer_panel <- function(primer_panel_path, reference_path, input_type, use_anchors, anchor_length, unpadded_cut_margin = 4L,
                                 probe_type = "asp_primers", # "asp_primers" or "reference" 
                                 snp_position = 500L) { # for reference probes only, SNP position in every sequence
  stopifnot("probe_type must be \"asp_primers\" or \"reference\"" = probe_type %in% c("asp_primers", "reference"),
            "reference probes need padded reads (input_type: \"padded\")" = probe_type == "asp_primers" || input_type == "padded")
    
  #this imports the excel data, keeps only the columns I need, I changed the names
  panel <- read_excel(primer_panel_path, sheet = "Final panel", .name_repair = "unique_quiet")
  panel <- panel[, c("SNP name", "Allele 1", "Allele 2", "ASP1_SEQ (forward allele1)", "ASP2_SEQ (forward allele2)", "STA_SEQ (forward STA step)", "LSP_SEQ (reverse)")]
  names(panel) <- c("locus", "allele1", "allele2", "asp1_seq", "asp2_seq", "sta_seq", "lsp_seq")
  panel <- panel[!is.na(panel$asp1_seq) & !is.na(panel$asp2_seq), ]
  longest_asp_primer <- max(nchar(trimws(c(panel$asp1_seq, panel$asp2_seq)))) # before any trimming, lentgh of every reference probe (SNP included)
  
  asp_common <- mapply(trim_to_common, trimws(panel$asp1_seq), trimws(panel$asp2_seq), MoreArgs = list(input_type = "padded"))[1, ] #ASP at common length, before any unpadded cut
  panel$sta_overlap <- mapply(sta_overlap_f, trimws(panel$sta_seq), asp_common, USE.NAMES = FALSE) # ASP bases ObiTools removes at each locus
  
  # each locus is cut by its own overlap (+ margin), padded is not cut
  trimmed <- mapply(function(a, b, cut) trim_to_common(a, b, input_type = input_type, cut = cut),
                    trimws(panel$asp1_seq), trimws(panel$asp2_seq), panel$sta_overlap + unpadded_cut_margin) #trimws removes trailing white spaces
  panel$asp1_seq <- trimmed[1, ] #saves the trimmed primers
  panel$asp2_seq <- trimmed[2, ]
  
  # this part creates anchors for the ASP1 and ASP2 primers on the 3' end right after the SNP position. It needs them because after lowering the
  # indel gap penalty (same penalty as minimap2 6 for gap opening and 2 for gap extension) it can in some cases replace the random third base 
  # (not matching either of the two alleles) at the SNP position for a gap and grabs the next base after SNP position which can coincidently match 
  # the one of the two alleles, creating a false call. The added 3 base anchor hopefully fixes that using those bases to really anchor that SNP position.
  # The anchors are computed on every run with anchor_length bases (also for unpadded input), use_anchors decides if they are used.
  reference <- readDNAStringSet(reference_path)
  names(reference) <- sub("::.*$", "", names(reference)) # grabs and replaces everthing after the name with blank ""
  stopifnot(all(panel$locus %in% names(reference))) # every panel locus needs a reference sequence
  
  get_anchor <- function(locus, asp_seq, lsp_seq) { #inputs are locus name, trimmed primer sequence and the LSP primer
    ref <- reference[[locus]] # grabs the locuses reference sequence 1000 bases long
    body <- DNAString(substr(asp_seq, 1, nchar(asp_seq) - 1)) # primer without its SNP base, asp1 and asp2 are identical
    lsp_rc <- reverseComplement(DNAString(trimws(lsp_seq))) # LSP binds the other strand, so on the ASP strand it appears reverse complemented, downstream of the SNP
    for (strand in c("+", "-")) { #first searched for refenrece as its stored (+) and for the revese compliment (-). Some of the primers were designed on the opposite trand for reference
      s <- if (strand == "+") ref else reverseComplement(ref) # s is the sequence searched
      hit <- matchPattern(body, s, max.mismatch = 2) # tolerate up to 2 primer/reference differences
      if (length(hit) == 1) { #onlt one single match is accepted, it has to be unambigious, if its 0 or 2 it searched the rc
        snp <- end(hit) + 1 # SNP = base right after the primer body
        if (snp + anchor_length > length(s)) break  #if its longer than the refernce it stops searcing, return NAs. theoreticall t cant happen (1000bs)
        lsp_hit <- matchPattern(lsp_rc, s, max.mismatch = 2) # where the LSP primer starts on this strand
        lsp_hit <- lsp_hit[start(lsp_hit) > snp] # only a hit downstream of the SNP is the amplicon end
        return(c(anchor = as.character(subseq(s, snp + 1, snp + anchor_length)), #cuts the anchor_length bases after SNP
                 ref_snp_base = as.character(subseq(s, snp, snp)), #refenrece base at snp position
                 ref_snp_pos = if (strand == "+") snp else length(s) - snp + 1, #converts the positon if its on rc
                 ref_strand = strand, #records orinetation
                 bases_to_lsp = if (length(lsp_hit) == 1) start(lsp_hit) - snp - 1 else NA)) # bases between the SNP and the LSP primer = what an unpadded read still has after the SNP
      }
    }
    c(anchor = NA, ref_snp_base = NA, ref_snp_pos = NA, ref_strand = NA, bases_to_lsp = NA) # not found
  }
  
  anchors <- t(mapply(get_anchor, panel$locus, asp_common, panel$lsp_seq)) #call the function once per panel row for both priemrs; uses the uncut primers so the reference search also works for unpadded input
  panel$anchor <- anchors[, "anchor"] #copies everyhing into the panel
  panel$ref_snp_base <- anchors[, "ref_snp_base"]
  panel$ref_snp_pos <- as.integer(anchors[, "ref_snp_pos"])
  panel$ref_strand <- anchors[, "ref_strand"]
  panel$bases_to_lsp <- as.integer(anchors[, "bases_to_lsp"])
  
  # unpadded reads end where the LSP primer starts, so the anchor is shortened per locus to the bases the read still has after the SNP (padded reads keep the full anchor)
  if (input_type == "unpadded") panel$anchor <- substr(panel$anchor, 1, pmin(anchor_length, panel$bases_to_lsp, na.rm = TRUE))
  if (!use_anchors) panel$anchor <- "" # no anchors get empty string
  panel$anchor[is.na(panel$anchor)] <- ""
  
  # Reference probes are cut from the reference + strand, centred on snp_position, as long as the longest ASP primer of the panel (SNP included).
  # Centred, because one alignment has no second primer to disagree with: a short side next to the SNP lets the aligner do the hopping again. use_anchors and anchor_length are not used here, 
  #the bases after the SNP are part of the window. asp1_seq and asp2_seq are two copies that differ only at the SNP (panel allele1 and allele2, 
  #which are written for the + strand), the scorer merges them into one probe with the IUPAC code at the SNP. Which way each read runs is found in 
  # the scorer, so no strand information is needed here.
  if (probe_type == "reference") {
    left_length <- longest_asp_primer %/% 2L # reference bases before the SNP (15 for a 30-base primer)
    right_length <- longest_asp_primer - 1L - left_length # reference bases after the SNP (14), stored in the anchor column
    allele1 <- toupper(trimws(panel$allele1))
    allele2 <- toupper(trimws(panel$allele2))
    reference_lengths <- width(reference[panel$locus])
    stopifnot("every allele must be one of A, C, G, T" = allele1 %in% c("A", "C", "G", "T") & allele2 %in% c("A", "C", "G", "T"),
              "allele1 and allele2 must differ" = allele1 != allele2,
              "the reference probe window runs past the end of a reference sequence" = snp_position - left_length >= 1 & snp_position + right_length <= reference_lengths)
    window_text <- as.character(subseq(reference[panel$locus], snp_position - left_length, snp_position + right_length)) # reference + strand around the SNP
    reference_base <- substr(window_text, left_length + 1L, left_length + 1L)
    stopifnot("a reference window contains a base other than A, C, G or T" = !grepl("[^ACGT]", window_text),
              "the reference base at snp_position is not one of the two alleles" = reference_base == allele1 | reference_base == allele2)
    panel$asp1_seq <- paste0(substr(window_text, 1, left_length), allele1) # copy with allele 1 at the SNP
    panel$asp2_seq <- paste0(substr(window_text, 1, left_length), allele2) # copy with allele 2 at the SNP
    panel$anchor <- substr(window_text, left_length + 2L, nchar(window_text)) # reference bases after the SNP
  }
  
  print(data.frame(locus = panel$locus, sta_overlap = panel$sta_overlap,
                   length_before_cut = nchar(asp_common), length_after_cut = nchar(panel$asp1_seq)), row.names = FALSE)
  panel
}

### Read scoring helpers for score_reads_against_primers_2pwalign_v8
## match_lookup and mismatch_cost_lookup are built in the Rmd from primer_quality: match_lookup <- build_match_lookup(primer_quality) ...
# Match score of one matching base for every read quality from 0 to 93 against a primer base quality (primer_quality).
# match_lookup[Q + 1] = score at Phred Q. Result it match_lookup with 94 scores for each quality maximum at 1.98.
build_match_lookup <- function(primer_quality = 22L) {
  phred_value <- 0:93 #every possible phred value
  score(pairwiseAlignment(DNAStringSet(rep("A", length(phred_value))), DNAString("A"), type = "local-global", #matches string set of 94 As to A
                          patternQuality = PhredQuality(vapply(phred_value, function(q) intToUtf8(q + 33), "")), #turns each phred_value to .fastq QS (ascii)
                          subjectQuality = PhredQuality(intToUtf8(primer_quality + 33)))) #each base aligned to chosen primer quality
}

### Calucaltes score difference between a matching and a mismatching read base of specific Q against a primer base with primer_quality.
## It calcualtes what a read loses when one base mismatches instead of matches at specific quality Q of a primer (higher the primer Q, higher the loss). Its used in error_cost adn for SNP filter
mismatch_cost_f <- function(q, primer_quality) {
  score_vs_A <- function(base) {
    score(pairwiseAlignment(DNAString(base), DNAString("A"), type = "global",
                                              patternQuality = PhredQuality(intToUtf8(q + 33)),
                                              subjectQuality = PhredQuality(intToUtf8(primer_quality + 33))))
  }
  score_vs_A("A") - score_vs_A("C") #score of a matching base at certain Q - score of a mismatch at Q.
}

### Function for computing exact perfect score per read, which is the  sum of match scores of the primer_len read bases starting at pattern_start. 
## It substitutes the previous calcualtion of mean perfect scores which werent accurate.
perfect_from_read_qualities <- function(read_quals, pattern_start, primer_len, match_lookup) {
  read_len <- nchar(read_quals) #length of each reads quality string
  s <- pmax(1L, pmin(pattern_start, read_len - primer_len + 1L)) #pmax stops it starting before base 1 if the read is shorter than primer, pmin pulls the window back if the primer would tretch over the end of the read
  region <- substring(read_quals, s, s + primer_len - 1L) #region of priemrs quality characters
  full <- nchar(region) == primer_len #false for reads shorter than the primer
  out <- numeric(length(read_quals)) #empty result vector
  
  if (any(full)) { #for reads that are long enough (true)
    codes <- utf8ToInt(paste(region[full], collapse = "")) #glues the windows together and converts them to Q numbers, only one call
    out[full] <- rowSums(matrix(match_lookup[codes - 32L], ncol = primer_len, byrow = TRUE)) #creates a matrix with primer_lentgh column and each read in a row with its matchlookup scores.
    #It fill it row by row, which guarantees each read is in one row (shorter one are treated separately). Rows sums calcualtes perfect score
  }
  if (any(!full)) { #for shorter reads (false)
    out[!full] <- vapply(region[!full], function(r) sum(match_lookup[utf8ToInt(r) - 32L]), numeric(1), USE.NAMES = FALSE) #calcualtes one at a time
  }
  out
}

### This function gathers what went wrong in each read alignemnt, mismatcheds (how many, hwo much they cost, ...) and indels
mismatch_summary <- function(alignment, n, mismatch_cost_lookup, high_quality = 20, # takes object from pairwiseAlignemnt, number of reads, coost of mismatches, Phred of a confident mismatch
                             skip_subject_position = NULL) { # takes object from pairwiseAlignemnt, number of reads, coost of mismatches, Phred of a confident mismatch
  summary_table <- data.frame(mismatch_cost = rep(0, n), #builds the empty results data frame, one read per row
                              max_mismatch_quality = NA_real_,
                              high_quality_mismatches = 0L,
                              mean_mismatch_quality = NA_real_)
  mt <- mismatchTable(alignment) #object from pairwiseAlignment
  if (!is.null(skip_subject_position)) mt <- mt[mt$SubjectStart != skip_subject_position[mt$PatternId], ] # reference probes: an allele at the IUPAC SNP code is a half match, not a mismatch
  if (nrow(mt) > 0) {
    quality <- utf8ToInt(paste(as.character(mt$PatternQuality), collapse = "")) - 33 # Phred of every mismatched read base
    read_id <- mt$PatternId
    cost_per_read <- tapply(mismatch_cost_lookup[quality + 1], read_id, sum) # cost of each mismatch at its own quality, summed per read
    rows <- as.integer(names(cost_per_read)) # reads that have mismatches
    summary_table$mismatch_cost[rows] <- cost_per_read
    summary_table$max_mismatch_quality[rows] <- tapply(quality, read_id, max) #eachtapply groups by the same read_id
    summary_table$high_quality_mismatches[rows] <- tapply(quality >= high_quality, read_id, sum)
    summary_table$mean_mismatch_quality[rows] <- tapply(quality, read_id, mean)
  }
  indels <- nindel(alignment) #runs outside if, indel can exist when read has no mismatches. It runs for every locus as well
  summary_table$indel_events <- insertion(indels)[, "Length"] + deletion(indels)[, "Length"] #saves howmany indel events hapenned in a read
  summary_table
}

### All aligned reads / primers are glued into one long column vector; cumsum() counts real (non "-") bases, and subtracting each read's count
## before its first column restarts the counters per read. Function takes, aligned strings with "-" for gas (same width per read), start of the 
# alignemnt and primer position of the SNP and checks if there is an indel close to the SNP postion.
map_snp_positions <- function(read_aln, primer_aln, read_start, target, window = 3L, primer_start = 1L) { 
  n <- length(read_aln)
  width <- nchar(primer_aln) # alignment columns per read, gaps included
  last <- cumsum(width) # global index of each read's last column
  first <- last - width + 1L # first column
  id <- rep.int(seq_len(n), width) # wcreates id for every column for a read based on vector of width (1, 1, 1, 1, 1, 2, 2, 2, 2, ...)
  
  read <- strsplit(paste(read_aln, collapse = ""), "", fixed = TRUE)[[1]] #glues aligned strings into one long and splits by each character
  primer <- strsplit(paste(primer_aln, collapse = ""), "", fixed = TRUE)[[1]] #same here
  read_gap_check <- read != "-" #true when there is a real base, false when a gap, turns all the column into either T or F
  primer_gap_check <- primer != "-"
  
  primer_cumulative <- cumsum(primer_gap_check) #counts primer bases form first to alst column, where T is 1 and F (gap) is 0. When there is a gap the column number and cmualtove primer sum dont match up
  primer_pos <- primer_cumulative - (primer_cumulative[first] - primer_gap_check[first])[id] +
    (primer_start - 1L)[id] #absolute primer position, for primers whos first bases hang of the read
  read_cumulative <- cumsum(read_gap_check) 
  read_pos <- read_cumulative - (read_cumulative[first] - read_gap_check[first])[id] # same for read bases
  
  hit <- which(primer_gap_check & primer_pos == target[id]) # the SNP column (exactly one per read)
  col <- rep(NA_integer_, n)
  col[id[hit]] <- hit
  deleted <- !read_gap_check[col] # read has "-" at the SNP column
  read_base <- ifelse(deleted, NA_character_, read[col])
  primer_base <- primer[col]
  
  gap <- !(read_gap_check & primer_gap_check) # a column with a gap on either side
  gap_cumulative <- cumsum(gap)
  low <- pmax(first, col - window)
  high <- pmin(last, col + window) # window, kept inside the read's own columns
  
  list(read_position = ifelse(deleted, NA_integer_, read_start - 1L + read_pos[col]),
       read_base = read_base,
       primer_base = primer_base,
       deleted_in_read = deleted,
       base_matches = ifelse(deleted, NA, read_base == primer_base),
       indel_near_snp = (gap_cumulative[high] - gap_cumulative[low] + gap[low]) > 0)
}

### The main scorer function which scores all the alignemnts and combines all the previous helper functions, results are used for the genotype caller function
score_reads_against_primers_2pwalign_v8 <- function(read_set, read_quals, read_direction,
                                                    asp1_seq, asp2_seq, allele1, allele2,
                                                    min_score_fraction = 0.65,
                                                    primer_quality = 22L,
                                                    match_lookup = NULL,
                                                    gap_opening = 10,
                                                    gap_extension = 4,
                                                    anchor = "", # panel$anchor[i] or default = "" no anchor
                                                    snp_indel_window = 3L, # an indel three bases near the snp
                                                    require_position_agreement = TRUE,
                                                    input_type = "unpadded", #or padded
                                                    threshold_type = "errors", # "fraction" = min_score_fraction * perfect, "errors" = error allowance
                                                    error_rate = 0.08, # "errors" only: allowed Q20 mismatches per around 12 bases + anchor
                                                    min_errors = 1, # "errors" only: every primer allows at least this many errors of any kind 
                                                    min_snp_quality = 0, #reads whose SNP base quality is below this Phred are flagged (compared directly)
                                                    max_mismatches = Inf, #default is off
                                                    mismatch_cost_lookup = NULL,
                                                    probe_type = "asp_primers") { # reference genome or primers
  
  if (is.null(match_lookup)) match_lookup <- build_match_lookup(primer_quality) #build matchlookup if it wasnt passed
  if (is.null(mismatch_cost_lookup)) mismatch_cost_lookup <- vapply(0:93, function(q) mismatch_cost_f(q, primer_quality), numeric(1))
  if (length(anchor) == 0 || is.na(anchor)) anchor <- "" # if no anchor column or NA it runs run without anchor
  anchor <- toupper(anchor)
  
  # SNP position inside each primer, before the anchor is added (for forward orientation its the last ASP base)
  snp_pos1 <- nchar(asp1_seq)
  snp_pos2 <- nchar(asp2_seq)
  
  asp1 <- DNAString(paste0(asp1_seq, anchor)) # creates primers with anchor and its reverse complements
  asp2 <- DNAString(paste0(asp2_seq, anchor))
  asp1_rc <- reverseComplement(asp1)
  asp2_rc <- reverseComplement(asp2)
  
  read_quality <- PhredQuality(read_quals) #object for pwalign of quality character strings
  primer_quality_obj <- PhredQuality(as.integer(primer_quality)) #turns preasigned primer quality to single phredscore
  n <- length(read_set) # number of reads a specific locus
  is_rev <- read_direction == "reverse" & input_type == "padded" # unpadded reads are already turned forward by ObiTools, so every read is scored with the forward primers,
  #padded reads keep the the original orientation. this just checks with which primer should read be aligned based on direction and padding. UNpadded reads are not reversed even though the tag says so.
  
  asp1_variant_chr <- rep(as.character(asp1), n) #repeats it  n times, one copy per read
  asp1_variant_chr[is_rev] <- as.character(asp1_rc) #replaces the copies for reverse reads, based on where is_rev is true
  asp1_variant <- DNAStringSet(asp1_variant_chr) #turns it into pairwisealignemnt object
  
  asp2_variant_chr <- rep(as.character(asp2), n) #same for asp2
  asp2_variant_chr[is_rev] <- as.character(asp2_rc)
  asp2_variant <- DNAStringSet(asp2_variant_chr)
  
  # two alignment for asp1 and asp2 usong localglobal (primer can sit anywhere insdie read but it must be aligned from first to last base), each read uses own quality score, primers use preasigned 19Q. Also preasgined gap penalties
  if (probe_type == "reference") {
    # Reference probe: the two copies merged at the SNP into one IUPAC probe (e.g. C and T -> Y). fuzzyMatrix makes the IUPAC code a half match for both alleles.
    stopifnot("reference copies must have the SNP at the same position" = snp_pos1 == snp_pos2,
              "reference copies must differ only at the SNP" = substr(asp1_seq, 1, snp_pos1 - 1) == substr(asp2_seq, 1, snp_pos2 - 1))
    iupac_fuzzy <- nucleotideSubstitutionMatrix(match = 1, mismatch = 0, baseOnly = FALSE)
    probe <- DNAString(paste0(substr(asp1_seq, 1, snp_pos1 - 1),
                              mergeIUPACLetters(paste0(substr(asp1_seq, snp_pos1, snp_pos1), substr(asp2_seq, snp_pos2, snp_pos2))),
                              anchor))
    probe_rc <- reverseComplement(probe)
    # Which way each read runs: the read is scored against the probe and its reverse complement (score only, fast), the better one is its direction
    score_probe <- pairwiseAlignment(read_set, probe, type = "local-global", patternQuality = read_quality, subjectQuality = primer_quality_obj,
                                     gapOpening = gap_opening, gapExtension = gap_extension, fuzzyMatrix = iupac_fuzzy, scoreOnly = TRUE)
    score_probe_rc <- pairwiseAlignment(read_set, probe_rc, type = "local-global", patternQuality = read_quality, subjectQuality = primer_quality_obj,
                                        gapOpening = gap_opening, gapExtension = gap_extension, fuzzyMatrix = iupac_fuzzy, scoreOnly = TRUE)
    is_rev <- score_probe_rc > score_probe # replaces the direction label: from here on is_rev means "aligned to the reverse complemented probe"
    probe_variant_chr <- rep(as.character(probe), n)
    probe_variant_chr[is_rev] <- as.character(probe_rc)
    asp1_alignment <- pairwiseAlignment(read_set, DNAStringSet(probe_variant_chr), type = "local-global",
                                        patternQuality = read_quality, subjectQuality = primer_quality_obj,
                                        gapOpening = gap_opening, gapExtension = gap_extension, fuzzyMatrix = iupac_fuzzy)
    asp2_alignment <- asp1_alignment # one alignment stands in for both copies
  } else {
    asp1_alignment <- pairwiseAlignment(read_set, asp1_variant, type = "local-global",
                                        patternQuality = read_quality, subjectQuality = primer_quality_obj,
                                        gapOpening = gap_opening, gapExtension = gap_extension)
    asp2_alignment <- pairwiseAlignment(read_set, asp2_variant, type = "local-global",
                                        patternQuality = read_quality, subjectQuality = primer_quality_obj,
                                        gapOpening = gap_opening, gapExtension = gap_extension)
  }
  
  asp1_score <- score(asp1_alignment) #gets the scores, computes margin, winners index...
  asp2_score <- score(asp2_alignment)
  
  # Reference probe: the one score is turned into the two scores the two copies would give. They differ only at the SNP column:
  # the copy with the read's base gets a full match there instead of the half match, the other copy a mismatch. Third base or deleted SNP: both equal (tie).
  if (probe_type == "reference") {
    snp_ref <- map_snp_positions(as.character(pattern(asp1_alignment)), as.character(subject(asp1_alignment)),
                                 start(pattern(asp1_alignment)), ifelse(is_rev, nchar(anchor) + 1L, snp_pos1), snp_indel_window,
                                 primer_start = start(subject(asp1_alignment)))
    snp_quality_ref <- rep(NA_integer_, n)
    has_base <- which(!is.na(snp_ref$read_position))
    if (length(has_base) > 0) snp_quality_ref[has_base] <- utf8ToInt(paste(substring(read_quals[has_base], snp_ref$read_position[has_base], snp_ref$read_position[has_base]), collapse = "")) - 33L
    snp_base_probe_strand <- ifelse(is_rev, chartr("ACGT", "TGCA", snp_ref$read_base), snp_ref$read_base) # reads on the reverse complemented probe: complemented back
    is_allele1 <- snp_base_probe_strand %in% substr(asp1_seq, snp_pos1, snp_pos1)
    is_allele2 <- snp_base_probe_strand %in% substr(asp2_seq, snp_pos2, snp_pos2)
    half_match_lookup <- score(pairwiseAlignment(DNAStringSet(rep("C", 94)), DNAString("Y"), type = "local-global",
                                                 patternQuality = PhredQuality(vapply(0:93, function(q) intToUtf8(q + 33), "")),
                                                 subjectQuality = PhredQuality(intToUtf8(primer_quality + 33)), fuzzyMatrix = iupac_fuzzy)) # a base against an IUPAC code, read quality 0..93
    matching_score <- asp1_score + match_lookup[snp_quality_ref + 1] - half_match_lookup[snp_quality_ref + 1] # the copy with the read's base
    other_score <- matching_score - mismatch_cost_lookup[snp_quality_ref + 1] # the copy with the other allele
    asp1_score <- ifelse(is_allele1, matching_score, ifelse(is_allele2, other_score, asp1_score))
    asp2_score <- ifelse(is_allele2, matching_score, ifelse(is_allele1, other_score, asp2_score))
  }
  
  score_margin <- abs(asp1_score - asp2_score)
  winner_index <- ifelse(asp1_score >= asp2_score, 1L, 2L) #picks the best score, when it ties it picks the first one, but its later fixed (called allele NA)
  pick <- cbind(seq_len(n), winner_index) #creates a two column matrix, every read gets a row and a column with winner index
  
  winning_orientation <- paste0(c("asp1", "asp2")[winner_index], ifelse(is_rev, "_rc", "_fwd")) # unpadded: always _fwd, the real strand is in read_direction
  winning_nmismatch <- cbind(nmismatch(asp1_alignment), nmismatch(asp2_alignment))[pick] #takes the winners column for number of mismatches

  snp_column_to_skip <- NULL
  if (probe_type == "reference") {
    winning_nmismatch <- winning_nmismatch - (is_allele1 | is_allele2) # an allele at the IUPAC code is not a mismatch (nmismatch counts it), a third base stays one
    snp_column_to_skip <- ifelse(is_allele1 | is_allele2, ifelse(is_rev, nchar(anchor) + 1L, snp_pos1), 0L) # 0 matches no probe position
  }
  mismatch_summary1 <- mismatch_summary(asp1_alignment, n, mismatch_cost_lookup, skip_subject_position = snp_column_to_skip) #creates mimsatch summamry for asp1
  mismatch_summary2 <- mismatch_summary(asp2_alignment, n, mismatch_cost_lookup, skip_subject_position = snp_column_to_skip) #grabs all the info about mimsatches saved in the helper
  
  winning_mismatch_quality <- cbind(mismatch_summary1$mean_mismatch_quality, mismatch_summary2$mean_mismatch_quality)[pick]
  winning_mismatch_cost <- cbind(mismatch_summary1$mismatch_cost, mismatch_summary2$mismatch_cost)[pick]
  winning_max_mismatch_quality <- cbind(mismatch_summary1$max_mismatch_quality, mismatch_summary2$max_mismatch_quality)[pick]
  winning_high_quality_mismatches <- cbind(mismatch_summary1$high_quality_mismatches, mismatch_summary2$high_quality_mismatches)[pick]
  winning_indel_events <- cbind(mismatch_summary1$indel_events, mismatch_summary2$indel_events)[pick]
  
  asp1_pattern_start <- start(pattern(asp1_alignment)) #gives the read position where each alignments starts
  asp2_pattern_start <- start(pattern(asp2_alignment))
  
  #Computes exact perfect scores with helper function, the window is primer + anchor long (length(asp1) includes the anchor)
  perfect1 <- perfect_from_read_qualities(read_quals, asp1_pattern_start, length(asp1), match_lookup)
  perfect2 <- perfect_from_read_qualities(read_quals, asp2_pattern_start, length(asp2), match_lookup)
  
  ### SNP position
  # target is the position of the SNP on primer. For forward primers its the last ASP base (the anchor comes after it).
  # Reverse-complement primer: the anchor comes first, so the SNP sits at nchar(anchor) + 1.
  target1 <- ifelse(is_rev, nchar(anchor) + 1L, snp_pos1)
  target2 <- ifelse(is_rev, nchar(anchor) + 1L, snp_pos2)
  
  # locate the SNP in both alignemnts. it gives the function aligned read and primer with gaps "-", where the asp starts, which primer base is the snp adn a window
  snp1 <- map_snp_positions(as.character(pattern(asp1_alignment)), as.character(subject(asp1_alignment)),
                            asp1_pattern_start, target1, snp_indel_window,
                            primer_start = start(subject(asp1_alignment))) 
  snp2 <- map_snp_positions(as.character(pattern(asp2_alignment)), as.character(subject(asp2_alignment)),
                            asp2_pattern_start, target2, snp_indel_window,
                            primer_start = start(subject(asp2_alignment)))
  
  # picks the winners values and data
  snp_read_position <- cbind(snp1$read_position, snp2$read_position)[pick]
  snp_read_base <- cbind(snp1$read_base, snp2$read_base)[pick]
  snp_deleted <- cbind(snp1$deleted_in_read, snp2$deleted_in_read)[pick]
  snp_indel_flag <- cbind(snp1$indel_near_snp, snp2$indel_near_snp)[pick] # reported only
  
  # at least one primer found no SNP base in the read (deleted, or the read ends before it)
  snp_missing_any <- is.na(snp1$read_position) | is.na(snp2$read_position)
  # both primers found a SNP base, but at different read positions (homopolymers)
  snp_disagree_any <- !snp_missing_any & snp1$read_position != snp2$read_position
  
  snp_position_quality <- rep(NA_real_, n) #empty vector
  has_pos <- which(!is.na(snp_read_position)) # finds which reads have SNP base
  if (length(has_pos) > 0) {
    snp_q_chars <- substring(read_quals[has_pos], snp_read_position[has_pos], snp_read_position[has_pos]) #cuts out quality of the SNP position
    snp_position_quality[has_pos] <- utf8ToInt(paste(snp_q_chars, collapse = "")) - 33
  }
  
  ### end of SNP position
  
  if (threshold_type == "errors") { #uses error budget rule preasigned in settings
    error_cost <- mismatch_cost_lookup[21] # cost of one Q20 mismatch, typical error (lookup holds Phred 0..93, so Q20 is element 21)
    indel_cost <- gap_opening + gap_extension + max(match_lookup) # cost of one indel: the gap plus the best possible lost match
    min_score1 <- perfect1 - max(min_errors * indel_cost, error_rate * length(asp1) * error_cost) # length(asp1) includes the anchor. min score is perfect minus an allowed cost, which is the larger (max) of the two penalties or costs. 
    #Short primers get the floor usualy that is one indel or one good mismatch. longer primers get the per lentgh budget which is one Q20 mismatch per cca 12 bases (it wins usualy at primer longer than 18)
    min_score2 <- perfect2 - max(min_errors * indel_cost, error_rate * length(asp2) * error_cost)
  } else {
    min_score1 <- min_score_fraction * perfect1 #if the setting is fraction
    min_score2 <- min_score_fraction * perfect2
  }
  
  relative_margin <- score_margin / ((perfect1 + perfect2) / 2) #only reported
  
  ### CALLING AND FLAGS
  
  #grabs the better primer, if it passes its own score threshold (ties count as NA)
  called_allele_raw <- ifelse(asp1_score - asp2_score >= 1e-4 & asp1_score >= min_score1, allele1,
                              ifelse(asp2_score - asp1_score >= 1e-4 & asp2_score >= min_score2, allele2, NA))
  
  tie_flag <- abs(asp1_score - asp2_score) < 1e-4 # tolerance for floating-point noise
  below_threshold_flag <- !tie_flag & is.na(called_allele_raw) # winner below its score threshold
  has_raw_call <- !is.na(called_allele_raw) 
  
  winning_snp_base <- ifelse(winner_index == 1L, substr(asp1_seq, snp_pos1, snp_pos1), substr(asp2_seq, snp_pos2, snp_pos2))
  winning_snp_base[is_rev] <- chartr("ACGT", "TGCA", winning_snp_base[is_rev])
  
  snp_missing <- has_raw_call & snp_missing_any #flag for SNP base not found by one of the primers
  snp_position_disagree <- has_raw_call & snp_disagree_any # only when both found a base
  snp_mismatch_flag <- has_raw_call & !is.na(snp_read_base) & toupper(snp_read_base) != toupper(winning_snp_base)
  snp_low_quality <- has_raw_call & min_snp_quality > 0 & !is.na(snp_position_quality) & snp_position_quality < min_snp_quality
  too_many_mismatches <- has_raw_call & winning_nmismatch > max_mismatches
  
  # removal: a read votes only if it has a raw call and no flags!
  remove_read <- !has_raw_call | snp_mismatch_flag | snp_low_quality | too_many_mismatches |
    (require_position_agreement & (snp_missing | snp_position_disagree))
  
  called_allele <- ifelse(remove_read, NA, called_allele_raw) # final call
  
  ### END OF CALLING AND FLAGS
  
  #this is the output. it saves all the scores for winning primer and is later saved into all_reads.tsv
  list(score1 = asp1_score, #alignments scores
       score2 = asp2_score,
       perfect1 = perfect1,#best scores possible
       perfect2 = perfect2,
       called_allele = called_allele, #final called allele (NA if any flag is set)
       called_allele_raw = called_allele_raw, # the score-based call before any flag (NA only for ties and below threshold)
       winning_orientation = winning_orientation,
       winning_nmismatch = winning_nmismatch, #number of mismatches, checked by cap
       winning_mismatch_quality = winning_mismatch_quality, #mean Phred of mismatches
       winning_mismatch_cost = winning_mismatch_cost, #score lost to mismatches
       winning_max_mismatch_quality = winning_max_mismatch_quality, #highest quality mismatch
       winning_high_quality_mismatches = winning_high_quality_mismatches, #number of mismatches over Q20
       winning_indel_events = winning_indel_events, #number of indels
       snp_position_quality = snp_position_quality, #info about SNPs
       snp_deleted = snp_deleted,
       snp_read_position = snp_read_position,
       snp_read_base = snp_read_base,
       tie_flag = tie_flag, # flag, removes the vote
       below_threshold_flag = below_threshold_flag, # flag, removes the vote
       snp_mismatch_flag = snp_mismatch_flag, #flag, removes the vote  (for raw calls only)
       snp_indel_flag = snp_indel_flag, #report only
       snp_missing = snp_missing, #flag, SNP base not found by one primer (removes the vote with require_position_agreement)
       snp_position_disagree = snp_position_disagree, # both primers found a SNP base but at different positions, removes the vote
       snp_low_quality = snp_low_quality, #flag, removes the vote (now quality-based)
       too_many_mismatches = too_many_mismatches, #flag, removes the vote 
       score_margin = score_margin,
       relative_margin = relative_margin)
}


### Genotype caller per sample and locus pair from the reads' called alleles, with the strand check for heterozygotes
call_genotype <- function(alleles, #vector of every read's called allele per locus and sample pair (includes NA if the read didnt pass a threshold)
                          directions, # read direction of every read
                          minimum_required_depth, # number of usable reads needed for calling a genotype, editable
                          HT_min_frequence, # minimal fraction of the second allele for it to call a HT, editable
                          strand_flag_prob = 0.01, # accepted chance of wrongly flagging a real HT per direction, minimum reads needed to flag
                          expected_minor_fraction = 0.44, # expected share of the second allele in a real HT 
                          strand_downgrade = FALSE) { #TRUE = turns HT into HM, FALSE = only flagged
  
  stopifnot(length(alleles) == length(directions)) # every read needs its direction
  min_other_strand <- ceiling(log(strand_flag_prob) / log(1 - expected_minor_fraction)) # number of reads a direction needs before its flagged
  
  keep <- !is.na(alleles) 
  alleles <- alleles[keep] # strips out NAs, reads that didnt produce confident calls
  directions <- directions[keep] # directions stay paired with their reads
  n <- length(alleles) # number of actual usable alleles
  
  out <- c(genotype = NA_character_, depth = as.character(n), # empty result vector, always text so every pair returns the same type (bind_rows needs that)
           strand_check = NA_character_, minor_fwd = NA_character_, minor_rev = NA_character_)
  if (n < minimum_required_depth)
    return(out) # if it doesnt pass the required min depth its NA, and saves the depth
  
  tab <- sort(table(alleles), decreasing = TRUE) # table counts each allele and sorts them
  top_allele <- names(tab)[1] #pulls out the name of top allele
  second_fraction <- if (length(tab) >= 2) { tab[[2]] / n } else { 0 } #computes a fraction of the second allele
  
  # actual genotyping decision! If every read agreed or second fraction is below the min_frequence the locus is called HM
  out["genotype"] <- paste0(top_allele, top_allele)
  if (length(tab) == 1 || second_fraction < HT_min_frequence) {
    return(out)
  }
  
  # strand check for HT. Count the second allele and all usable reads per direction.
  second_allele <- names(tab)[2]
  is_fwd <- directions %in% "forward" 
  is_rev <- directions %in% "reverse"
  minor_fwd <- sum(alleles == second_allele & is_fwd)
  minor_rev <- sum(alleles == second_allele & is_rev)
  n_fwd <- sum(is_fwd)
  n_rev <- sum(is_rev)
  out["minor_fwd"] <- minor_fwd
  out["minor_rev"] <- minor_rev
  
  allelescore2 <- sort(c(top_allele, second_allele)) # calls using both alleles and alphabetically sorted for consistency
  missing_fwd <- n_fwd >= min_other_strand && minor_fwd == 0 #enough forward reads, but none for second allele
  missing_rev <- n_rev >= min_other_strand && minor_rev == 0 #same for reverse
  out["strand_check"] <- if (missing_fwd || missing_rev) "failed" else "passed" # second allele missing from one direction or seen in both directions 
  
  if (out["strand_check"] == "failed" && strand_downgrade) return(out) # with strand_downgrade = TRUE the genotype stays HM for the top allele
  out["genotype"] <- paste0(allelescore2[1], allelescore2[2]) #flag only, HT is kept
  out
}
