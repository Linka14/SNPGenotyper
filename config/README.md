# Configuration

All the settings form the pipeline are inside "config/config.yaml". Snakefile reads it, decides what needs to be rerun based on changed
settings, and build appropriate folder naemes based on settings, so different run get saved inside separate folders.
R script saved insid `Code/R_scripts/render_stage.R` passes those settings to each .Rmd stage as params, each stage only gets the
settings which are declared inside the YAML header, changing those only reruns stages whose setting inside YAML has been changed.

Before everything runs, Snakemake checks whether the settings chosen inside config.yaml match those declared and saved inside
`config/config.schema.yaml`. Any misspelled string, weirds value, out of bounds value stops the run with a message naming that specific 
setting. 

All paths are written from te project folder (Full_workflow)

## SHARED SETTINGS (used by several stages)

| Setting -> [] -> Meaning

| `dataset_name` -> ["big", "small"] -> Any other dataset name (folder names)
| `n_barcodes` -> [24] -> Any number of barcodes (6 for small, 24 for big)
| `dorado_bin_folder`, `obitools_bin_folder`, `minimap2_bin_folder`, `samtools_bin_folder` -> [""] -> Folders with the installed programs, empty means they are found automatically
| `barcode_map_path` -> ["Data/Supporting_data/Samples_ONT_barcodes.tsv"] -> Sample to barcode table
| `primer_panel_path` -> "["Data/Supporting_data/Primer_panels/...xlsx"] -> Primer panel with alleles and allele frequencies
| `reference_path` -> [Data/Supporting_data/Genome_reference_data/snps_500bp_sequences.fa] -> 1000 bp SNP reference sequences
| `minimum_required_depth` -> [3] -> Minimum depth per sample locus pair to call a genotype
| `HT_min_frequence` -> [0.20] -> MAF needed to call a HT
| `min_snp_quality` -> [10] -> Reads whose SNP base has a lower Phred quality do not vote (ASP) or are not counted (Minimap2), 0 switches it off
| `n_cores` -> [4] -> Number of cores for ASP scoring, Stan chains and downsampling
| `input_type` -> ["unpadded", "padded"] -> Reads for ASP genotyping
| `use_anchors` -> [true] -> Reference anchor bases after the SNP in the ASP primers
| `anchor_length` -> [3] -> 1 to 5 anchor bases
| `probe_type` -> ["asp_primers", "asp_primers" or "reference"] ->  Needs `input_type: "padded"`
| `threshold_type` -> ["errors" or "fraction"] -> Alignment threshold
| `min_score_fraction` -> [0.65] -> For "fraction" setting only, share of the perfect primer alignment score
| `error_rate` -> [0.08] -> For "errors" setting only, allowed errors per primer base
| `strand_flag_prob`, `expected_minor_fraction` -> [0.01, 0.44] -> Strand bias check of heterozgyotes
| `strand_downgrade` -> ["true", "false"] -> "true" turns flagged heterozygotes into homozygotes, "false"" only flags them

# ASP settings `input_type` to `error_rate` make up the ASP results folder name, so every combination gets its own folder.


## settings_01_demultiplexing

| `kit_name` -> ["SQK-NBD114-24""] -> Nanopore barcoding kit for dorado demultiplexing
| `ngsfilter_path` -> ["Data/Supporting_data/NGSfilter_file/ngsfilter.csv"] -> NGSfilter file for obimultiplex (prepared by hand)



## settings_03_minimap2

| `pad` -> [50] -> bp of extra context restored on each side of the reads
| `snp_position` -> [500] -> SNP position within every 1000 bp reference
| `minimum_mapping_quality_threshold` -> [20] -> Minimum mapping quality for a read to count



## settings_04_genotyping

| `primer_quality` -> [19] -> Phred quality given to every primer base in the alignment
| `gap_opening` -> [6] -> Gap cost of opening
| `gap_extension` -> [2] -> Gap cost of extension
| `min_errors` -> [1] -> For "errors" setting only, minimum allowed errors per primer
| `max_mismatches` -> [2] -> reads with more primer mismatches do not vote, ".inf"" switches it off



## settings_05_bayesian (The prior settings make up the Bayes results folder name.)

| `prior_type` -> ["flat", "het_share", "hwe" or "ancestry"] -> Genotype prior of the main Bayes run
| `prior_het_share` -> [0.15] -> for "het_share" setting only
| `prior_population` -> ["FreqWOLF_EastEurope"] -> Allele frequency column used by "hwe" and "ancestry"
| `prior_floor` -> [0.001] -> Lower limit of allele frequencies, no genotype gets a zero prior
| `ancestry_prior_shape` -> [1] -> Beta prior of the wolf ancestry, 1 is uniform
| `phi_hom_structure` -> ["shared", "locus", "sample" or "locus_sample"] -> HM overdispersion
| `stan_iter` -> [4000] -> Stan sampling settings
| `stan_warmup` -> [1500] -> Number of removed steps
| `stan_chains` -> [4] -> Number of chains
| `stan_adapt_delta` -> [0.95] -> Stans sampler step tuning parameter
| `max_error` -> [0.1] -> Cap on any error rate in the model
| `read_depths_to_test` -> [3, 4, 5, 6, 7, 8, 10, 12, 15, 20, 25, 30, 40, 50, 60, 75] -> Read depths to test fo nrdownsampling
| `draws_per_depth` -> [300] -> Number of draws per downsampling calibration
| `minimum_reads_for_truth_pair` -> [75] -> Depth cut-off for pairs which are tested
| `bootstrap_rounds` -> [1000] -> NUmber of bootstraping rounds



## settings_07_failed_genotypes

| `failed_genotypes_source` -> ["asp_genotypes" or "bayes_genotypes"] -> Which calls are compared with the reference genotypes
| `failed_max_reads_plotted` -> [20] -> Reads drawn per amplicon plot, .tsv keeps every read though



## settings_08_ancestry_validation

| `genotype_prior_type` -> ["flat"] -> prior of the second Bayes run whose genotypes go into PCA, DAPC, NewHybrids and STRUCTURE. Must differ from `prior_type` (normally "flat", to avoid circularity)
| `gq_cutoff` -> [20] -> Genotypes with a lower GQ are set to missing
| `reference_source` -> ["real" or "simulated"] -> Reference populations from the .ped/.map files or simulated from frequencies
| `reference_ped_path`  -> ["Data/Supporting_data/Reference_populations/GTIC-D-20-00312R1.ped"] -> Real reference populations
| `reference_map_path` -> ["Data/Supporting_data/Reference_populations/GTIC-D-20-00312R1.ped"] -> Real reference populations
| `reference_wolf_groups` -> ["wolGer", "wolRom", "wolFin", "wolRus"] -> EastEuropean sub populations
| `reference_dog_groups` -> ["dog_Illumina", "dog_nonpedigree"] -> Popultions used as reference for dogs
| `wolf_frequency_column` -> ["FreqWOLF_EastEurope"] -> Panel frequency columns for simulated references: FreqWOLF_all, FreqWOLF_EastEurope, FreqWOLF_IT or FreqWOLF_IB
| `dog_frequency_column` -> ["FreqDOG"] -> Panel frequency columns for simulated references for dogs
| `n_reference_wolves`, `n_reference_dogs` -> [186 or 300] -> Reference animals per group (at least 100, max 186 for the real East European wolves)
| `n_simulated_per_hybrid_class` -> [186 or 300] -> Simulated F1 and backcrosses per class (at least 100)
| `simulation_seed` -> [1] -> Random seed
| `structure_runs` -> [5] -> Independent STRUCTURE runs



## pipeline_settings 

| `timer_command` -> ["/usr/bin/time -l"or "/usr/bin/time -v"] -> for macOS or Linux, "" switches it off
| `newhybrids_command` -> path to the NewHybrids program
| `structure_command` -> path to the STRUCTURE program

Run options of Snakemake itself (for example the limit of one NewHybrids or STRUCTURE job at a time) are in `profiles/default/config.yaml`.

