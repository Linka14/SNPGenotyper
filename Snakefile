import glob
import yaml
import subprocess

configfile: "config.yaml"

# Settings each stage actually uses. Each rmd stage has settings in the yaml header 
# which each one uses specifically, those are taken from the config.yaml sections. Otherwise changing the settings from Bayes stage
# would rerun whole pipeline. It reads the rmd stage, splits it by "---", takes just the yaml header, turns it into dictionary
# and grabs just the params block. Then for every setting name inside Rmd it finds the config.yaml value and and updates the settings.

def stage_settings(stage_rmd, section_names):
    declared_settings = yaml.safe_load(open(stage_rmd).read().split("---")[1])["params"]
    combined_settings = {}
    for section_name in section_names:
        combined_settings.update(config[section_name])
    return {name: combined_settings[name] for name in declared_settings if name in combined_settings}

# Settings or values used to build the file paths are taken from config.yaml 
dataset_name = config["shared_settings"]["dataset_name"]
barcode_numbers = [f"{barcode_number:02d}" for barcode_number in range(1, config["shared_settings"]["n_barcodes"] + 1)]
demultiplexed_folder = f"Data/Processed_data/{dataset_name}/Demultiplexed_sequence_data"
raw_reads = glob.glob(f"Data/Raw_data/{dataset_name}/*.fastq") or glob.glob(f"Data/Raw_data/{dataset_name}/*.bam")
dorado_fastqs = [f"{demultiplexed_folder}/Dorado_demultiplexed_by_sample/barcode{barcode}.fastq" for barcode in barcode_numbers]
identified_fastqs = [f"{demultiplexed_folder}/Identified_sample_sequences_Obi4/barcode{barcode}_identified.fastq" for barcode in barcode_numbers]
unidentified_fastqs = [f"{demultiplexed_folder}/Unidentified_sample_sequences_Obi4/barcode{barcode}_unidentified.fastq" for barcode in barcode_numbers]
padded_fastqs = [f"{demultiplexed_folder}/Identified_padded_sample_sequences/barcode{barcode}_padded.fastq" for barcode in barcode_numbers]
minimap2_folder = f"Results/{dataset_name}/Genotypes_results/MiniMap2_genotypes"
locus_summary = f"Results/{dataset_name}/Graphics/Demultiplexing_and_Filtration_Visualization/locus_indel_homopolymer_summary.tsv"
### Stage 04 prep: reads different input fastq files based on input_type setting, those are either "padded" or "unpadded" reads.
asp_input_fastqs = padded_fastqs if config["shared_settings"]["input_type"] == "padded" else identified_fastqs


### Stage 05 prep: same folder name building idea as with the 04 Stage, based on settings uploaded function from .R script
run_names = subprocess.run(["Rscript", "-e",
    'source("Code/R_scripts/asp_functions.R"); config <- yaml::read_yaml("config.yaml"); s <- c(config$shared_settings, config$settings_05_bayesian); '
    'run_names <- build_run_names(s$input_type, s$use_anchors, s$anchor_length, s$threshold_type, s$error_rate, s$min_score_fraction, '
    's$prior_type, s$prior_het_share, s$prior_population, s$prior_floor, probe_type = s$probe_type, '
    'ancestry_prior_shape = s$ancestry_prior_shape, phi_hom_structure = s$phi_hom_structure); '
    'cat(run_names$asp_run_name, run_names$bayes_run_name, run_names$bayes_settings_label, sep = "\\n")'],
    capture_output = True, text = True, check = True).stdout
    
# Splits the names to each variable    
asp_run_name, bayes_run_name, bayes_settings_label = run_names.strip().split("\n")
asp_folder = f"Results/{dataset_name}/Genotypes_results/{asp_run_name}"
bayes_folder = f"Results/{dataset_name}/Genotypes_results/{bayes_run_name}"

# Downsampling folder, built from the downsampling settings the same way as in the 05 setup chunk
settings_05 = config["settings_05_bayesian"]
downsampling_folder = (f"Results/{dataset_name}/Graphics/STAN_Hierarchical/{bayes_settings_label}/Downsampling_calibration"
                       f"{'_sdown' if config['shared_settings']['strand_downgrade'] else '_sflag'}"
                       f"_r{settings_05['minimum_reads_for_truth_pair']}_d{max(settings_05['read_depths_to_test'])}_n{settings_05['draws_per_depth']}")
                      
rule all:
    input:
        "Code/HTMLs/02_ObiTools4_Basic_Visualisation.html",
        "Code/HTMLs/03_SNP_Genotyping_MiniMap2.html",
        "Code/HTMLs/06b_Genotyping_Results.html"

rule stage_01a_dorado: 
    input:
      raw_reads = raw_reads
    output:
      html = "Code/HTMLs/01a_Dorado_Demultiplexing.html",
      dorado = dorado_fastqs
    params:
      settings = stage_settings("Code/Rmd_Stages/01a_Dorado_Demultiplexing.Rmd", ["shared_settings", "settings_01_demultiplexing"]),
      timer = config["workflow_settings"]["timer_command"]
    log:
      "logs/stage_01a.log"
    benchmark:
      "benchmarks/stage_01a.tsv"
    shell:
      "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/01a_Dorado_Demultiplexing.Rmd config.yaml shared_settings settings_01_demultiplexing 2>&1 | tee {log}"
      
rule stage_01b_obitools:
    input:
      dorado = dorado_fastqs,
      ngsfilter = config["settings_01_demultiplexing"]["ngsfilter_path"]
    output:
      html = "Code/HTMLs/01b_ObiTools4_Demultiplexing.html",
      identified = identified_fastqs,
      unidentified = unidentified_fastqs
    params:
      settings = stage_settings("Code/Rmd_Stages/01b_ObiTools4_Demultiplexing.Rmd", ["shared_settings", "settings_01_demultiplexing"]),
      timer = config["workflow_settings"]["timer_command"]
    log:
      "logs/stage_01b.log"
    benchmark:
      "benchmarks/stage_01b.tsv"
    shell:
      "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/01b_ObiTools4_Demultiplexing.Rmd config.yaml shared_settings settings_01_demultiplexing 2>&1 | tee {log}"  

rule stage_02_obitools_summaries:
    input: 
      identified = identified_fastqs,
      unidentified = unidentified_fastqs,
      barcode_map = config["shared_settings"]["barcode_map_path"]
    output:
      "Code/HTMLs/02_ObiTools4_Basic_Visualisation.html"
    params:
      settings = stage_settings("Code/Rmd_Stages/02_ObiTools4_Basic_Visualisation.Rmd",["shared_settings"]),
      timer = config["workflow_settings"]["timer_command"]
    log:
      "logs/stage_02.log"
    benchmark:
     "benchmarks/stage_02.tsv"
    shell:
      "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/02_ObiTools4_Basic_Visualisation.Rmd config.yaml shared_settings 2>&1 | tee {log}"

rule stage_03_minimap2:
  input:
    dorado = dorado_fastqs,
    identified = identified_fastqs,
    reference = config["shared_settings"]["reference_path"],
    primer_panel = config["shared_settings"]["primer_panel_path"],
    barcode_map = config["shared_settings"]["barcode_map_path"]
  output:
    html = "Code/HTMLs/03_SNP_Genotyping_MiniMap2.html",
    padded = padded_fastqs,
    genotypes = f"{minimap2_folder}/all_barcodes_Minimap2_genotypes.tsv",
    mapping = f"{minimap2_folder}/mapping_summary_03.tsv"
  params:
    settings = stage_settings("Code/Rmd_Stages/03_SNP_Genotyping_MiniMap2.Rmd", ["shared_settings", "settings_03_minimap2"]),
    timer = config["workflow_settings"]["timer_command"]
  log:
    "logs/stage_03.log"
  benchmark:
     "benchmarks/stage_03.tsv"
  shell: 
    "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/03_SNP_Genotyping_MiniMap2.Rmd config.yaml shared_settings settings_03_minimap2 2>&1 | tee {log}"
    
rule stage_04_asp_genotyping:
  input:
    reads = asp_input_fastqs,
    functions = "Code/R_scripts/asp_functions.R",
    reference = config["shared_settings"]["reference_path"],
    primer_panel = config["shared_settings"]["primer_panel_path"],
    barcode_map = config["shared_settings"]["barcode_map_path"]
  output:
    html = "Code/HTMLs/04_SNP_Genotyping_ASP_primer.html",
    genotypes = f"{asp_folder}/all_barcodes_asp_genotypes.tsv",
    reads = f"{asp_folder}/all_barcodes_asp_reads.tsv"
  params:
    settings = stage_settings("Code/Rmd_Stages/04_SNP_Genotyping_ASP_primer.Rmd",["shared_settings", "settings_04_genotyping"]),
    timer = config["workflow_settings"]["timer_command"]
  threads:
    config["shared_settings"]["n_cores"]
  log:
    "logs/stage_04.log"
  benchmark:
    "benchmarks/stage_04.tsv"
  shell:
    "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/04_SNP_Genotyping_ASP_primer.Rmd config.yaml shared_settings settings_04_genotyping 2>&1 | tee {log}"

rule stage_05a_bayesian:
  input:
    asp_reads = f"{asp_folder}/all_barcodes_asp_reads.tsv",
    asp_genotypes = f"{asp_folder}/all_barcodes_asp_genotypes.tsv",
    stan_model = "Stan_config/base_mixture_model_hierarchical.stan",
    primer_panel = config["shared_settings"]["primer_panel_path"]
  output:
    html = "Code/HTMLs/05a_SNP_Bayesian_Genotyping.html",
    genotypes = f"{bayes_folder}/all_barcodes_base_mixture_genotypes_hierarchical.tsv",
    samples = f"{bayes_folder}/sample_parameters_hierarchical.tsv",
    stan_results = f"{bayes_folder}/stan_plot_inputs.rds"
  params:
    settings = stage_settings("Code/Rmd_Stages/05a_SNP_Bayesian_Genotyping.Rmd", ["shared_settings", "settings_05_bayesian"]),
    timer = config["workflow_settings"]["timer_command"]
  threads:
    config["shared_settings"]["n_cores"]
  log:
    "logs/stage_05a.log"
  benchmark:
    "benchmarks/stage_05a.tsv"
  shell:
    "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/05a_SNP_Bayesian_Genotyping.Rmd config.yaml shared_settings settings_05_bayesian 2>&1 | tee {log}"

rule stage_05b_downsampling:
  input:
    stan_results = f"{bayes_folder}/stan_plot_inputs.rds",
    asp_reads = f"{asp_folder}/all_barcodes_asp_reads.tsv",
    reference_genotypes = f"Data/Supporting_data/Reference_genotypes/genotypes_{dataset_name}.tsv",
    functions = "Code/R_scripts/asp_functions.R"
  output:
    html = "Code/HTMLs/05b_Downsampling_Calibration.html",
    downsampling_results = f"{downsampling_folder}/downsampling_plot_inputs.rds"
  params:
    settings = stage_settings("Code/Rmd_Stages/05b_Downsampling_Calibration.Rmd", ["shared_settings", "settings_05_bayesian"]),
    timer = config["workflow_settings"]["timer_command"]
  threads:
    config["shared_settings"]["n_cores"]
  log:
    "logs/stage_05b.log"
  benchmark:
    "benchmarks/stage_05b.tsv"
  shell:
    "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/05b_Downsampling_Calibration.Rmd config.yaml shared_settings settings_05_bayesian 2>&1 | tee {log}"

rule stage_06a_diagnostics:
    input:
        asp_reads = f"{asp_folder}/all_barcodes_asp_reads.tsv",
        asp_genotypes = f"{asp_folder}/all_barcodes_asp_genotypes.tsv",
        bayes_genotypes = f"{bayes_folder}/all_barcodes_base_mixture_genotypes_hierarchical.tsv",
        reference_genotypes = f"Data/Supporting_data/Reference_genotypes/genotypes_{dataset_name}.tsv",
        raw_reads = raw_reads,
        dorado = dorado_fastqs,
        identified = identified_fastqs,
        unidentified = unidentified_fastqs,
        padded = padded_fastqs,
        reference = config["shared_settings"]["reference_path"],
        primer_panel = config["shared_settings"]["primer_panel_path"],
        functions = "Code/R_scripts/asp_functions.R"
    output:
        html = "Code/HTMLs/06a_Pipeline_Diagnostics.html",
        locus_summary = locus_summary
    params:
        settings = stage_settings("Code/Rmd_Stages/06a_Pipeline_Diagnostics.Rmd", ["shared_settings", "settings_04_genotyping", "settings_05_bayesian"]),
        timer = config["workflow_settings"]["timer_command"]
    threads:
        config["shared_settings"]["n_cores"]
    log:
        "logs/stage_06a.log"
    benchmark:
        "benchmarks/stage_06a.tsv"
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/06a_Pipeline_Diagnostics.Rmd config.yaml shared_settings settings_04_genotyping settings_05_bayesian 2>&1 | tee {log}"

rule stage_06b_results:
    input:
        asp_genotypes = f"{asp_folder}/all_barcodes_asp_genotypes.tsv",
        bayes_genotypes = f"{bayes_folder}/all_barcodes_base_mixture_genotypes_hierarchical.tsv",
        minimap2_genotypes = f"{minimap2_folder}/all_barcodes_Minimap2_genotypes.tsv",
        reference_genotypes = f"Data/Supporting_data/Reference_genotypes/genotypes_{dataset_name}.tsv",
        stan_results = f"{bayes_folder}/stan_plot_inputs.rds",
        downsampling_results = f"{downsampling_folder}/downsampling_plot_inputs.rds",
        locus_summary = locus_summary,
        functions = "Code/R_scripts/asp_functions.R",
        primer_panel = config["shared_settings"]["primer_panel_path"],
        reference = config["shared_settings"]["reference_path"]
    output:
        "Code/HTMLs/06b_Genotyping_Results.html"
    params:
        settings = stage_settings("Code/Rmd_Stages/06b_Genotyping_Results.Rmd", ["shared_settings", "settings_05_bayesian"]),
        timer = config["workflow_settings"]["timer_command"]
    log:
        "logs/stage_06b.log"
    benchmark:
        "benchmarks/stage_06b.tsv"
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/06b_Genotyping_Results.Rmd config.yaml shared_settings settings_05_bayesian 2>&1 | tee {log}"

rule stage_06c_speed:
    input:
        asp_reads = f"{asp_folder}/all_barcodes_asp_reads.tsv",
        reference_genotypes = f"Data/Supporting_data/Reference_genotypes/genotypes_{dataset_name}.tsv",
        functions = "Code/R_scripts/asp_functions.R"
    output:
        "Code/HTMLs/06c_Speed_Optimisation.html"
    params:
        settings = stage_settings("Code/Rmd_Stages/06c_Speed_Optimisation.Rmd", ["shared_settings"]),
        timer = config["workflow_settings"]["timer_command"]
    log:
        "logs/stage_06c.log"
    benchmark:
        "benchmarks/stage_06c.tsv"
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/06c_Speed_Optimisation.Rmd config.yaml shared_settings 2>&1 | tee {log}"

