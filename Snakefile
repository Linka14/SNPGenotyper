import glob
import yaml
from snakemake.utils import validate
import subprocess


configfile: "config/config.yaml"
validate(config, "config/config.schema.yaml")

### ==================================
# Settings each stage actually uses. Each rmd stage has settings in the yaml header
# which each one uses specifically, those are taken from the config.yaml sections. Otherwise changing the settings from Bayes stage
# would rerun whole pipeline. It reads the rmd stage, splits it by "---", takes just the yaml header, turns it into dictionary
# and grabs just the params block. Then for every setting name inside Rmd it finds the config.yaml value and updates the settings.

def stage_settings(stage_rmd, section_names):
    declared_settings = yaml.safe_load(open(stage_rmd).read().split("---")[1])["params"]
    combined_settings = {}
    for section_name in section_names:
        combined_settings.update(config[section_name])
    return {
        name: combined_settings[name]
        for name in declared_settings
        if name in combined_settings
    }

### ==================================
### Settings or values used to build the file paths are taken from config.yaml

dataset_name = config["shared_settings"]["dataset_name"]
barcode_numbers = [
    f"{barcode_number:02d}"
    for barcode_number in range(1, config["shared_settings"]["n_barcodes"] + 1)
]
settings_05 = config["settings_05_bayesian"]
settings_08 = config["settings_08_ancestry_validation"]
genotype_prior_type = settings_08["genotype_prior_type"]
if genotype_prior_type == settings_05["prior_type"]:
    raise ValueError(
        "genotype_prior_type (settings_08) must differ from prior_type (settings_05), otherwise both 05a rules make the same files"
    )

raw_reads = glob.glob(f"Data/Raw_data/{dataset_name}/*.fastq") or glob.glob(
    f"Data/Raw_data/{dataset_name}/*.bam"
)
reference_genotypes = (
    f"Data/Supporting_data/Reference_genotypes/genotypes_{dataset_name}.tsv"
)
asp_functions = "Code/R_scripts/asp_functions.R"

### ==================================
### Folder name for every run, its built in R with build_runnames(), uses shared function from asp_funtcions.R script

run_names = subprocess.run(
    [
        "Rscript",
        "-e",
   'source("Code/R_scripts/asp_functions.R"); config <- yaml::read_yaml("config/config.yaml"); s <- c(config$shared_settings, config$settings_05_bayesian); '
        "run_names <- build_run_names(s$input_type, s$use_anchors, s$anchor_length, s$threshold_type, s$error_rate, s$min_score_fraction, "
        "s$prior_type, s$prior_het_share, s$prior_population, s$prior_floor, probe_type = s$probe_type, "
        "ancestry_prior_shape = s$ancestry_prior_shape, phi_hom_structure = s$phi_hom_structure); "
        "genotype_run_names <- build_run_names(s$input_type, s$use_anchors, s$anchor_length, s$threshold_type, s$error_rate, s$min_score_fraction, "
        "config$settings_08_ancestry_validation$genotype_prior_type, s$prior_het_share, s$prior_population, s$prior_floor, probe_type = s$probe_type, "
        "ancestry_prior_shape = s$ancestry_prior_shape, phi_hom_structure = s$phi_hom_structure); "
        "cat(run_names$asp_run_name, run_names$bayes_run_name, run_names$bayes_settings_label, "
        'genotype_run_names$bayes_run_name, genotype_run_names$bayes_settings_label, sep = "\\n")',
    ],
    capture_output=True,
    text=True,
    check=True,
).stdout

# Splits the names to each variable
(
    asp_run_name,
    bayes_run_name,
    bayes_settings_label,
    genotype_bayes_run_name,
    genotype_settings_label,
) = run_names.strip().split("\n")

### ==================================
### Files for stages 01 to 03, demultiplexing reads, obitools summaries and Minimap2

demultiplexed_folder = f"Data/Processed_data/{dataset_name}/Demultiplexed_sequence_data"
dorado_fastqs = [
    f"{demultiplexed_folder}/Dorado_demultiplexed_by_sample/barcode{barcode}.fastq"
    for barcode in barcode_numbers
]
identified_fastqs = [
    f"{demultiplexed_folder}/Identified_sample_sequences_Obi4/barcode{barcode}_identified.fastq"
    for barcode in barcode_numbers
]
unidentified_fastqs = [
    f"{demultiplexed_folder}/Unidentified_sample_sequences_Obi4/barcode{barcode}_unidentified.fastq"
    for barcode in barcode_numbers
]
padded_fastqs = [
    f"{demultiplexed_folder}/Identified_padded_sample_sequences/barcode{barcode}_padded.fastq"
    for barcode in barcode_numbers
]
minimap2_folder = f"Results/{dataset_name}/Genotypes_results/MiniMap2_genotypes"

### ==================================
### Stage 04, ASP genotyping
# Reads prep from different input fastq files based on input_type setting, those are either "padded" or "unpadded" reads.

asp_input_fastqs = (
    padded_fastqs
    if config["shared_settings"]["input_type"] == "padded"
    else identified_fastqs
)

asp_folder = f"Results/{dataset_name}/Genotypes_results/{asp_run_name}"
asp_genotypes_table = f"{asp_folder}/all_barcodes_asp_genotypes.tsv"
asp_reads_table = f"{asp_folder}/all_barcodes_asp_reads.tsv"

### ==================================
### Stage 05 for Bayesian genotyping, ancestry and flat prior runs, downsampling

bayes_folder = f"Results/{dataset_name}/Genotypes_results/{bayes_run_name}"
bayes_genotypes_table = (
    f"{bayes_folder}/all_barcodes_base_mixture_genotypes_hierarchical.tsv"
)
genotype_bayes_folder = (
    f"Results/{dataset_name}/Genotypes_results/{genotype_bayes_run_name}"
)

# Downsampling folder, built from the downsampling settings the same way as in the 05 setup chunk
downsampling_folder = (
    f"Results/{dataset_name}/Graphics/STAN_Hierarchical/{bayes_settings_label}/Downsampling_calibration"
    f"{'_sdown' if config['shared_settings']['strand_downgrade'] else '_sflag'}"
    f"_r{settings_05['minimum_reads_for_truth_pair']}_d{max(settings_05['read_depths_to_test'])}_n{settings_05['draws_per_depth']}"
)

### ==================================
### Stage 06 for diagnostics and results

locus_summary = f"Results/{dataset_name}/Graphics/Demultiplexing_and_Filtration_Visualization/locus_indel_homopolymer_summary.tsv"

### ==================================
### Stage 07 for failed genotypes analysis and plots
# Checking which setting (genotypes table from which specific method to check) is set in config.yaml, based on that it chooses the
# input and output folder

failed_from_asp = (
    config["settings_07_failed_genotypes"]["failed_genotypes_source"] == "asp_genotypes"
)
failed_genotypes_table = (
    asp_genotypes_table if failed_from_asp else bayes_genotypes_table
)
failed_folder = f"Results/{dataset_name}/Graphics/Failed_Genotypes_Analysis/{asp_run_name if failed_from_asp else bayes_run_name}"

### ==================================
### Stage 08 for ancestry analysis, NewHybrids and STRUCTURE
# the genotype run (flat prior) used by the ancestry methods, and the folders of 08a, NewHybrids and STRUCTURE (same as ancestry_out_dir in the 08 setup chunks)
ancestry_folder = f"Results/{dataset_name}/Graphics/STAN_Hierarchical/{genotype_settings_label}/Ancestry_validation_{settings_08['reference_source']}"
newhybrids_folder = f"{ancestry_folder}/NewHybrids"
newhybrids_power_folder = f"{ancestry_folder}/NewHybrids_power"
structure_folder = f"{ancestry_folder}/STRUCTURE"
structure_outputs = [
    f"{structure_folder}/structure_K2_run{run}_f"
    for run in range(1, settings_08["structure_runs"] + 1)
]

### ==================================


rule all:
    input:
        "Code/HTMLs/02_ObiTools4_Basic_Visualisation.html",
        "Code/HTMLs/03_SNP_Genotyping_MiniMap2.html",
        "Code/HTMLs/06b_Genotyping_Results.html",
        "Code/HTMLs/07_Failed_Genotypes_Analysis.html",
        "Code/HTMLs/08b_Ancestry_Results.html",


rule stage_01a_dorado:
    input:
        raw_reads=raw_reads,
    output:
        html="Code/HTMLs/01a_Dorado_Demultiplexing.html",
        dorado=dorado_fastqs,
    log:
        "logs/stage_01a.log",
    benchmark:
        "benchmarks/stage_01a.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/01a_Dorado_Demultiplexing.Rmd",
            ["shared_settings", "settings_01_demultiplexing"],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/01a_Dorado_Demultiplexing.Rmd config/config.yaml shared_settings settings_01_demultiplexing 2>&1 | tee {log}"


rule stage_01b_obitools:
    input:
        dorado=dorado_fastqs,
        ngsfilter=config["settings_01_demultiplexing"]["ngsfilter_path"],
    output:
        html="Code/HTMLs/01b_ObiTools4_Demultiplexing.html",
        identified=identified_fastqs,
        unidentified=unidentified_fastqs,
    log:
        "logs/stage_01b.log",
    benchmark:
        "benchmarks/stage_01b.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/01b_ObiTools4_Demultiplexing.Rmd",
            ["shared_settings", "settings_01_demultiplexing"],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/01b_ObiTools4_Demultiplexing.Rmd config/config.yaml shared_settings settings_01_demultiplexing 2>&1 | tee {log}"


rule stage_02_obitools_summaries:
    input:
        identified=identified_fastqs,
        unidentified=unidentified_fastqs,
        barcode_map=config["shared_settings"]["barcode_map_path"],
    output:
        "Code/HTMLs/02_ObiTools4_Basic_Visualisation.html",
    log:
        "logs/stage_02.log",
    benchmark:
        "benchmarks/stage_02.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/02_ObiTools4_Basic_Visualisation.Rmd", ["shared_settings"]
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/02_ObiTools4_Basic_Visualisation.Rmd config/config.yaml shared_settings 2>&1 | tee {log}"


rule stage_03_minimap2:
    input:
        dorado=dorado_fastqs,
        identified=identified_fastqs,
        reference=config["shared_settings"]["reference_path"],
        primer_panel=config["shared_settings"]["primer_panel_path"],
        barcode_map=config["shared_settings"]["barcode_map_path"],
    output:
        html="Code/HTMLs/03_SNP_Genotyping_MiniMap2.html",
        padded=padded_fastqs,
        genotypes=f"{minimap2_folder}/all_barcodes_Minimap2_genotypes.tsv",
        mapping=f"{minimap2_folder}/mapping_summary_03.tsv",
    log:
        "logs/stage_03.log",
    benchmark:
        "benchmarks/stage_03.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/03_SNP_Genotyping_MiniMap2.Rmd",
            ["shared_settings", "settings_03_minimap2"],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/03_SNP_Genotyping_MiniMap2.Rmd config/config.yaml shared_settings settings_03_minimap2 2>&1 | tee {log}"


rule stage_04_asp_genotyping:
    input:
        reads=asp_input_fastqs,
        functions=asp_functions,
        reference=config["shared_settings"]["reference_path"],
        primer_panel=config["shared_settings"]["primer_panel_path"],
        barcode_map=config["shared_settings"]["barcode_map_path"],
    output:
        html="Code/HTMLs/04_SNP_Genotyping_ASP_primer.html",
        asp_genotypes=asp_genotypes_table,
        asp_reads=asp_reads_table,
    log:
        "logs/stage_04.log",
    benchmark:
        "benchmarks/stage_04.tsv"
    threads: config["shared_settings"]["n_cores"]
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/04_SNP_Genotyping_ASP_primer.Rmd",
            ["shared_settings", "settings_04_genotyping"],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/04_SNP_Genotyping_ASP_primer.Rmd config/config.yaml shared_settings settings_04_genotyping 2>&1 | tee {log}"


rule stage_05a_bayesian:
    input:
        asp_reads=asp_reads_table,
        asp_genotypes=asp_genotypes_table,
        stan_model="Stan_config/base_mixture_model_hierarchical.stan",
        primer_panel=config["shared_settings"]["primer_panel_path"],
    output:
        html="Code/HTMLs/05a_SNP_Bayesian_Genotyping.html",
        bayes_genotypes=bayes_genotypes_table,
        samples=f"{bayes_folder}/sample_parameters_hierarchical.tsv",
        stan_results=f"{bayes_folder}/stan_plot_inputs.rds",
    log:
        "logs/stage_05a.log",
    benchmark:
        "benchmarks/stage_05a.tsv"
    threads: config["shared_settings"]["n_cores"]
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/05a_SNP_Bayesian_Genotyping.Rmd",
            ["shared_settings", "settings_05_bayesian"],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/05a_SNP_Bayesian_Genotyping.Rmd config/config.yaml shared_settings settings_05_bayesian 2>&1 | tee {log}"


rule stage_05a_bayesian_genotype_prior:
    input:
        asp_reads=asp_reads_table,
        asp_genotypes=asp_genotypes_table,
        stan_model="Stan_config/base_mixture_model_hierarchical.stan",
        primer_panel=config["shared_settings"]["primer_panel_path"],
    output:
        html=f"Code/HTMLs/05a_SNP_Bayesian_Genotyping_{genotype_prior_type}.html",
        genotypes=f"{genotype_bayes_folder}/all_barcodes_base_mixture_genotypes_hierarchical.tsv",
        samples=f"{genotype_bayes_folder}/sample_parameters_hierarchical.tsv",
        stan_results=f"{genotype_bayes_folder}/stan_plot_inputs.rds",
    log:
        "logs/stage_05a_genotype_prior.log",
    benchmark:
        "benchmarks/stage_05a_genotype_prior.tsv"
    threads: config["shared_settings"]["n_cores"]
    params:
        settings={
            **stage_settings(
                "Code/Rmd_Stages/05a_SNP_Bayesian_Genotyping.Rmd",
                ["shared_settings", "settings_05_bayesian"],
            ),
            "prior_type": genotype_prior_type,
        },
        prior=genotype_prior_type,
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/05a_SNP_Bayesian_Genotyping.Rmd config/config.yaml shared_settings settings_05_bayesian prior_type={params.prior} 2>&1 | tee {log}"


rule stage_05b_downsampling:
    input:
        stan_results=f"{bayes_folder}/stan_plot_inputs.rds",
        asp_reads=asp_reads_table,
        reference_genotypes=reference_genotypes,
        functions=asp_functions,
    output:
        html="Code/HTMLs/05b_Downsampling_Calibration.html",
        downsampling_results=f"{downsampling_folder}/downsampling_plot_inputs.rds",
    log:
        "logs/stage_05b.log",
    benchmark:
        "benchmarks/stage_05b.tsv"
    threads: config["shared_settings"]["n_cores"]
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/05b_Downsampling_Calibration.Rmd",
            ["shared_settings", "settings_05_bayesian"],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/05b_Downsampling_Calibration.Rmd config/config.yaml shared_settings settings_05_bayesian 2>&1 | tee {log}"


rule stage_06a_diagnostics:
    input:
        asp_reads=asp_reads_table,
        asp_genotypes=asp_genotypes_table,
        bayes_genotypes=bayes_genotypes_table,
        reference_genotypes=reference_genotypes,
        raw_reads=raw_reads,
        dorado=dorado_fastqs,
        identified=identified_fastqs,
        unidentified=unidentified_fastqs,
        padded=padded_fastqs,
        reference=config["shared_settings"]["reference_path"],
        primer_panel=config["shared_settings"]["primer_panel_path"],
        functions=asp_functions,
    output:
        html="Code/HTMLs/06a_Pipeline_Diagnostics.html",
        locus_summary=locus_summary,
    log:
        "logs/stage_06a.log",
    benchmark:
        "benchmarks/stage_06a.tsv"
    threads: config["shared_settings"]["n_cores"]
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/06a_Pipeline_Diagnostics.Rmd",
            ["shared_settings", "settings_04_genotyping", "settings_05_bayesian"],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/06a_Pipeline_Diagnostics.Rmd config/config.yaml shared_settings settings_04_genotyping settings_05_bayesian 2>&1 | tee {log}"


rule stage_06b_results:
    input:
        asp_genotypes=asp_genotypes_table,
        bayes_genotypes=bayes_genotypes_table,
        minimap2_genotypes=f"{minimap2_folder}/all_barcodes_Minimap2_genotypes.tsv",
        reference_genotypes=reference_genotypes,
        stan_results=f"{bayes_folder}/stan_plot_inputs.rds",
        downsampling_results=f"{downsampling_folder}/downsampling_plot_inputs.rds",
        locus_summary=locus_summary,
        functions=asp_functions,
        primer_panel=config["shared_settings"]["primer_panel_path"],
        reference=config["shared_settings"]["reference_path"],
    output:
        "Code/HTMLs/06b_Genotyping_Results.html",
    log:
        "logs/stage_06b.log",
    benchmark:
        "benchmarks/stage_06b.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/06b_Genotyping_Results.Rmd",
            ["shared_settings", "settings_05_bayesian"],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/06b_Genotyping_Results.Rmd config/config.yaml shared_settings settings_05_bayesian 2>&1 | tee {log}"


rule stage_06c_speed:
    input:
        asp_reads=asp_reads_table,
        reference_genotypes=reference_genotypes,
        functions=asp_functions,
    output:
        "Code/HTMLs/06c_Speed_Optimisation.html",
    log:
        "logs/stage_06c.log",
    benchmark:
        "benchmarks/stage_06c.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/06c_Speed_Optimisation.Rmd", ["shared_settings"]
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/06c_Speed_Optimisation.Rmd config/config.yaml shared_settings 2>&1 | tee {log}"


rule stage_07_failed_genotypes:
    input:
        genotypes=failed_genotypes_table,
        asp_reads=asp_reads_table,
        reference_genotypes=reference_genotypes,
        reads=asp_input_fastqs,
        barcode_map=config["shared_settings"]["barcode_map_path"],
        primer_panel=config["shared_settings"]["primer_panel_path"],
        reference=config["shared_settings"]["reference_path"],
        functions=asp_functions,
    output:
        html="Code/HTMLs/07_Failed_Genotypes_Analysis.html",
        summary=f"{failed_folder}/Failed_markers_summary.tsv",
    log:
        "logs/stage_07.log",
    benchmark:
        "benchmarks/stage_07.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/07_Failed_Genotypes_Analysis.Rmd",
            [
                "shared_settings",
                "settings_05_bayesian",
                "settings_07_failed_genotypes",
            ],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/07_Failed_Genotypes_Analysis.Rmd config/config.yaml shared_settings settings_05_bayesian settings_07_failed_genotypes 2>&1 | tee {log}"


rule stage_08a_ancestry_inputs:
    input:
        genotypes=f"{genotype_bayes_folder}/all_barcodes_base_mixture_genotypes_hierarchical.tsv",
        sample_parameters=f"{bayes_folder}/sample_parameters_hierarchical.tsv",
        reference_genotypes=reference_genotypes,
        reference_ped=settings_08["reference_ped_path"],
        reference_map=settings_08["reference_map_path"],
        primer_panel=config["shared_settings"]["primer_panel_path"],
        functions=asp_functions,
    output:
        html="Code/HTMLs/08a_Ancestry_Inputs.html",
        ancestry_inputs=f"{ancestry_folder}/ancestry_inputs.rds",
        newhybrids_input=f"{newhybrids_folder}/newhybrids_input.txt",
        newhybrids_key=f"{newhybrids_folder}/newhybrids_key.tsv",
        power_input=f"{newhybrids_power_folder}/newhybrids_power_input.txt",
        power_key=f"{newhybrids_power_folder}/newhybrids_power_key.tsv",
        structure_input=f"{structure_folder}/structure_input.txt",
        mainparams=f"{structure_folder}/mainparams",
        extraparams=f"{structure_folder}/extraparams",
    log:
        "logs/stage_08a.log",
    benchmark:
        "benchmarks/stage_08a.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/08a_Ancestry_Inputs.Rmd",
            [
                "shared_settings",
                "settings_05_bayesian",
                "settings_08_ancestry_validation",
            ],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/08a_Ancestry_Inputs.Rmd config/config.yaml shared_settings settings_05_bayesian settings_08_ancestry_validation 2>&1 | tee {log}"


rule newhybrids:
    input:
        f"{newhybrids_folder}/newhybrids_input.txt",
    output:
        f"{newhybrids_folder}/aa-PofZ.txt",
    log:
        "logs/newhybrids.log",
    benchmark:
        "benchmarks/newhybrids.tsv"
    resources:
        heavy_jobs=1,
    params:
        newhybrids=config["workflow_settings"]["newhybrids_command"],
        timer=config["workflow_settings"]["timer_command"],
    # Newhybrids exits with 1, which is now ignored, results is checked by the files itself
    shell:
        "(cd $(dirname {input}) && ({params.timer} {params.newhybrids} -d newhybrids_input.txt --burn-in 100000 --num-sweeps 500000 "
        "--theta-prior Jeffreys --pi-prior Jeffreys -s 1 2 --no-gui || true) && test -s aa-PofZ.txt) > {log} 2>&1"


rule newhybrids_power:
    input:
        f"{newhybrids_power_folder}/newhybrids_power_input.txt",
    output:
        f"{newhybrids_power_folder}/aa-PofZ.txt",
    log:
        "logs/newhybrids_power.log",
    benchmark:
        "benchmarks/newhybrids_power.tsv"
    resources:
        heavy_jobs=1,
    params:
        newhybrids=config["workflow_settings"]["newhybrids_command"],
        timer=config["workflow_settings"]["timer_command"],
    # Newhybrids exits with 1, which is now ignored, results is checked by the files itself
    shell:
        "(cd $(dirname {input}) && ({params.timer} {params.newhybrids} -d newhybrids_power_input.txt --burn-in 100000 --num-sweeps 500000 "
        "--theta-prior Jeffreys --pi-prior Jeffreys -s 1 2 --no-gui || true) && test -s aa-PofZ.txt) > {log} 2>&1"


rule structure_run:
    input:
        structure_input=f"{structure_folder}/structure_input.txt",
        mainparams=f"{structure_folder}/mainparams",
        extraparams=f"{structure_folder}/extraparams",
    output:
        f"{structure_folder}/structure_K2_run{{run}}_f",
    log:
        "logs/structure_run{run}.log",
    benchmark:
        "benchmarks/structure_run{run}.tsv"
    wildcard_constraints:
        run="[0-9]+",
    resources:
        heavy_jobs=1,
    params:
        structure=config["workflow_settings"]["structure_command"],
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "(cd $(dirname {input.structure_input}) && {params.timer} {params.structure} -m mainparams -e extraparams -K 2 -D {wildcards.run} "
        "-o structure_K2_run{wildcards.run}) > {log} 2>&1"


rule stage_08b_ancestry_results:
    input:
        ancestry_inputs=f"{ancestry_folder}/ancestry_inputs.rds",
        sample_parameters=f"{bayes_folder}/sample_parameters_hierarchical.tsv",
        newhybrids_key=f"{newhybrids_folder}/newhybrids_key.tsv",
        power_key=f"{newhybrids_power_folder}/newhybrids_power_key.tsv",
        newhybrids=f"{newhybrids_folder}/aa-PofZ.txt",
        newhybrids_power=f"{newhybrids_power_folder}/aa-PofZ.txt",
        structure=structure_outputs,
        functions=asp_functions,
    output:
        "Code/HTMLs/08b_Ancestry_Results.html",
    log:
        "logs/stage_08b.log",
    benchmark:
        "benchmarks/stage_08b.tsv"
    params:
        settings=stage_settings(
            "Code/Rmd_Stages/08b_Ancestry_Results.Rmd",
            [
                "shared_settings",
                "settings_05_bayesian",
                "settings_08_ancestry_validation",
            ],
        ),
        timer=config["workflow_settings"]["timer_command"],
    shell:
        "{params.timer} Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/08b_Ancestry_Results.Rmd config/config.yaml shared_settings settings_05_bayesian settings_08_ancestry_validation 2>&1 | tee {log}"
