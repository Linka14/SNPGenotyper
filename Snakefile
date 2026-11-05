import glob
import yaml

configfile: "config.yaml"

# Settings or values used to build the file paths are taken from config.yaml
dataset_name = config["shared_settings"]["dataset_name"]
barcode_numbers = [f"{barcode_number:02d}" for barcode_number in range(1, config["shared_settings"]["n_barcodes"] + 1)]
demultiplexed_folder = f"Data/Processed_data/{dataset_name}/Demultiplexed_sequence_data"
raw_reads = glob.glob(f"Data/Raw_data/{dataset_name}/*.fastq") or glob.glob(f"Data/Raw_data/{dataset_name}/*.bam")
dorado_fastqs = expand(demultiplexed_folder + "/Dorado_demultiplexed_by_sample/barcode{barcode}.fastq", barcode = barcode_numbers)
identified_fastqs = expand(demultiplexed_folder + "/Identified_sample_sequences_Obi4/barcode{barcode}_identified.fastq", barcode = barcode_numbers)
unidentified_fastqs = expand(demultiplexed_folder + "/Unidentified_sample_sequences_Obi4/barcode{barcode}_unidentified.fastq", barcode = barcode_numbers)

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
  
rule all:
    input:
        "Code/HTMLs/02_ObiTools4_Basic_Visualisation.html"

rule stage_01a_dorado: 
    input:
      raw_reads = raw_reads
    output:
      html = "Code/HTMLs/01a_Dorado_Demultiplexing.html",
      dorado = dorado_fastqs
    params:
      settings = stage_settings("Code/Rmd_Stages/01a_Dorado_Demultiplexing.Rmd", ["shared_settings", "settings_01_demultiplexing"])
    log:
      "logs/stage_01a.log"
    shell:
        "Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/01a_Dorado_Demultiplexing.Rmd config.yaml shared_settings settings_01_demultiplexing 2>&1 | tee {log}"
      
rule stage_01b_obitools:
    input:
      dorado = dorado_fastqs,
      ngsfilter = config["settings_01_demultiplexing"]["ngsfilter_path"]
    output:
      html = "Code/HTMLs/01b_ObiTools4_Demultiplexing.html",
      identified = identified_fastqs,
      unidentified = unidentified_fastqs
    params:
      settings = stage_settings("Code/Rmd_Stages/01b_ObiTools4_Demultiplexing.Rmd", ["shared_settings", "settings_01_demultiplexing"])
    log:
      "logs/stage_01b.log"
    shell:
      "Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/01b_ObiTools4_Demultiplexing.Rmd config.yaml shared_settings settings_01_demultiplexing 2>&1 | tee {log}"  

rule stage_02_obitools_summaries:
    input: 
      identified = identified_fastqs,
      unidentified = unidentified_fastqs,
      barcode_map = config["shared_settings"]["barcode_map_path"]
    output:
        "Code/HTMLs/02_ObiTools4_Basic_Visualisation.html"
    params:
      settings = stage_settings("Code/Rmd_Stages/02_ObiTools4_Basic_Visualisation.Rmd",["shared_settings"])
    log:
        "logs/stage_02.log"
    shell:
        "Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/02_ObiTools4_Basic_Visualisation.Rmd config.yaml shared_settings 2>&1 | tee {log}"
