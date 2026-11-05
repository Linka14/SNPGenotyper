configfile: "config.yaml"

rule all:
    input:
        "Code/HTMLs/02_ObiTools4_Basic_Visualisation.html"

rule stage_02_obitools_summaries:
    output:
        "Code/HTMLs/02_ObiTools4_Basic_Visualisation.html"
    params:
        shared_settings = config["shared_settings"]
    log:
        "logs/stage_02.log"
    shell:
        "Rscript Code/R_scripts/render_stage.R Code/Rmd_Stages/02_ObiTools4_Basic_Visualisation.Rmd config.yaml shared_settings 2>&1 | tee {log}"