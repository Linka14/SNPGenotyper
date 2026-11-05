
##### Runner for one pipeline stage, started by Master_Pipeline.Rmd or snakemake in a fresh R process.

#commandArgs grabs all the arguments, but keeps only the last two from run_stage(), path for the .Rmd stage and temporary .rds file
runner_arguments <- commandArgs(trailingOnly = TRUE)
stage_rmd <- runner_arguments[1]

# Changed to work with both snakemake or ordinary rendering
settings_file <- runner_arguments[2]
if (grepl("\\.rds$", settings_file)) {
  stage_settings <- readRDS(settings_file) # Checks for .rds, if it exists its the old master_pipeline run
} else {
  section_names <- runner_arguments[-(1:2)] # For snakemake, reads the yaml file and converts them to a list
  all_settings <- yaml::read_yaml(settings_file)
  if (!all(section_names %in% names(all_settings))) stop("Section not found in ", settings_file, ": ", paste(setdiff(section_names, names(all_settings)), collapse = ", "))
  stage_settings <- do.call(c, unname(all_settings[section_names]))
}

# keep only the settings this stage declares in its YAML header (shared settings not used by this stage are dropped)
declared_settings <- names(rmarkdown::yaml_front_matter(stage_rmd)$params)

#all the unused settings, which are printed to the console
unused_settings <- setdiff(names(stage_settings), declared_settings)
if (length(unused_settings) > 0) message("Not used by this stage: ", paste(unused_settings, collapse = ", "))

#final stage settings, gives names and values of the settings
stage_settings <- stage_settings[intersect(names(stage_settings), declared_settings)]

#knits every stage 
rmarkdown::render(stage_rmd,
  params = stage_settings, #replaces default settings in YAML
  output_dir = normalizePath("Code/HTMLs", mustWork = FALSE), # full path to Code/HTMLs
  knit_root_dir = getwd(), # knit from Full_workflow, so all paths in the Rmd start there
  envir = new.env(), #fresh environment and printed progress
  quiet = FALSE)