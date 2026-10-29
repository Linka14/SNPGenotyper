
##### Runner for one pipeline stage, started by Master_Pipeline.Rmd in a fresh R process.

#commandArgs grabs all the arguments, but keeps only the last two from run_stage(), path for the .Rmd stage and temporary .rds file
runner_arguments <- commandArgs(trailingOnly = TRUE)
stage_rmd <- runner_arguments[1]
stage_settings <- readRDS(runner_arguments[2])

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