## README

Repository for Nanopore SNP Genotyper project using simulated data.

This project is used for genotyping SNP nanopore simulated data on two different datasets ("small" and "big").
Main directory is divided into Data, Nanopore_SNP_genotyper.Rproj, Code (R script and Rmd stages),
Results and Stan_config files.

CODE:
- folder contains all the R scripts and Rmd stages for SNP analysis. Scripts are written in such a way that minimum change is needed when
switching to the bigger dataset. All the output paths and input paths are depended on two parameters dataset_name and n_barcodes
(same for all scripts, always at the top of the script). Changing those two to (small -> big) or (6 -> 24) when running scripts in 
order automatically changes all the paths and creates new folders. Everything is executed from the Master_pipieline.Rmd where 
all the settings are gathered and stages can be run in the background

DATA: 
- data includes supporting data with all the reference genotypes, primer panels, ngs filter files, CanFam3 reference data
and other files which stay unchanged during the analysis.
- Raw data holds the data for small and big dataset. Path to this needs to be changes only in the first (01) script.
- Processed data is seperated into small and big folders but essentialy golds the same folder structure. When switching between dataset,
all the paths (input and ouput) automatically change based on the dataset you are working on

RESULTS:
- are seperated into two folders again, small and big. It contains Genotypes_results, Graphics of those results and simple summaries.

STAN_CONFIG:
- this folder contains all the necessary .stan files for Rstan bayesian genotyping analysis.

 

