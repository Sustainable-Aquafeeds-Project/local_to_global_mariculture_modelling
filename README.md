# Local to global nutrient input modelling

This repository holds all the code used to analyse data and produce figures for the publication <link>.

## Project structure

Folder/file name |  Description    
-----------------|--------------------------------------------
_quarto.yml          | 
_targets.yaml        | 
archive              | Contains old files no longer in use
data                 | Input data for model
 └ _general_data     | 
 └ atlantic_salmon   | 
docs-rendered        | Documentation files rendered (.html)
docs-working         | Documentation files in progress (.qmd) 
outputs              | Modelled outputs and other prepared files
 └ data_analysis     | 
 └ farm_data         | 
 └ model_run_outputs | 
 └ sensitivity_data  | 
 └ species_data      | 
renv                 | Project package library
src                  | Global scripts used across the project
targets              | Information and outputs from targets pipelines (only used for model sensitivity analysis)
renv.lock            | Project package lockfile
