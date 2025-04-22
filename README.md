
# Summary
Contained within this repository are:
1. Scripts for a literature view  which organize papers from the scientific literature which have applied the Rota et al. 2016 multi-species occupancy model for interacting species.
2. Simulation study 1. power analyses exploring the prevelance of type M errors in a series of scenarios examining different sampling regimes and underlying parameter values for co-occurence models.
3. Simulation study 2. a simulation study exploring the impacts of shared responses to habitat covariates on estimating species interactions. 
4. Simulation study 3. a simulation study exploring the impact of extra-pair interactions when estimating species interactions from co-occurence data. 
5. Simulation study 4. a simulation study exploring how indirect interactions between species can arise indirectly. 

## The working directory

Below you will find descriptions of each folder in this repository and files contained within them.

## The folder directory (./main/...)

This folder has five subfolders

## 1. The literature review (./main/Literature review)

### 1.1. 01_identify_papers_to_include.R 

This script takes the full set of papers citing Rota et al. 2016 (`All_Citing_Papers.csv`), removes papers that didn't actually apply the model, and adds blank columns to be filled in during the literature review.
The output is a new CSV `included_papers.csv`.

### 1.2. 02_included_papers_summary_stats.R

This script calculates key parameter values of interest and summary statistics included in the results section of the paper from the raw data in `included_papers.csv`.

## 2. Simulation study 1 (./main/Simulation study 1)

### 2.1. 01_create_effect_size_scenarios.R 

This script creates parameter values for the various interaction effect sizes used in the simulations, conditional on the baseline occupancy and the number of species.
It also defines other parameter values for the occupancy and detection models.
The parameter values are saved in the files `state_scenarios.Rds` and `det_scenarios.Rds`.

### 2.2. 02_power_simulations_interactions.R 

This script provides code for the power analysis focused on power to detect interactions as a function of sample size and effect size, using the parameter values/effect sizes created by the previous script.
The script also builds the associated figures.

### 2.3. 03_power_simulations_covariates.R 

This script provides code for the power analysis focused on power to detect covariate effects on interactions as a function of sample size and effect size, using the parameter values/effect sizes created by the previous script.
The script also builds the associated figures.

## 3. Simulation study 2 (./main/Simulation study 2)

### 3.1. occuMulti-sim_twospecies_missinghabitatcovs_annotated.R

This script simulates data and fits models for a series of scenarios where two species do not interact directly but have either shared, or opposing responses to a shared environmental covariate

## 4. Simulation study 3 (./main/Simulation study 3)

### 4.1. occuMulti-sim_threespecies_sharedpredator_annotated.R

This script simulates data and fit models for a scenario where there are extra-pair interactions in the system (both modelled and unmodelled). Our simulation is framed around a case model of a generalist 
predator with multiple prey species (but is also relevant to cases of herbivores with multiple resources, or parasites/pathogens with multiple hosts etc.).

## 5. Simulation study 4 (./main/Simulation study 4)

### 5.1. occuMulti-sim_threespecies_cascadinginteractions_annotated.R

This script simulates data and fits models to examine how interactions can arise indirectly. We considered a typical tri-trophic cascade scenario involving three interacting species whereby there are negative
interactions between the three species, with direct interactions between species S1 and S2 (f12 = -1), and between S2 and S3 (f23 = -1), but no direct interactions between S1 and S3 (f13 = 0).

