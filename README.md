
# Summary
Contained within this repo are 1) a review of the literature identifying papers which have applied the Rota et al. 2016 multi-species occupancy model for interacting species. 2) Simulation study 1. a power analyses exploring the prevelance of type M and type S errors in a series of scenarios examining different sampling regimes and underlying parameter values. 
3) Simulation study 2. a simulation study exploring the impacts of shared responses to habitat covariates on estimating species interactions, 4) Simulation study 3. a simulation study exploring the impact of extra-pair interactions when estimating species interactions from co-occurence data. and 
5) Simulation study 4. a simulation study exploring how indirect interactions between species can arise indirectly. 

# The working directory

Below you fill find descriptions of each folder in this repository and files contained within them.

# The folder directory (./scripts/...)

This folder has five subfolders

# 1. The literature review (./scripts/Literature review)

**1.1. 01_identify_papers_to_include.R **

This script identifies xx

**1.2. 02_included_papers_summary_stats.R **

This script summarizes key parameter values of interest from all the reviewed papers

# 2. Simulation study 1 (./scripts/Simulation study 1)

**2.1. 03_create_effect_size_scenarios.R **

**2.2. 04_power_simulations_interactions.R **

**2.3. 05_power_simulations_covariates.R **

# 3. Simulation study 2 (./scripts/Simulation study 2)

** 3.1. occuMulti-sim_twospecies_missinghabitatcovs_annotated.R

This script simulates data and fits models for a series of scenarios where two species do not interact directly but have either shared, or opposing responses to a shared environmental covariate

# 4. Simulation study 3 (./scripts/Simulation study 3)

** 4.1. occuMulti-sim_threespecies_sharedpredator_annotated.R

This script simulates data and fit models for a scenario where there are extra-pair interactions in the system (both modelled and unmodelled). Our simulation is framed around a case model of a generalist 
predator with multiple prey species (but is also relevant to cases of herbivores with multiple resources, or parasites/pathogens with multiple hosts etc.).

# 5. Simulation study 4 (./scripts/Simulation study 4)

** 5.1. occuMulti-sim_threespecies_cascadinginteractions_annotated.R

This script simulates data and fits models to examine how interactions can arise indirectly. We considered a typical tri-trophic cascade scenario involving three interacting species whereby there are negative
interactions between the three species, with direct interactions between species S1 and S2 (f12 = -1), and between S2 and S3 (f23 = -1), but no direct interactions between S1 and S3 (f13 = 0).

