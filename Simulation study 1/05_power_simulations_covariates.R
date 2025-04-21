library(unmarked)

# species: 2, 3, 5
# sites  : 50, 100, 200, 500, 1000
# occasions: 10
# psi: 0.6
# p  : 0.6
# big, small, and 0 interaction intercept
# big (0.6) and small (0.2) covariate effect on interaction; with standardized (z-score) covariate
# chosen so 2 SD would be 1.2 (as big as the 'big' intecept) and 2 SD would be 0.4 (small intercept)

# Read in saved scenarios
state_scenarios <- readRDS("state_scenarios.Rds")
# We only look at the "high occupancy" scenarios when including a covariate
state_scenarios <- state_scenarios[grepl("highocc", names(state_scenarios))]
saveRDS(state_scenarios, "state_scenarios_covariate.Rds")

det_scenarios <- readRDS("det_scenarios.Rds")

# Get inputs
get_inputs <- function(species = 2, sites = 50){
  
  ymat <- matrix(NA, sites, 10) # fixed occasions at 10
  ylist <- replicate(species, ymat, simplify=FALSE)
  names(ylist) <- paste0("sp", 1:species)
  
  # Site covariate x
  sc <- data.frame(x = rnorm(sites))

  umf <- unmarkedFrameOccuMulti(ylist, siteCovs = sc, maxOrder=2)

  nform <- ncol(umf@fDesign)
  stateformulas <- rep("~1", nform)
  stateformulas[species + 1] <- "~x" # first interaction has a covariate
  detformulas <- rep("~1", species)

  list(umf = umf, stateformulas = stateformulas, detformulas = detformulas)
}

# Get final effect sizes (adding covariate effect)
get_effect_size <- function(species = 2, effectsize = 0, covariate_effect = 0.2){
  stopifnot(species %in% c(2,3,5))
  stopifnot(effectsize %in% c(0, 0.1, 0.25)) # interaction value
  stopifnot(covariate_effect %in% c(0.2, 0.6))
  psi <- "high"
  p <- "highdet"
  if(effectsize == 0){
    effectsize <- "noeffect"
  } else if(effectsize == 0.1){
    effectsize <- "smalleffect"
  } else if(effectsize == 0.25) {
    effectsize <- "bigeffect"
  }
  state_scenario <- paste0(psi, "occ_", species, "sp_", effectsize)
  if(!state_scenario %in% names(state_scenarios)) stop("Incorrect scenario")
  state_effects <- state_scenarios[[state_scenario]]

  # Insert the covariate effect
  if(species == 2){
    # when species = 2 the covariate effect is at the end
    state_effects_new <- c(state_effects, covariate_effect)
  } else {
    split_idx <- species + 1 # index of first interaction term
    # covariate effect is right after the first interaction term
    state_effects_new <- c(state_effects[1:split_idx], covariate_effect,
                           state_effects[(split_idx+1):length(state_effects)])
  }

  det_effects <- rep(det_scenarios[[p]], species)
  effects <- list(state = state_effects_new, det = det_effects)
  effects
}

# Run a power analysis simulation under specific settings
run_scenario <- function(species = 2, sites = 50, effectsize = 0, 
                         covariate_effect = 0.2, nsim = 100){

  inp <- get_inputs(species = species, sites = sites)
  eff <- get_effect_size(species = species, effectsize = effectsize,
                         covariate_effect = covariate_effect)

  suppressMessages(
  powerAnalysis(inp$umf, 
                effects=eff,
                stateformulas=inp$stateformulas,
                detformulas=inp$detformulas, maxOrder=2, nsim=nsim)
  )

}

# Run all settings. Do not run this in parallel due to Armadillo issues.
scenarios <- expand.grid(species = c(2,3,5),
                sites = c(50, 100, 200, 500, 1000),
                effectsize = c(0, 0.1, 0.25),
                covariate_effect = c(0.2, 0.6)
              ) 

power_analyses_covariate <- lapply(1:nrow(scenarios), function(i){
  scen <- scenarios[i,]
  print(scen)
  run_scenario(species = scen$species,
               sites = scen$sites,
               effectsize = scen$effectsize,
               covariate_effect = scen$covariate_effect,
               nsim = 100)
})

saveRDS(power_analyses_covariate, "power_analyses_covariate.Rds")

power_analyses_covariate <- readRDS("power_analyses_covariate.Rds")


power_est <- lapply(power_analyses_covariate, function(x){
  s <- summary(x)
  s[s$Parameter == "[sp1:sp2] x", c("Power", "Type S", "Type M")]
})
power_est <- do.call(rbind, power_est)

power_est <- cbind(scenarios, power_est)

# Figures
library(ggplot2)
library(cowplot)

# Species = 2, interaction = 0
dat_sub <- power_est#[power_est$effectsize == 0&power_est$species ==2,]
dat_sub$species <- factor(paste0("Species = ", dat_sub$species))
dat_sub$effectsize <- factor(paste0("Interact intercept = ", dat_sub$effectsize))
dat_sub$covariate_effect <- factor(dat_sub$covariate_effect)
pl <- ggplot(data = dat_sub, aes(x = sites, y=Power)) +
  geom_point(aes(col=covariate_effect)) +
  geom_line(aes(col=covariate_effect)) +
  facet_grid(rows=vars(species), cols=vars(effectsize)) +
  geom_hline(yintercept=0.8, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  labs(col="Covariate effect") +
  ggtitle("psi = 0.6, p = 0.6")

pl

png("Covariate_power.png", height=7, width=9, units='in', res=300)
pl
dev.off()

pl <- ggplot(data = dat_sub, aes(x = sites, y=`Type M`)) +
  geom_point(aes(col=covariate_effect)) +
  geom_line(aes(col=covariate_effect)) +
  facet_grid(rows=vars(species), cols=vars(effectsize)) +
  geom_hline(yintercept=0.8, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  labs(col="Covariate effect") +
  ggtitle("psi = 0.6, p = 0.6")

pl

png("Covariate_TypeM.png", height=7, width=9, units='in', res=300)
pl
dev.off()
