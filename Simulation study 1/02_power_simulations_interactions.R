library(unmarked)

# species: 2, 3, 5
# sites  : 50, 100, 200, 500, 1000
# occasions: 5, 10, 30
# psi: 0.2, 0.6
# p  : 0.2, 0.6
# big, small, and no effect size (=absolute value of interaction term)

# Read in saved effect size scenarios
state_scenarios <- readRDS("state_scenarios.Rds")
det_scenarios <- readRDS("det_scenarios.Rds")

# Function to create simulation inputs
get_inputs <- function(species = 2, sites = 50, occasions = 5){
  
  ymat <- matrix(NA, sites, occasions)
  ylist <- replicate(species, ymat, simplify=FALSE)
  names(ylist) <- paste0("sp", 1:species)

  # Unmarked frame
  umf <- unmarkedFrameOccuMulti(ylist, maxOrder=2)

  # Formulas
  nform <- ncol(umf@fDesign)
  stateformulas <- rep("~1", nform)
  detformulas <- rep("~1", species)

  list(umf = umf, stateformulas = stateformulas, detformulas = detformulas)
}

# Retrieve true effect sizes based on simulation settings
get_effect_size <- function(species = 2, psi = 0.2, p = 0.2, effectsize = 0){
  stopifnot(species %in% c(2,3,5))
  stopifnot(psi %in% c(0.2, 0.6))
  stopifnot(p %in% c(0.05, 0.2, 0.6))
  stopifnot(effectsize %in% c(0, 0.1, 0.25))
  psi <- ifelse(psi == 0.2, "low", "high")
  p <- ifelse(p == 0.2, "lowdet", ifelse(p == 0.6, "highdet", "tinydet"))
  if(effectsize == 0){
    effectsize <- "noeffect"
  } else if(effectsize == 0.10){
    effectsize <- "smalleffect"
  } else if(effectsize == 0.25) {
    effectsize <- "bigeffect"
  }
  state_scenario <- paste0(psi, "occ_", species, "sp_", effectsize)
  if(!state_scenario %in% names(state_scenarios)) stop("Incorrect scenario")
  state_effects <- state_scenarios[[state_scenario]]
  det_effects <- rep(det_scenarios[[p]], species)
  effects <- list(state = state_effects, det = det_effects)
  effects
}

# Run a power analysis scenario
run_scenario <- function(species = 2, sites = 50, occasions = 5,
                         psi = 0.2, p = 0.2, effectsize = 0, nsim = 100){

  inp <- get_inputs(species = species, sites = sites, occasions = occasions)
  eff <- get_effect_size(species = species, psi = psi, p = p,
                             effectsize = effectsize)

  suppressMessages(
  powerAnalysis(inp$umf, 
                effects=eff,
                stateformulas=inp$stateformulas,
                detformulas=inp$detformulas, maxOrder=2, nsim=nsim)
  )

}

# Test
#pa50 <- run_scenario(sites = 50, effectsize=0.1)
#pa100 <- run_scenario(sites = 100, effectsize=0.1)
#pa200 <- run_scenario(sites = 200, effectsize=0.1)
#pa500 <- run_scenario(sites = 500, effectsize=0.1)
#pa1000 <- run_scenario(sites = 1000, effectsize=0.1)
#pa_tiny <- run_scenario(sites = 50, effectsize=0.1, p=0.05)
#pl <- unmarkedPowerList(pa50, pa100, pa200, pa500, pa1000)
#plot(pl, param="[sp1:sp2] (Intercept)")

# Full set of scenarios
scenarios <- expand.grid(species = c(2,3,5),
                sites = c(50, 100, 200, 500, 1000),
                occasions = c(5, 10, 30),
                psi = c(0.2, 0.6),
                p = c(0.2, 0.6),
                effectsize = c(0, 0.1, 0.25)
              )

# Don't run this in parallel, Armadillo is doing stuff in parallel at the 
# same time and it causes hangs
power_analyses <- lapply(1:nrow(scenarios), function(i){
  scen <- scenarios[i,]
  print(scen)
  run_scenario(species = scen$species,
               sites = scen$sites,
               occasions = scen$occasions,
               psi = scen$psi,
               p = scen$p,
               effectsize = scen$effectsize,
               nsim = 100)
})

saveRDS(power_analyses, "power_analyses.Rds")
power_analyses <- readRDS("power_analyses.Rds")

# p = 0.05 scenarios (additional scenarios requested by reviewer)
# since we don't end up using the 0 effect size scenarios, they are
# not included here
scenarios_p005 <- expand.grid(species = c(2,3,5),
                    sites = c(50, 100, 200, 500, 1000),
                    occasions = c(5, 10, 30),
                    psi = c(0.2, 0.6),
                    p = 0.05,
                    effectsize = c(0.1, 0.25)
                  )

scenarios_p005 <- scenarios_p005[order(scenarios_p005$species),]
rownames(scenarios_p005) <- NULL

power_analyses_p005 <- lapply(1:nrow(scenarios_p005), function(i){
  scen <- scenarios_p005[i,]
  print(scen)
  run_scenario(species = scen$species,
               sites = scen$sites,
               occasions = scen$occasions,
               psi = scen$psi,
               p = scen$p,
               effectsize = scen$effectsize,
               nsim = 100)
})

saveRDS(power_analyses_p005, "power_analyses_p005.Rds")

# Figures
power_est <- lapply(power_analyses, function(x){
  s <- summary(x)
  s[s$Parameter == "[sp1:sp2] (Intercept)", c("Power", "Type S", "Type M")]
})
power_est <- do.call(rbind, power_est)

power_est <- cbind(scenarios, power_est)

library(ggplot2)
library(cowplot)

# Species = 2
dat_sub <- power_est[power_est$effectsize == 0.1&power_est$species ==2,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl2_1 <- ggplot(data = dat_sub, aes(x = sites, y=Power)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=0.8, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ggtitle("Species = 2, Effect size = 0.1")

dat_sub <- power_est[power_est$effectsize == 0.25&power_est$species ==2,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl2_2 <- ggplot(data = dat_sub, aes(x = sites, y=Power)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=0.8, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ggtitle("Species = 2, Effect size = 0.25")

png("Power_2_species.png", height=9, width=6, units='in', res=300)
plot_grid(pl2_1, pl2_2, nrow=2)
dev.off()

# Species = 3
dat_sub <- power_est[power_est$effectsize == 0.1&power_est$species ==3,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl3_1 <- ggplot(data = dat_sub, aes(x = sites, y=Power)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=0.8, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ggtitle("Species = 3, Effect size = 0.1")

dat_sub <- power_est[power_est$effectsize == 0.25&power_est$species ==3,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl3_2 <- ggplot(data = dat_sub, aes(x = sites, y=Power)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=0.8, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ggtitle("Species = 3, Effect size = 0.25")

png("Power_3_species.png", height=9, width=6, units='in', res=300)
plot_grid(pl3_1, pl3_2, nrow=2)
dev.off()

# Species = 5
dat_sub <- power_est[power_est$effectsize == 0.1&power_est$species ==5,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl5_1 <- ggplot(data = dat_sub, aes(x = sites, y=Power)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=0.8, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank(), strip.background=element_rect("white")) +
  ggtitle("Species = 5, Effect size = 0.1")

dat_sub <- power_est[power_est$effectsize == 0.25&power_est$species ==5,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl5_2 <- ggplot(data = dat_sub, aes(x = sites, y=Power)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=0.8, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank(), strip.background=element_rect("white")) +
  ggtitle("Species = 5, Effect size = 0.25")

png("Power_5_species.png", height=9, width=6, units='in', res=300)
plot_grid(pl5_1, pl5_2, nrow=2)
dev.off()

# Type M

# Species = 2
dat_sub <- power_est[power_est$effectsize == 0.1&power_est$species ==2,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl2_1 <- ggplot(data = dat_sub, aes(x = sites, y=`Type M`)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=1, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ylab("Type M error") +
  ggtitle("Species = 2, Effect size = 0.1")

dat_sub <- power_est[power_est$effectsize == 0.25&power_est$species ==2,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl2_2 <- ggplot(data = dat_sub, aes(x = sites, y=`Type M`)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=1, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ylab("Type M error") +
  ggtitle("Species = 2, Effect size = 0.25")

png("TypeM_2_species.png", height=9, width=6, units='in', res=300)
plot_grid(pl2_1, pl2_2, nrow=2)
dev.off()

# Species = 2
dat_sub <- power_est[power_est$effectsize == 0.1&power_est$species ==3,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl3_1 <- ggplot(data = dat_sub, aes(x = sites, y=`Type M`)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=1, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ylab("Type M error") +
  ylim(1, 3.5) +
  ggtitle("Species = 3, Effect size = 0.1")

dat_sub <- power_est[power_est$effectsize == 0.25&power_est$species ==3,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl3_2 <- ggplot(data = dat_sub, aes(x = sites, y=`Type M`)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=1, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ylab("Type M error") +
  ylim(1, 3.5) +
  ggtitle("Species = 3, Effect size = 0.25")

png("TypeM_3_species.png", height=9, width=6, units='in', res=300)
plot_grid(pl3_1, pl3_2, nrow=2)
dev.off()

# Species = 5
dat_sub <- power_est[power_est$effectsize == 0.1&power_est$species ==5,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl5_1 <- ggplot(data = dat_sub, aes(x = sites, y=`Type M`)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=1, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ylab("Type M error") +
  ylim(1, 5) +
  ggtitle("Species = 5, Effect size = 0.1")

dat_sub <- power_est[power_est$effectsize == 0.25&power_est$species ==5,]
dat_sub$occasions <- factor(dat_sub$occasions)
dat_sub$psi <- factor(paste0("psi = ",dat_sub$psi))
dat_sub$p <- factor(paste0("p = ", dat_sub$p))
pl5_2 <- ggplot(data = dat_sub, aes(x = sites, y=`Type M`)) +
  geom_point(aes(col=occasions)) +
  geom_line(aes(col=occasions)) +
  facet_grid(rows=vars(p), cols=vars(psi)) +
  geom_hline(yintercept=1, linetype=2) + 
  theme_bw(base_size=14) +
  theme(panel.grid=element_blank()) +
  ylab("Type M error") +
  ylim(1, 5) +
  ggtitle("Species = 5, Effect size = 0.25")

png("TypeM_5_species.png", height=9, width=6, units='in', res=300)
plot_grid(pl5_1, pl5_2, nrow=2)
dev.off()
