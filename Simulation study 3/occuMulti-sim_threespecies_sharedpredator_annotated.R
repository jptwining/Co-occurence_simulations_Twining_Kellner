library(unmarked)
library(ggplot2)
# set working directory and seed
setwd('C:/Users/twininjo/Documents/R/Multi-species Rota simulations')
set.seed(123)


########### Multi-species occupancy (Rota et al.) / species interactions models ###########
########################## SIMULATION STUDY 3 ##################################
###################### THREE SPECIES INTERACTION - WITH EXTRA-PAIR INTERACTIONS ############## 

library('unmarked')
# True parameter values on logit-scale
f1 <- 0  # mean occupancy of species 1
f2 <- 1.1 # mean occupancy of species 2 
f3 <- 1.1  # mean occupancy of species 3
f12 <- -1  # second order interaction term s1 - s2
f13 <- -1 # second order interaction term s1 - s3
f23 <- 0 # second order interaction term s2 - s3
s1b1 <- 1 # beta_1 (slope term) for s1
s2b1 <- 1 # beta_1 (slope term) for s2
s3b1 <- 1 # beta_1 (slope term) for s3
s1a0 <- 0.4075 # mean detection probability of species 1 
s2a0 <- 0.4075  # mean detection probability of species 2 
s3a0 <- 0.4075 # mean detection probability of species 3


# Calculate psi to show results of parameter values
2^3
psi <- numeric(8)
psi[1] <- exp(f1 + f2 + f3 + f12 + f13 + f23)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- exp(f3)
psi[5] <- exp(f1 + f2 + f12)
psi[6] <- exp(f1 + f3 + f13)
psi[7] <- exp(f2 + f3 + f23)
psi[8] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[111]","[100]","[010]","[001]", "[110]", "[101]", "[011]", "[000]")
round(psi, 3)

# Organize effect sizes
ef <- list(state = c(f1,s1b1, f2, s2b1, f3, s3b1, f12, f13), det = c(s1a0, s2a0, s3a0))
# set sites and occasions
J <- 3000 # sites
K <- 5   # occasions

# create some covariates
site_covs <- data.frame(grassland = rnorm(J))
# Blank unmarkedFrame with study design
y <- matrix(NA, J, K)
temp <- unmarkedFrameOccuMulti(list(y,y,y), siteCovs = site_covs)

# Formulas
sf <- c("~1+grassland", "~1+grassland", "~1+grassland", "~1", "~1", "~0", "~0") 
df <- c("~1", "~1", "~1")

# Simulate datasets
umf <- simulate(temp, stateformulas=sf, detformulas=df, coefs=ef, nsim=100)

# Fit true model 
truemod1_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=sf, detformulas=df, data=x)
})

# Get CIs for true model (state model)
truemod1_cis <- lapply(1:length(truemod1_fits), function(i){
  data.frame(confint(truemod1_fits[[i]], type='state'), pars = names(coef(truemod1_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get CIs for true model (observation model)
truemod1_det_cis <- lapply(1:length(truemod1_fits), function(i){
  data.frame(confint(truemod1_fits[[i]], type='det'), pars = names(coef(truemod1_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format and bind
true_out <- do.call(rbind, truemod1_cis)
true_det_out <- do.call(rbind, truemod1_det_cis)
rownames(true_out) <- NULL
rownames(true_det_out) <- NULL
true_out


##### Process results########
# s1_psi/f1
truemod1_f1_est <- sapply(truemod1_fits, function(x) coef(x)[1])
truemod1_f1_CIs <- subset(true_out, pars == "psi([sp1] (Intercept))")
truemod1_s1psi_est <- plogis(truemod1_f1_est)
# bias
f1truemod1_bias = (truemod1_f1_est-f1)
s1psitruemod1_bias = (truemod1_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1truemod1_bias)
mean(s1psitruemod1_bias)
hist(truemod1_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_truemod1_coverage = 1*(truemod1_f1_CIs[2]>=f1&truemod1_f1_CIs[1]<=f1)
f1_truemod1_coverage_table <- table(f1_truemod1_coverage)
f1_truemod1_coverage_table[2]/100

# s1b1
truemod1_s1b1_est <- sapply(truemod1_fits, function(x) coef(x)[2])
truemod1_s1b1_CIs <- subset(true_out, pars == "psi([sp1] grassland)")
# bias
s1b1truemod1_bias = (truemod1_s1b1_est-s1b1)/s1b1
mean(s1b1truemod1_bias)
hist(truemod1_s1b1_est, main = "s1b1", xlab = "s1b1 est")
abline(v=s1b1, col='red', lwd=2)
# coverage
s1b1_truemod1_coverage = 1*(truemod1_s1b1_CIs[2]>=s1b1&truemod1_s1b1_CIs[1]<=s1b1)
s1b1_truemod1_coverage_table <- table(s1b1_truemod1_coverage)
s1b1_truemod1_coverage_table[2]/100

# f2/s2_psi
truemod1_f2_est <- sapply(truemod1_fits, function(x) coef(x)[3])
truemod1_f2_CIs <- subset(true_out, pars == "psi([sp2] (Intercept))")
truemod1__s2psi_est <- plogis(truemod1_f2_est)
# bias
f2truemod1_bias = (truemod1_f2_est-f2)/f2
s2_psitruemod1_bias = (truemod1__s2psi_est-(plogis(f2)))/plogis(f2)
mean(f2truemod1_bias)
mean(s2_psitruemod1_bias)
hist(truemod1_f2_est, main="f2", xlab="f2 estimate")
abline(v=f2, col='red', lwd=2)
# coverage
f2_truemod1_coverage = 1*(truemod1_f2_CIs[2]>=f2&truemod1_f2_CIs[1]<=f2)
f2_truemod1_coverage_table <- table(f2_truemod1_coverage)

# s2b1
truemod1_s2b1_est <- sapply(truemod1_fits, function(x) coef(x)[4])
truemod1_s2b1_CIs <- subset(true_out, pars == "psi([sp2] grassland)")
# bias
s2b1truemod1_bias = (truemod1_s2b1_est-s2b1)/s2b1
mean(s2b1truemod1_bias)
hist(truemod1_s2b1_est, main="s2b1", xlab="s2b1 est")
abline(v=s2b1, col='red', lwd=2)
# coverage
s2b1_truemod1_coverage = 1*(truemod1_s2b1_CIs[2]>=s2b1&truemod1_s2b1_CIs[1]<=s2b1)
s2b1_truemod1_coverage_table <- table(s2b1_truemod1_coverage)


# f3/s3_psi
truemod1_f3_est <- sapply(truemod1_fits, function(x) coef(x)[5])
truemod1_f3_CIs <- subset(true_out, pars == "psi([sp3] (Intercept))")
truemod1__s3psi_est <- plogis(truemod1_f3_est)
# bias
f3truemod1_bias = (truemod1_f3_est-f3)/f3
s3_psitruemod1_bias = (truemod1__s3psi_est-(plogis(f3)))/plogis(f3)
mean(f3truemod1_bias)
mean(s3_psitruemod1_bias)
hist(truemod1_f3_est, main="f3", xlab="f3 estimate")
abline(v=f3, col='red', lwd=2)
# coverage
f3_truemod1_coverage = 1*(truemod1_f3_CIs[2]>=f3&truemod1_f3_CIs[1]<=f3)
f3_truemod1_coverage_table <- table(f3_truemod1_coverage)

# s3b1
truemod1_s3b1_est <- sapply(truemod1_fits, function(x) coef(x)[6])
truemod1_s3b1_CIs <- subset(true_out, pars == "psi([sp3] grassland)")
# bias
s3b1truemod1_bias = (truemod1_s3b1_est-s3b1)/s3b1
mean(s3b1truemod1_bias)
hist(truemod1_s3b1_est, main="s3b1", xlab="s3b1 est")
abline(v=s3b1, col='red', lwd=2)
# coverage
s3b1_truemod1_coverage = 1*(truemod1_s3b1_CIs[2]>=s3b1&truemod1_s3b1_CIs[1]<=s3b1)
s3b1_truemod1_coverage_table <- table(s3b1_truemod1_coverage)


# interaction term f12
truemod1_f12_est <- sapply(truemod1_fits, function(x) coef(x)[7])
truemod1_f12_CIs <- subset(true_out, pars == "psi([sp1:sp2] (Intercept))")
# bias
f12truemod1_bias = (truemod1_f12_est-f12)/f12
mean(f12truemod1_bias)
hist(truemod1_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_truemod1_coverage = 1*(truemod1_f12_CIs[2]>=f12&truemod1_f12_CIs[1]<=f12)
f12_truemod1_coverage_table <- table(f12_truemod1_coverage)

# interaction term f13
truemod1_f13_est <- sapply(truemod1_fits, function(x) coef(x)[8])
truemod1_f13_CIs <- subset(true_out, pars == "psi([sp1:sp3] (Intercept))")
# bias
f13truemod1_bias = (truemod1_f13_est-f13)/f13
mean(f13truemod1_bias)
hist(truemod1_f13_est, main="f13", xlab="f13 estimate")
abline(v=f13, col='red', lwd=2)
# coverage
f13_truemod1_coverage = 1*(truemod1_f13_CIs[2]>=f13&truemod1_f13_CIs[1]<=f13)
f13_truemod1_coverage_table <- table(f13_truemod1_coverage)

# p_s1
# Process results
truemod1_a0_s1_est <- sapply(truemod1_fits, function(x) coef(x)[9])
truemod1_a0_s1_CIs <- subset(true_det_out, pars == "p([sp1] (Intercept))")
truemod1_p_s1_est <- plogis(truemod1_a0_s1_est)
# bias
s1a0truemod1_bias = (truemod1_a0_s1_est-s1a0)/s1a0
s1ptruemod1_bias = (truemod1_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0truemod1_bias)
mean(s1ptruemod1_bias)
hist(truemod1_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=s1a0, col='red', lwd=2)
# coverage
s1a0truemod1_coverage = 1*(truemod1_a0_s1_CIs[2]>=s1a0&truemod1_a0_s1_CIs[1]<s1a0)
s1a0_truemod1_coverage_table <- table(s1a0truemod1_coverage)

# p_s2
# Process results
truemod1_a0_s2_est <- sapply(truemod1_fits, function(x) coef(x)[10])
truemod1_a0_s2_CIs <- subset(true_det_out, pars == "p([sp2] (Intercept))")
truemod1_p_s2_est <- plogis(truemod1_a0_s2_est)
# bias
s2a0truemod1_bias = (truemod1_a0_s2_est-s2a0)/s2a0
s2ptruemod1_bias = (truemod1_p_s2_est-(plogis(s2a0)))/plogis(s2a0)
mean(s2a0truemod1_bias)
mean(s2ptruemod1_bias)

hist(truemod1_p_s2_est, main="p_s2", xlab="p_s2 estimate")
abline(v=plogis(s2a0), col='red', lwd=2)
# coverage
a0s2truemod1_coverage = 1*(truemod1_a0_s2_CIs[2]>=s2a0&truemod1_a0_s2_CIs[1]<s2a0)
a0s2_truemod1_coverage_table <- table(a0s2truemod1_coverage)

# p_s3
# Process results
truemod1_a0_s3_est <- sapply(truemod1_fits, function(x) coef(x)[11])
truemod1_a0_s3_CIs <- subset(true_det_out, pars == "p([sp3] (Intercept))")
truemod1_p_s3_est <- plogis(truemod1_a0_s3_est)
# bias
s3a0truemod1_bias = (truemod1_a0_s3_est-s3a0)/s3a0
s3ptruemod1_bias = (truemod1_p_s3_est-(plogis(s3a0)))/plogis(s3a0)
mean(s3a0truemod1_bias)
mean(s3ptruemod1_bias)

hist(truemod1_p_s3_est, main="p_s3", xlab="p_s3 estimate")
abline(v=plogis(s3a0), col='red', lwd=2)
# coverage
a0s3truemod1_coverage = 1*(truemod1_a0_s3_CIs[2]>=s3a0&truemod1_a0_s3_CIs[1]<s3a0)
a0s3truemod1_coverage_table <- table(a0s3truemod1_coverage)

## # # # # # # #  fit the global model (all pairwise interactions possible) # # # # # # # # # # # 
ef <- list(state = c(f1,s1b1, f2, s2b1, f3, s3b1, f12, f13, f23), det = c(s1a0, s2a0, s3a0))

# Formulas
sf <- c("~1+grassland", "~1+grassland", "~1+grassland", "~1", "~1", "~1", "~0") 
df <- c("~1", "~1", "~1")

# Simulate datasets
umf <- simulate(temp, stateformulas=sf, detformulas=df, coefs=ef, nsim=100)

# Fit global model 
globalmod1_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=sf, detformulas=df, data=x)
})

# Get CIs for global model (state model)
globalmod1_cis <- lapply(1:length(globalmod1_fits), function(i){
  data.frame(confint(globalmod1_fits[[i]], type='state'), pars = names(coef(globalmod1_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get CIs for global model (observation model)
globalmod1_det_cis <- lapply(1:length(globalmod1_fits), function(i){
  data.frame(confint(globalmod1_fits[[i]], type='det'), pars = names(coef(globalmod1_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format and bind 
global_out <- do.call(rbind, globalmod1_cis)
global_det_out <- do.call(rbind, globalmod1_det_cis)
rownames(global_out) <- NULL
rownames(global_det_out) <- NULL
global_out


##### Process results########
# s1_psi/f1
globalmod1_f1_est <- sapply(globalmod1_fits, function(x) coef(x)[1])
globalmod1_f1_CIs <- subset(global_out, pars == "psi([sp1] (Intercept))")
globalmod1_s1psi_est <- plogis(globalmod1_f1_est)
# bias
f1globalmod1_bias = (globalmod1_f1_est-f1)
s1psiglobalmod1_bias = (globalmod1_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1globalmod1_bias)
mean(s1psiglobalmod1_bias)
hist(globalmod1_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_globalmod1_coverage = 1*(globalmod1_f1_CIs[2]>=f1&globalmod1_f1_CIs[1]<=f1)
f1_globalmod1_coverage_table <- table(f1_globalmod1_coverage)


# s1b1
globalmod1_s1b1_est <- sapply(globalmod1_fits, function(x) coef(x)[2])
globalmod1_s1b1_CIs <- subset(global_out, pars == "psi([sp1] grassland)")
# bias
s1b1globalmod1_bias = (globalmod1_s1b1_est-s1b1)/s1b1
mean(s1b1globalmod1_bias)
hist(globalmod1_s1b1_est, main = "s1b1", xlab = "s1b1 est")
abline(v=s1b1, col='red', lwd=2)
# coverage
s1b1_globalmod1_coverage = 1*(globalmod1_s1b1_CIs[2]>=s1b1&globalmod1_s1b1_CIs[1]<=s1b1)
s1b1_globalmod1_coverage_table <- table(s1b1_globalmod1_coverage)

# f2/s2_psi
globalmod1_f2_est <- sapply(globalmod1_fits, function(x) coef(x)[3])
globalmod1_f2_CIs <- subset(global_out, pars == "psi([sp2] (Intercept))")
globalmod1__s2psi_est <- plogis(globalmod1_f2_est)
# bias
f2globalmod1_bias = (globalmod1_f2_est-f2)/f2
s2_psiglobalmod1_bias = (globalmod1__s2psi_est-(plogis(f2)))/plogis(f2)
mean(f2globalmod1_bias)
mean(s2_psiglobalmod1_bias)
hist(globalmod1_f2_est, main="f2", xlab="f2 estimate")
abline(v=f2, col='red', lwd=2)
# coverage
f2_globalmod1_coverage = 1*(globalmod1_f2_CIs[2]>=f2&globalmod1_f2_CIs[1]<=f2)
f2_globalmod1_coverage_table <- table(f2_globalmod1_coverage)

# s2b1
globalmod1_s2b1_est <- sapply(globalmod1_fits, function(x) coef(x)[4])
globalmod1_s2b1_CIs <- subset(global_out, pars == "psi([sp2] grassland)")
# bias
s2b1globalmod1_bias = (globalmod1_s2b1_est-s2b1)/s2b1
mean(s2b1globalmod1_bias)
hist(globalmod1_s2b1_est, main="s2b1", xlab="s2b1 est")
abline(v=s2b1, col='red', lwd=2)
# coverage
s2b1_globalmod1_coverage = 1*(globalmod1_s2b1_CIs[2]>=s2b1&globalmod1_s2b1_CIs[1]<=s2b1)
s2b1_globalmod1_coverage_table <- table(s2b1_globalmod1_coverage)


# f3/s3_psi
globalmod1_f3_est <- sapply(globalmod1_fits, function(x) coef(x)[5])
globalmod1_f3_CIs <- subset(global_out, pars == "psi([sp3] (Intercept))")
globalmod1__s3psi_est <- plogis(globalmod1_f3_est)
# bias
f3globalmod1_bias = (globalmod1_f3_est-f3)/f3
s3_psiglobalmod1_bias = (globalmod1__s3psi_est-(plogis(f3)))/plogis(f3)
mean(f3globalmod1_bias)
mean(s3_psiglobalmod1_bias)
hist(globalmod1_f3_est, main="f3", xlab="f3 estimate")
abline(v=f3, col='red', lwd=2)
# coverage
f3_globalmod1_coverage = 1*(globalmod1_f3_CIs[2]>=f3&globalmod1_f3_CIs[1]<=f3)
f3_globalmod1_coverage_table <- table(f3_globalmod1_coverage)

# s3b1
globalmod1_s3b1_est <- sapply(globalmod1_fits, function(x) coef(x)[6])
globalmod1_s3b1_CIs <- subset(global_out, pars == "psi([sp3] grassland)")
# bias
s3b1globalmod1_bias = (globalmod1_s3b1_est-s3b1)/s3b1
mean(s3b1globalmod1_bias)
hist(globalmod1_s3b1_est, main="s3b1", xlab="s3b1 est")
abline(v=s3b1, col='red', lwd=2)
# coverage
s3b1_globalmod1_coverage = 1*(globalmod1_s3b1_CIs[2]>=s3b1&globalmod1_s3b1_CIs[1]<=s3b1)
s3b1_globalmod1_coverage_table <- table(s3b1_globalmod1_coverage)


# interaction term f12
globalmod1_f12_est <- sapply(globalmod1_fits, function(x) coef(x)[7])
globalmod1_f12_CIs <- subset(global_out, pars == "psi([sp1:sp2] (Intercept))")
# bias
f12globalmod1_bias = (globalmod1_f12_est-f12)/f12
mean(f12globalmod1_bias)
hist(globalmod1_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_globalmod1_coverage = 1*(globalmod1_f12_CIs[2]>=f12&globalmod1_f12_CIs[1]<=f12)
f12_globalmod1_coverage_table <- table(f12_globalmod1_coverage)

# interaction term f13
globalmod1_f13_est <- sapply(globalmod1_fits, function(x) coef(x)[8])
globalmod1_f13_CIs <- subset(global_out, pars == "psi([sp1:sp3] (Intercept))")
# bias
f13globalmod1_bias = (globalmod1_f13_est-f13)/f13
mean(f13globalmod1_bias)
hist(globalmod1_f13_est, main="f13", xlab="f13 estimate")
abline(v=f13, col='red', lwd=2)
# coverage
f13_globalmod1_coverage = 1*(globalmod1_f13_CIs[2]>=f13&globalmod1_f13_CIs[1]<=f13)
f13_globalmod1_coverage_table <- table(f13_globalmod1_coverage)

# interaction term f23
globalmod1_f23_est <- sapply(globalmod1_fits, function(x) coef(x)[9])
globalmod1_f23_CIs <- subset(global_out, pars == "psi([sp2:sp3] (Intercept))")
# bias
f23globalmod1_bias = (globalmod1_f23_est-f23)
mean(f23globalmod1_bias)
hist(globalmod1_f23_est, main="f23", xlab="f23 estimate")
abline(v=f23, col='red', lwd=2)
# coverage
f23_globalmod1_coverage = 1*(globalmod1_f23_CIs[2]>=f23&globalmod1_f23_CIs[1]<=f23)
f23_globalmod1_coverage_table <- table(f23_globalmod1_coverage)

# p_s1
# Process results
globalmod1_a0_s1_est <- sapply(globalmod1_fits, function(x) coef(x)[10])
globalmod1_a0_s1_CIs <- subset(global_det_out, pars == "p([sp1] (Intercept))")
globalmod1_p_s1_est <- plogis(globalmod1_a0_s1_est)
# bias
s1a0globalmod1_bias = (globalmod1_a0_s1_est-s1a0)/s1a0
s1pglobalmod1_bias = (globalmod1_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0globalmod1_bias)
mean(s1pglobalmod1_bias)
hist(globalmod1_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=s1a0, col='red', lwd=2)
# coverage
s1a0globalmod1_coverage = 1*(globalmod1_a0_s1_CIs[2]>=s1a0&globalmod1_a0_s1_CIs[1]<s1a0)
s1a0_globalmod1_coverage_table <- table(s1a0globalmod1_coverage)

# p_s2
# Process results
globalmod1_a0_s2_est <- sapply(globalmod1_fits, function(x) coef(x)[11])
globalmod1_a0_s2_CIs <- subset(global_det_out, pars == "p([sp2] (Intercept))")
globalmod1_p_s2_est <- plogis(globalmod1_a0_s2_est)
# bias
s2a0globalmod1_bias = (globalmod1_a0_s2_est-s2a0)/s2a0
s2pglobalmod1_bias = (globalmod1_p_s2_est-(plogis(s2a0)))/plogis(s2a0)
mean(s2a0globalmod1_bias)
mean(s2pglobalmod1_bias)

hist(globalmod1_p_s2_est, main="p_s2", xlab="p_s2 estimate")
abline(v=plogis(s2a0), col='red', lwd=2)
# coverage
a0s2globalmod1_coverage = 1*(globalmod1_a0_s2_CIs[2]>=s2a0&globalmod1_a0_s2_CIs[1]<s2a0)
a0s2_globalmod1_coverage_table <- table(a0s2globalmod1_coverage)

# p_s3
# Process results
globalmod1_a0_s3_est <- sapply(globalmod1_fits, function(x) coef(x)[12])
globalmod1_a0_s3_CIs <- subset(global_det_out, pars == "p([sp3] (Intercept))")
globalmod1_p_s3_est <- plogis(globalmod1_a0_s3_est)
# bias
s3a0globalmod1_bias = (globalmod1_a0_s3_est-s3a0)/s3a0
s3pglobalmod1_bias = (globalmod1_p_s3_est-(plogis(s3a0)))/plogis(s3a0)
mean(s3a0globalmod1_bias)
mean(s3pglobalmod1_bias)

hist(globalmod1_p_s3_est, main="p_s3", xlab="p_s3 estimate")
abline(v=plogis(s3a0), col='red', lwd=2)
# coverage
a0s3globalmod1_coverage = 1*(globalmod1_a0_s3_CIs[2]>=s3a0&globalmod1_a0_s3_CIs[1]<s3a0)
a0s3globalmod1_coverage_table <- table(a0s3globalmod1_coverage)


##################### Misspecify model 1.0 - modelling only interaction between prey species which do not interact ##################################
# misspecify model (missing species interactions between f12 [gamma0=-1], and f13 [gamma0=-1], model only f23 [which = 0])
modelstateform <- c("~1+grassland", "~1+grassland", "~1+grassland", "~0", "~0", "~1", "~0")

# Fit models to check
mismod1_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=modelstateform, detformulas=df, data=x)
})

# Get misspecified model 1's CIs (state model)
mismod1_cis <- lapply(1:length(mismod1_fits), function(i){
  data.frame(confint(mismod1_fits[[i]], type='state'), pars = names(coef(mismod1_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get misspecified model 1's CIs (observation model)
mismod1_det_cis <- lapply(1:length(mismod1_fits), function(i){
  data.frame(confint(mismod1_fits[[i]], type='det'), pars = names(coef(mismod1_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format and bind
mismod1_out <- do.call(rbind, mismod1_cis)
mismod1_det_out <- do.call(rbind, mismod1_det_cis)
rownames(mismod1_out) <- NULL
rownames(mismod1_det_out) <- NULL

# bias and coverage for each parameter
# f1/s1_psi
mismod1_f1_est <- sapply(mismod1_fits, function(x) coef(x)[1])
mismod1_f1_CIs <- subset(mismod1_out, pars == "psi([sp1] (Intercept))")
mismod1_s1psi_est <- plogis(mismod1_f1_est)
# bias
f1mismod1_bias = (mismod1_f1_est-f1)
s1psimismod1_bias = (mismod1_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1mismod1_bias)
mean(s1psimismod1_bias)
hist(mismod1_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_mismod1_coverage = 1*(mismod1_f1_CIs[2]>=f1&mismod1_f1_CIs[1]<=f1)
f1_mismod1_coverage_table <- table(f1_mismod1_coverage)


# s1b1
mismod1_s1b1_est <- sapply(mismod1_fits, function(x) coef(x)[2])
mismod1_s1b1_CIs <- subset(mismod1_out, pars == "psi([sp1] grassland)")
# bias
s1b1mismod1_bias = (mismod1_s1b1_est-s1b1)/s1b1
mean(s1b1mismod1_bias)
hist(mismod1_s1b1_est, main = "s1b1", xlab = "s1b1 est")
abline(v=s1b1, col='red', lwd=2)
# coverage
s1b1_mismod1_coverage = 1*(mismod1_s1b1_CIs[2]>=s1b1&mismod1_s1b1_CIs[1]<=s1b1)
s1b1_mismod1_coverage_table <- table(s1b1_mismod1_coverage)


# f2
mismod1_f2_est <- sapply(mismod1_fits, function(x) coef(x)[3])
mismod1_f2_CIs <- subset(mismod1_out, pars == "psi([sp2] (Intercept))")
mismod1_s2psi_est <- plogis(mismod1_f2_est)
# bias
f2mismod1_bias = (mismod1_f2_est-f2)/f2
s2psimismod1_bias = (mismod1_s2psi_est-(plogis(f2)))/plogis(f2)
# bias
mean(f2mismod1_bias)
mean(s2psimismod1_bias)
hist(mismod1_f2_est, main="f2", xlab="f2 estimate")
abline(v=f2, col='red', lwd=2)
# coverage
f2_mismod1_coverage = 1*(mismod1_f2_CIs[2]>=f2&mismod1_f2_CIs[1]<=f2)
f2_mismod1_coverage_table <- table(f2_mismod1_coverage)

# s2b1
mismod1_s2b1_est <- sapply(mismod1_fits, function(x) coef(x)[4])
mismod1_s2b1_CIs <- subset(mismod1_out, pars == "psi([sp2] grassland)")
# bias
s2b1mismod1_bias = (mismod1_s2b1_est-s2b1)/s2b1
mean(s2b1mismod1_bias)
hist(mismod1_s2b1_est, main="s2b1", xlab="s2b1 est")
abline(v=s2b1, col='red', lwd=2)
# coverage
s2b1_mismod1_coverage = 1*(mismod1_s2b1_CIs[2]>=s2b1&mismod1_s2b1_CIs[1]<=s2b1)
s2b1_mismod1_coverage_table <- table(s2b1_mismod1_coverage)


# f3
mismod1_f3_est <- sapply(mismod1_fits, function(x) coef(x)[5])
mismod1_f3_CIs <- subset(mismod1_out, pars == "psi([sp3] (Intercept))")
mismod1_s3psi_est <- plogis(mismod1_f3_est)
# bias
f3mismod1_bias = (mismod1_f3_est-f3)/f3
s3psimismod1_bias = (mismod1_s3psi_est-(plogis(f3)))/plogis(f3)
# bias
mean(f3mismod1_bias)
mean(s3psimismod1_bias)
# bias
hist(mismod1_f3_est, main="f3", xlab="f3 estimate")
abline(v=f3, col='red', lwd=2)
# coverage
f3_mismod1_coverage = 1*(mismod1_f3_CIs[2]>=f3&mismod1_f3_CIs[1]<=f3)
f3_mismod1_coverage_table <- table(f3_mismod1_coverage)

# s3b1
mismod1_s3b1_est <- sapply(mismod1_fits, function(x) coef(x)[6])
mismod1_s3b1_CIs <- subset(mismod1_out, pars == "psi([sp3] grassland)")
# bias
s3b1mismod1_bias = (mismod1_s3b1_est-s3b1)/s3b1
mean(s3b1mismod1_bias)
hist(mismod1_s3b1_est, main="s3b1", xlab="s3b1 est")
abline(v=s3b1, col='red', lwd=2)
# coverage
s3b1_mismod1_coverage = 1*(mismod1_s3b1_CIs[2]>=s3b1&mismod1_s3b1_CIs[1]<=s3b1)
s3b1_mismod1_coverage_table <- table(s3b1_mismod1_coverage)


# interaction term f23
mismod1_f23_est <- sapply(mismod1_fits, function(x) coef(x)[7])
mismod1_f23_CIs <- subset(mismod1_out, pars == "psi([sp2:sp3] (Intercept))")
# bias
f23mismod1_bias = (mismod1_f23_est-f23)
mean(f23mismod1_bias)
hist(mismod1_f23_est, main="f23", xlab="f23 estimate")
abline(v=f23, col='red', lwd=2)
# coverage
f23_mismod1_coverage = 1*(mismod1_f23_CIs[2]>=f23&mismod1_f23_CIs[1]<=f23)
f23_mismod1_coverage_table <- table(f23_mismod1_coverage)

# Process results
mismod1_a0_s1_est <- sapply(mismod1_fits, function(x) coef(x)[8])
mismod1_a0_s1_CIs <- subset(mismod1_det_out, pars == "p([sp1] (Intercept))")
mismod1_p_s1_est <- plogis(mismod1_a0_s1_est)
# bias
s1a0mismod1_bias = (mismod1_a0_s1_est-s1a0)/s1a0
s1pmismod1_bias = (mismod1_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0mismod1_bias)
mean(s1pmismod1_bias)
hist(mismod1_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=s1a0, col='red', lwd=2)
# coverage
s1a0mismod1_coverage = 1*(mismod1_a0_s1_CIs[2]>=s1a0&mismod1_a0_s1_CIs[1]<s1a0)
s1a0_mismod1_coverage_table <- table(s1a0mismod1_coverage)

# p_s2
# Process results
mismod1_a0_s2_est <- sapply(mismod1_fits, function(x) coef(x)[9])
mismod1_a0_s2_CIs <- subset(mismod1_det_out, pars == "p([sp2] (Intercept))")
mismod1_p_s2_est <- plogis(mismod1_a0_s2_est)
# bias
s2a0mismod1_bias = (mismod1_a0_s2_est-s2a0)/s2a0
s2pmismod1_bias = (mismod1_p_s2_est-(plogis(s2a0)))/plogis(s2a0)
mean(s2a0mismod1_bias)
mean(s2pmismod1_bias)

hist(mismod1_p_s2_est, main="p_s2", xlab="p_s2 estimate")
abline(v=plogis(s2a0), col='red', lwd=2)
# coverage
a0s2mismod1_coverage = 1*(mismod1_a0_s2_CIs[2]>=s2a0&mismod1_a0_s2_CIs[1]<s2a0)
a0s2_mismod1_coverage_table <- table(a0s2mismod1_coverage)

# p_s3
# Process results
mismod1_a0_s3_est <- sapply(mismod1_fits, function(x) coef(x)[10])
mismod1_a0_s3_CIs <- subset(mismod1_det_out, pars == "p([sp3] (Intercept))")
mismod1_p_s3_est <- plogis(mismod1_a0_s3_est)
# bias
s3a0mismod1_bias = (mismod1_a0_s3_est-s3a0)/s3a0
s3pmismod1_bias = (mismod1_p_s3_est-(plogis(s3a0)))/plogis(s3a0)
mean(s3a0mismod1_bias)
mean(s3pmismod1_bias)

hist(mismod1_p_s3_est, main="p_s3", xlab="p_s3 estimate")
abline(v=plogis(s3a0), col='red', lwd=2)
# coverage
a0s3mismod1_coverage = 1*(mismod1_a0_s3_CIs[2]>=s3a0&mismod1_a0_s3_CIs[1]<s3a0)
a0s3mismod1_coverage_table <- table(a0s3mismod1_coverage)


######################### Misspecify model 2.0 - modelling only interaction between predator and prey species 1 ##################################
# Misspecify model by modelling only interaction between predator and prey 1 (f12)
# misspecify model (missing species interactions between f13 [true f13=-1])
modelstateform <- c("~1+grassland", "~1+grassland", "~1+grassland", "~1", "~0", "~0", "~0")

# Fit models
mismod2_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=modelstateform, detformulas=df, data=x)
})

# Get misspecified model 2's CIs (state model)
mismod2_cis <- lapply(1:length(mismod2_fits), function(i){
  data.frame(confint(mismod2_fits[[i]], type='state'), pars = names(coef(mismod2_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})
mismod2_det_cis <- lapply(1:length(mismod2_fits), function(i){
  data.frame(confint(mismod2_fits[[i]], type='det'), pars = names(coef(mismod2_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

mismod2_out <- do.call(rbind, mismod2_cis)
mismod2_det_out <- do.call(rbind, mismod2_det_cis)
rownames(mismod2_out) <- NULL
rownames(mismod2_det_out) <- NULL

# f1/s1_psi
mismod2_f1_est <- sapply(mismod2_fits, function(x) coef(x)[1])
mismod2_f1_CIs <- subset(mismod2_out, pars == "psi([sp1] (Intercept))")
mismod2_s1psi_est <- plogis(mismod2_f1_est)
# bias
f1mismod2_bias = (mismod2_f1_est-f1)
s1psimismod2_bias = (mismod2_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1mismod2_bias)
mean(s1psimismod2_bias)
hist(mismod2_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_mismod2_coverage = 1*(mismod2_f1_CIs[2]>=f1&mismod2_f1_CIs[1]<=f1)
f1_mismod2_coverage_table <- table(f1_mismod2_coverage)


# s1b1
mismod2_s1b1_est <- sapply(mismod2_fits, function(x) coef(x)[2])
mismod2_s1b1_CIs <- subset(mismod2_out, pars == "psi([sp1] grassland)")
# bias
s1b1mismod2_bias = (mismod2_s1b1_est-s1b1)/s1b1
mean(s1b1mismod2_bias)
hist(mismod2_s1b1_est, main = "s1b1", xlab = "s1b1 est")
abline(v=s1b1, col='red', lwd=2)
# coverage
s1b1_mismod2_coverage = 1*(mismod2_s1b1_CIs[2]>=s1b1&mismod2_s1b1_CIs[1]<=s1b1)
s1b1_mismod2_coverage_table <- table(s1b1_mismod2_coverage)


# f2
mismod2_f2_est <- sapply(mismod2_fits, function(x) coef(x)[3])
mismod2_f2_CIs <- subset(mismod2_out, pars == "psi([sp2] (Intercept))")
mismod2_s2psi_est <- plogis(mismod2_f2_est)
# bias
f2mismod2_bias = (mismod2_f2_est-f2)/f2
s2psimismod2_bias = (mismod2_s2psi_est-(plogis(f2)))/plogis(f2)
# bias
mean(f2mismod2_bias)
mean(s2psimismod2_bias)
hist(mismod2_f2_est, main="f2", xlab="f2 estimate")
abline(v=f2, col='red', lwd=2)
# coverage
f2_mismod2_coverage = 1*(mismod2_f2_CIs[2]>=f2&mismod2_f2_CIs[1]<=f2)
f2_mismod2_coverage_table <- table(f2_mismod2_coverage)

# s2b1
mismod2_s2b1_est <- sapply(mismod2_fits, function(x) coef(x)[4])
mismod2_s2b1_CIs <- subset(mismod2_out, pars == "psi([sp2] grassland)")
# bias
s2b1mismod2_bias = (mismod2_s2b1_est-s2b1)/s2b1
mean(s2b1mismod2_bias)
hist(mismod2_s2b1_est, main="s2b1", xlab="s2b1 est")
abline(v=s2b1, col='red', lwd=2)
# coverage
s2b1_mismod2_coverage = 1*(mismod2_s2b1_CIs[2]>=s2b1&mismod2_s2b1_CIs[1]<=s2b1)
s2b1_mismod2_coverage_table <- table(s2b1_mismod2_coverage)


# f3
mismod2_f3_est <- sapply(mismod2_fits, function(x) coef(x)[5])
mismod2_f3_CIs <- subset(mismod2_out, pars == "psi([sp3] (Intercept))")
mismod2_s3psi_est <- plogis(mismod2_f3_est)
# bias
f3mismod2_bias = (mismod2_f3_est-f3)/f3
s3psimismod2_bias = (mismod2_s3psi_est-(plogis(f3)))/plogis(f3)
# bias
mean(f3mismod2_bias)
mean(s3psimismod2_bias)
# bias
hist(mismod2_f3_est, main="f3", xlab="f3 estimate")
abline(v=f3, col='red', lwd=2)
# coverage
f3_mismod2_coverage = 1*(mismod2_f3_CIs[2]>=f3&mismod2_f3_CIs[1]<=f3)
f3_mismod2_coverage_table <- table(f3_mismod2_coverage)

# s3b1
mismod2_s3b1_est <- sapply(mismod2_fits, function(x) coef(x)[6])
mismod2_s3b1_CIs <- subset(mismod2_out, pars == "psi([sp3] grassland)")
# bias
s3b1mismod2_bias = (mismod2_s3b1_est-s3b1)/s3b1
mean(s3b1mismod2_bias)
hist(mismod2_s3b1_est, main="s3b1", xlab="s3b1 est")
abline(v=s3b1, col='red', lwd=2)
# coverage
s3b1_mismod2_coverage = 1*(mismod2_s3b1_CIs[2]>=s3b1&mismod2_s3b1_CIs[1]<=s3b1)
s3b1_mismod2_coverage_table <- table(s3b1_mismod2_coverage)


# interaction term f12
mismod2_f12_est <- sapply(mismod2_fits, function(x) coef(x)[7])
mean(mismod2_f12_est)

mismod2_f12_CIs <- subset(mismod2_out, pars == "psi([sp1:sp2] (Intercept))")
# bias
f12mismod2_bias = (mismod2_f12_est-f12)
mean(f12mismod2_bias)
hist(mismod2_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_mismod2_coverage = 1*(mismod2_f12_CIs[2]>=f12&mismod2_f12_CIs[1]<=f12)
f12_mismod2_coverage_table <- table(f12_mismod2_coverage)

# Process results
mismod2_a0_s1_est <- sapply(mismod2_fits, function(x) coef(x)[8])
mismod2_a0_s1_CIs <- subset(mismod2_det_out, pars == "p([sp1] (Intercept))")
mismod2_p_s1_est <- plogis(mismod2_a0_s1_est)
# bias
s1a0mismod2_bias = (mismod2_a0_s1_est-s1a0)/s1a0
s1pmismod2_bias = (mismod2_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0mismod2_bias)
mean(s1pmismod2_bias)
hist(mismod2_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=s1a0, col='red', lwd=2)
# coverage
s1a0mismod2_coverage = 1*(mismod2_a0_s1_CIs[2]>=s1a0&mismod2_a0_s1_CIs[1]<s1a0)
s1a0_mismod2_coverage_table <- table(s1a0mismod2_coverage)

# p_s2
# Process results
mismod2_a0_s2_est <- sapply(mismod2_fits, function(x) coef(x)[9])
mismod2_a0_s2_CIs <- subset(mismod2_det_out, pars == "p([sp2] (Intercept))")
mismod2_p_s2_est <- plogis(mismod2_a0_s2_est)
# bias
s2a0mismod2_bias = (mismod2_a0_s2_est-s2a0)/s2a0
s2pmismod2_bias = (mismod2_p_s2_est-(plogis(s2a0)))/plogis(s2a0)
mean(s2a0mismod2_bias)
mean(s2pmismod2_bias)

hist(mismod2_p_s2_est, main="p_s2", xlab="p_s2 estimate")
abline(v=plogis(s2a0), col='red', lwd=2)
# coverage
a0s2mismod2_coverage = 1*(mismod2_a0_s2_CIs[2]>=s2a0&mismod2_a0_s2_CIs[1]<s2a0)
a0s2_mismod2_coverage_table <- table(a0s2mismod2_coverage)

# p_s3
# Process results
mismod2_a0_s3_est <- sapply(mismod2_fits, function(x) coef(x)[10])
mismod2_a0_s3_CIs <- subset(mismod2_det_out, pars == "p([sp3] (Intercept))")
mismod2_p_s3_est <- plogis(mismod2_a0_s3_est)
# bias
s3a0mismod2_bias = (mismod2_a0_s3_est-s3a0)/s3a0
s3pmismod2_bias = (mismod2_p_s3_est-(plogis(s3a0)))/plogis(s3a0)
mean(s3a0mismod2_bias)
mean(s3pmismod2_bias)

hist(mismod2_p_s3_est, main="p_s3", xlab="p_s3 estimate")
abline(v=plogis(s3a0), col='red', lwd=2)
# coverage
a0s3mismod2_coverage = 1*(mismod2_a0_s3_CIs[2]>=s3a0&mismod2_a0_s3_CIs[1]<s3a0)
a0s3mismod2_coverage_table <- table(a0s3mismod2_coverage)

# combine all bias from all simulations into a df
all_biasandcov_df <- as.data.frame(rbind(cbind(parameter = "s1_psi", spec = "correct specification", bias = f1truemod1_bias),
                                         cbind(parameter = "s2_psi", spec = "correct specification", bias = f2truemod1_bias),
                                         cbind(parameter = "s3_psi", spec = "correct specification", bias = f3truemod1_bias),
                                         cbind(parameter = "s1b1", spec = "correct specification", bias = s1b1truemod1_bias),
                                         cbind(parameter = "s2b1", spec = "correct specification", bias = s2b1truemod1_bias),
                                         cbind(parameter = "s3b1", spec = "correct specification", bias = s3b1truemod1_bias),
                                         cbind(parameter = "f12", spec = "correct specification", bias = f12truemod1_bias),
                                         cbind(parameter = "f13", spec = "correct specification", bias = f13truemod1_bias),
                                         cbind(parameter = "f23", spec = "correct specification", bias = NA),
                                         cbind(parameter = "s1_psi", spec = "all interaction terms", bias = f1globalmod1_bias),
                                         cbind(parameter = "s2_psi", spec = "all interaction terms", bias = f2globalmod1_bias),
                                         cbind(parameter = "s3_psi", spec = "all interaction terms", bias = f3globalmod1_bias),
                                         cbind(parameter = "s1b1", spec = "all interaction terms", bias = s1b1globalmod1_bias),
                                         cbind(parameter = "s2b1", spec = "all interaction terms", bias = s2b1globalmod1_bias),
                                         cbind(parameter = "s3b1", spec = "all interaction terms", bias = s3b1globalmod1_bias),
                                         cbind(parameter = "f12", spec = "all interaction terms", bias = f12globalmod1_bias),
                                         cbind(parameter = "f13", spec = "all interaction terms", bias = f13globalmod1_bias),
                                         cbind(parameter = "f23", spec = "all interaction terms", bias = f23globalmod1_bias),
                                         cbind(parameter = "s1_psi", spec = "prey interactions only (f23)", bias = f1mismod1_bias),
                                         cbind(parameter = "s2_psi", spec = "prey interactions only (f23)", bias = f2mismod1_bias),
                                         cbind(parameter = "s3_psi", spec = "prey interactions only (f23)", bias = f3mismod1_bias),
                                         cbind(parameter = "s1b1", spec = "prey interactions only (f23)", bias = s1b1mismod1_bias),
                                         cbind(parameter = "s2b1", spec = "prey interactions only (f23)", bias = s2b1mismod1_bias),
                                         cbind(parameter = "s3b1", spec = "prey interactions only (f23)", bias = s3b1mismod1_bias),
                                         cbind(parameter = "f12", spec = "prey interactions only (f23)", bias = NA),
                                         cbind(parameter = "f13", spec = "prey interactions only (f23)", bias = NA),
                                         cbind(parameter = "f23", spec = "prey interactions only (f23)", bias = f23mismod1_bias),
                                         cbind(parameter = "s1_psi", spec = "one predator-prey interaction (f12)", bias = f1mismod2_bias),
                                         cbind(parameter = "s2_psi", spec = "one predator-prey interaction (f12)", bias = f2mismod2_bias),
                                         cbind(parameter = "s3_psi", spec = "one predator-prey interaction (f12)", bias = f3mismod2_bias),
                                         cbind(parameter = "s1b1", spec = "one predator-prey interaction (f12)", bias = s1b1mismod2_bias),
                                         cbind(parameter = "s2b1", spec = "one predator-prey interaction (f12)", bias = s2b1mismod2_bias),
                                         cbind(parameter = "s3b1", spec = "one predator-prey interaction (f12)", bias = s3b1mismod2_bias),
                                         cbind(parameter = "f12", spec = "one predator-prey interaction (f12)", bias = f12mismod2_bias),
                                         cbind(parameter = "f13", spec = "one predator-prey interaction (f12)", bias = NA),
                                         cbind(parameter = "f23", spec = "one predator-prey interaction (f12)", bias = NA)))

all_biasandcov_df$bias <- as.numeric(all_biasandcov_df$bias)

# replace parameter names with parse format style for mathematical symbols in figures
all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_psi"),]$parameter  <- "beta[0]^~s1"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_psi"),]$parameter  <- "beta[0]^~s2"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s3_psi"),]$parameter  <- "beta[0]^~s3"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s1b1"),]$parameter  <- "beta[1]^~s1"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s2b1"),]$parameter  <- "beta[1]^~s2"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s3b1"),]$parameter  <- "beta[1]^~s3"
all_biasandcov_df[which(all_biasandcov_df$parameter == "f12"),]$parameter  <- "f[12]"
all_biasandcov_df[which(all_biasandcov_df$parameter == "f13"),]$parameter  <- "f[13]"
all_biasandcov_df[which(all_biasandcov_df$parameter == "f23"),]$parameter  <- "f[23]"


# set as factor and set levels
all_biasandcov_df$parameter <- factor(all_biasandcov_df$parameter, levels = c( 'f[12]', 'f[13]', 'f[23]', "beta[0]^~s1", "beta[0]^~s2","beta[0]^~s3", "beta[1]^~s1", "beta[1]^~s2", 'beta[1]^~s3'))
all_biasandcov_df$spec <- factor(all_biasandcov_df$spec, levels = c('correct specification', 'prey interactions only (f23)', 'one predator-prey interaction (f12)', 'all interaction terms'))

# data summary function for violin plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# custom colour palette
cbPalette <- c( "#E69F00", "#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# check where we are saving out the figure
getwd()
tiff("cooccu_sim_threespeciesinteraction_extrapairinteractions_violinplot_remake.tiff", units="in", width=4, height=12, res=300)

ggplot(all_biasandcov_df, aes(x = parameter, y = bias, fill = spec)) + geom_violin() + stat_summary(fun.data=data_summary, size = 0.2, position=position_dodge(width=0.9)) +scale_fill_manual(values=cbPalette, name = 'model specification') +
  theme_classic()  + xlab('model parameter')+ geom_hline(data=all_biasandcov_df, aes(yintercept=0)) +   scale_x_discrete(labels = scales::parse_format()) +
  theme(legend.position="bottom") + facet_grid(all_biasandcov_df$spec) + guides(fill = guide_legend(nrow = 2))


dev.off()

# save out / load the simulations
# save.image('occu_multi_sharedpredator_simulationsrun.RData')
load('occu_multi_sharedpredator_simulationsrun.RData')

# combine all bias and coverage
all_biasandcov_df <- as.data.frame(rbind(cbind(parameter = "s1_psi", spec = "correct specification", meanbias = mean(f1truemod1_bias), coverage = f1_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "s2_psi", spec = "correct specification", meanbias = mean(f2truemod1_bias), coverage =  f2_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "s3_psi", spec = "correct specification", meanbias = mean(f3truemod1_bias), coverage =  f3_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "s1b1", spec = "correct specification", meanbias = mean(s1b1truemod1_bias), coverage = s1b1_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "s2b1", spec = "correct specification", meanbias = mean(s2b1truemod1_bias), coverage = s2b1_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "s3b1", spec = "correct specification", meanbias = mean(s3b1truemod1_bias), coverage = s3b1_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "f12", spec = "correct specification", meanbias = mean(f12truemod1_bias), coverage = f12_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "f13", spec = "correct specification", meanbias = mean(f13truemod1_bias), coverage = f13_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "f23", spec = "correct specification", meanbias = NA, coverage = NA),
                                         cbind(parameter = "s1_p", spec = "correct specification", meanbias = mean(s1ptruemod1_bias), coverage = s1a0_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "s2_p", spec = "correct specification", meanbias = mean(s2ptruemod1_bias), coverage = a0s2_truemod1_coverage_table[2]/100),
                                         cbind(parameter = "s3_p", spec = "correct specification", meanbias = mean(s3ptruemod1_bias), coverage = a0s3truemod1_coverage_table[2]/100),
                                         cbind(parameter = "s1_psi", spec = "all interaction terms", meanbias = mean(f1globalmod1_bias), coverage = f1_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s2_psi", spec = "all interaction terms", meanbias = mean(f2globalmod1_bias), coverage = f2_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s3_psi", spec = "all interaction terms", meanbias = mean(f3globalmod1_bias), coverage = f3_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s1b1", spec = "all interaction terms", meanbias = mean(s1b1globalmod1_bias), coverage = s1b1_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s2b1", spec = "all interaction terms", meanbias = mean(s2b1globalmod1_bias), coverage = s2b1_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s3b1", spec = "all interaction terms", meanbias = mean(s3b1globalmod1_bias), coverage = s3b1_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "f12", spec = "all interaction terms", meanbias = mean(f12globalmod1_bias), coverage = f12_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "f13", spec = "all interaction terms", meanbias = mean(f13globalmod1_bias), coverage = f13_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "f23", spec = "all interaction terms", meanbias = mean(f23globalmod1_bias), coverage = f23_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s1_p", spec = "all interaction terms", meanbias = mean(s1pglobalmod1_bias), coverage = s1a0_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s2_p", spec = "all interaction terms", meanbias = mean(s2pglobalmod1_bias), coverage = a0s2_globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s3_p", spec = "all interaction terms", meanbias = mean(s3pglobalmod1_bias), coverage = a0s3globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s1_psi", spec = "prey interactions only (f23)", meanbias = mean(f1mismod1_bias), coverage = f1_mismod1_coverage_table[2]/100), 
                                         cbind(parameter = "s2_psi", spec = "prey interactions only (f23)", meanbias = mean(f2mismod1_bias), coverage = f2_mismod1_coverage_table[2]/100),
                                         cbind(parameter = "s3_psi", spec = "prey interactions only (f23)", meanbias = mean(f3mismod1_bias), coverage = f3_mismod1_coverage_table[2]/100),
                                         cbind(parameter = "s1b1", spec = "prey interactions only (f23)", meanbias = mean(s1b1mismod1_bias), coverage = s1b1_mismod1_coverage_table[2]/100),
                                         cbind(parameter = "s2b1", spec = "prey interactions only (f23)", meanbias = mean(s2b1mismod1_bias), coverage = s2b1_mismod1_coverage_table[2]/100), 
                                         cbind(parameter = "s3b1", spec = "prey interactions only (f23)", meanbias = mean(s3b1mismod1_bias), coverage = s3b1_mismod1_coverage_table[2]/100),
                                         cbind(parameter = "f12", spec = "prey interactions only (f23)", meanbias = NA, coverage = NA),
                                         cbind(parameter = "f13", spec = "prey interactions only (f23)", meanbias = NA, coverage =NA),
                                         cbind(parameter = "f23", spec = "prey interactions only (f23)", meanbias = mean(f23mismod1_bias), coverage = f23_mismod1_coverage_table[2]/100),
                                         cbind(parameter = "s1_p", spec = "prey interactions only (f23)", meanbias = mean(s1pmismod1_bias), coverage = s1a0_mismod1_coverage_table[2]/100),
                                         cbind(parameter = "s2_p", spec = "prey interactions only (f23)", meanbias = mean(s2pmismod1_bias), coverage = a0s2_mismod1_coverage_table[2]/100),
                                         cbind(parameter = "s3_p", spec = "prey interactions only (f23)", meanbias = mean(s3pmismod1_bias), coverage = a0s3globalmod1_coverage_table[2]/100),
                                         cbind(parameter = "s1_psi", spec = "one predator-prey interaction (f12)", meanbias = mean(f1mismod2_bias), coverage = f1_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "s2_psi", spec = "one predator-prey interaction (f12)", meanbias = mean(f2mismod2_bias), coverage = f2_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "s3_psi", spec = "one predator-prey interaction (f12)", meanbias = mean(f3mismod2_bias), coverage = f3_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "s1b1", spec = "one predator-prey interaction (f12)", meanbias = mean(s1b1mismod2_bias), coverage = s1b1_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "s2b1", spec = "one predator-prey interaction (f12)", meanbias = mean(s2b1mismod2_bias), coverage = s2b1_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "s3b1", spec = "one predator-prey interaction (f12)", meanbias = mean(s3b1mismod2_bias), coverage = s3b1_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "f12", spec = "one predator-prey interaction (f12)", meanbias = mean(f12mismod2_bias), coverage = f12_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "f13", spec = "one predator-prey interaction (f12)", meanbias = NA, coverage = NA),
                                         cbind(parameter = "f23", spec = "one predator-prey interaction (f12)", meanbias = NA, coverage = NA),
                                         cbind(parameter = "s1_p", spec = "one predator-prey interaction (f12)", meanbias = mean(s1pmismod2_bias), coverage = s1a0_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "s2_p", spec = "one predator-prey interaction (f12)", meanbias = mean(s2pmismod2_bias), coverage = a0s2_mismod2_coverage_table[2]/100),
                                         cbind(parameter = "s3_p", spec = "one predator-prey interaction (f12)", meanbias = mean(s3pmismod2_bias), coverage = a0s3mismod2_coverage_table[2]/100)))
                                         

all_biasandcov_df$meanbias <- as.numeric(all_biasandcov_df$meanbias)
all_biasandcov_df$meanbias <- round(all_biasandcov_df$meanbias, digits = 2)

write.csv(all_biasandcov_df, 'multispecies_coocc_extrapairinteractions.csv')
