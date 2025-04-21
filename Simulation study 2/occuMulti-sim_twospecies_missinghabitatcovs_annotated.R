library(unmarked)
library(ggplot2)
# set wd and seed
setwd('C:/Users/twininjo/Documents/R/Multi-species Rota simulations')
set.seed(123)


######## Multi-species occupancy (Rota et al.) / species interactions models #########
########################## SIMULATION STUDY 2 ##################################
###### TWO SPECIES SIMULATIONS - SHARED RESPONSES TO ENVIRONMENTAL COVARIATES ########


# # # # # #  SIMULATION STUDY 2.1 # # # # # # 
############ Shared positive response to habitat covariate ###########

# True parameter values on logit-scale
f1 <- 0 # mean occupancy of species 1
f2 <- 0 # mean occupancy of species 2 
f12 <- 0 # second order interaction term s1 - s2
s1b1 <- 1 # beta_1 (slope term) for s1
s2b1 <- 1 # beta_1 (slope term) for s2
s1a0 <- 0.4075  # mean detection probability of species 1 
s2a0 <- 0.4075 # mean detection probability of species 2 


# Calculate psi to show results of parameter values
psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)


# Organize effect sizes
ef <- list(state = c(f1,s1b1, f2, s2b1, f12), det = c(s1a0, s2a0))

# Blank unmarkedFrame with study design
J <- 1000
K <- 5   # occasions
y <- matrix(NA, J, K)


# create some covariates
site_covs <- data.frame(forest = rnorm(J))

# create unmarked multi dataframe
temp <- unmarkedFrameOccuMulti(list(y,y), siteCovs = site_covs)
str(temp)
head(temp)

# Formulas
sf <- c("~1+forest", "~1+forest", "~1") 
df <- c("~1", "~1")

# Simulate datasets
umf <- simulate(temp, stateformulas=sf, detformulas=df, coefs=ef, nsim=250)
str(umf)

# Fit models 
truemod1_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=sf, detformulas=df, data=x)
})

# Get CIs for true model (state model)
truemod1_cis <- lapply(1:length(truemod1_fits), function(i){
  data.frame(confint(truemod1_fits[[i]], type='state'), pars = names(coef(truemod1_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get CIs for true model (detection model)
truemod1_det_cis <- lapply(1:length(truemod1_fits), function(i){
  data.frame(confint(truemod1_fits[[i]], type='det'), pars = names(coef(truemod1_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format data
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
f1truemod1_bias = (truemod1_f1_est-f1)/f1
s1psitruemod1_bias = (truemod1_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1truemod1_bias)
mean(s1psitruemod1_bias)
hist(truemod1_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_truemod1_coverage = 1*(truemod1_f1_CIs[2]>=f1&truemod1_f1_CIs[1]<=f1)
f1_truemod1_coverage_table <- table(f1_truemod1_coverage)

# s1b1
truemod1_s1b1_est <- sapply(truemod1_fits, function(x) coef(x)[2])
truemod1_s1b1_CIs <- subset(true_out, pars == "psi([sp1] forest)")
# bias
s1b1truemod1_bias = (truemod1_s1b1_est-s1b1)/s1b1
mean(s1b1truemod1_bias)
hist(truemod1_s1b1_est, main = "s1b1", xlab = "s1b1 est")
abline(v=s1b1, col='red', lwd=2)
# coverage
s1b1_truemod1_coverage = 1*(truemod1_s1b1_CIs[2]>=s1b1&truemod1_s1b1_CIs[1]<=s1b1)
s1b1_truemod1_coverage_table <- table(s1b1_truemod1_coverage)

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
truemod1_s2b1_CIs <- subset(true_out, pars == "psi([sp2] forest)")
# bias
s2b1truemod1_bias = (truemod1_s2b1_est-s2b1)/s2b1
mean(s2b1truemod1_bias)
hist(truemod1_s2b1_est, main="s2b1", xlab="s2b1 est")
abline(v=s2b1, col='red', lwd=2)
# coverage
s2b1_truemod1_coverage = 1*(truemod1_s2b1_CIs[2]>=s2b1&truemod1_s2b1_CIs[1]<=s2b1)
s2b1_truemod1_coverage_table <- table(s2b1_truemod1_coverage)


# interaction term f12
truemod1_f12_est <- sapply(truemod1_fits, function(x) coef(x)[5])
truemod1_f12_CIs <- subset(true_out, pars == "psi([sp1:sp2] (Intercept))")
# bias
f12truemod1_bias = (truemod1_f12_est-f12)
mean(f12truemod1_bias)
hist(truemod1_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_truemod1_coverage = 1*(truemod1_f12_CIs[2]>=f12&truemod1_f12_CIs[1]<=f12)
f12_truemod1_coverage_table <- table(f12_truemod1_coverage)

# p_s1
# Process results
truemod1_a0_s1_est <- sapply(truemod1_fits, function(x) coef(x)[6])
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
truemod1_a0_s2_est <- sapply(truemod1_fits, function(x) coef(x)[7])
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

# misspecify model (missing habitat covariate)
modelstateform <- c("~1", "~1", "~1")
str(umf)

# Fit models to check
mismod1_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=modelstateform, detformulas=df, data=x)
})

# Get misspecified model CIs (state model)
mismod1_cis <- lapply(1:length(mismod1_fits), function(i){
  data.frame(confint(mismod1_fits[[i]], type='state'), pars = names(coef(mismod1_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get misspecified model CIs (observation model)
mismod1_det_cis <- lapply(1:length(mismod1_fits), function(i){
  data.frame(confint(mismod1_fits[[i]], type='det'), pars = names(coef(mismod1_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format simulation results 
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
f1mismod1_bias = (mismod1_f1_est-f1)/f1
# s1psimismod1_bias = (mismod1_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1mismod1_bias)
mean(s1psimismod1_bias)
hist(mismod1_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_mismod1_coverage = 1*(mismod1_f1_CIs[2]>=f1&mismod1_f1_CIs[1]<=f1)
f1_mismod1_coverage_table <- table(f1_mismod1_coverage)


# f2
mismod1_f2_est <- sapply(mismod1_fits, function(x) coef(x)[3])
mismod1_f2_CIs <- subset(mismod1_out, pars == "psi([sp2] (Intercept))")
mismod1_s2psi_est <- plogis(mismod1_f2_est)
# bias
f2mismod1_bias = (mismod1_f2_est-f2)/f2
# s2psimismod1_bias = (mismod1_s2psi_est-(plogis(f2)))/plogis(f2)
# bias
mean(f2mismod1_bias)
mean(s2psimismod1_bias)
hist(mismod1_f2_est, main="f2", xlab="f2 estimate")
abline(v=f2, col='red', lwd=2)
# coverage
f2_mismod1_coverage = 1*(mismod1_f2_CIs[2]>=f2&mismod1_f2_CIs[1]<=f2)
f2_mismod1_coverage_table <- table(f2_mismod1_coverage)


# f12
mismod1_f12_est <- sapply(mismod1_fits, function(x) coef(x)[3])
mismod1_f12_CIs <- subset(mismod1_out, pars == "psi([sp1:sp2] (Intercept))")
mean(mismod1_f12_est)
mean(as.matrix(mismod1_f12_CIs[1]))
mean(as.matrix(mismod1_f12_CIs[2]))
# bias
f12mismod1_bias = (mismod1_f12_est-f12)
mean(f12mismod1_bias)
hist(mismod1_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_coverage = 1*(mismod1_f12_CIs[2]>=f12&mismod1_f12_CIs[1]<=f12)
mismod1_f12_coverage <- table(f12_coverage)
(mismod1_f12_coverage[2]/250)*100

# s1_p/a0
mismod1_a0_s1_est <- sapply(mismod1_fits, function(x) coef(x)[4])
mismod1_a0_s1_CIs <- subset(mismod1_det_out, pars == "p([sp1] (Intercept))")
mismod1_p_s1_est <- plogis(mismod1_a0_s1_est)
# bias
s1a0mismod1_bias = (mismod1_a0_s1_est-s1a0)/s1a0
s1pmismod1_bias = (mismod1_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0mismod1_bias)
mean(s1pmismod1_bias)
hist(mismod1_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=plogis(s1a0), col='red', lwd=2)
# coverage
s1a0mismod1_coverage = 1*(mismod1_a0_s1_CIs[2]>=s1a0&mismod1_a0_s1_CIs[1]<s1a0)
s1a0_mismod1_coverage_table <- table(s1a0mismod1_coverage)

# s2_p/a0
mismod1_a0_s2_est <- sapply(mismod1_fits, function(x) coef(x)[5])
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


# combine all bias into a single dataframe
all_biasandcov_df <- as.data.frame(rbind(cbind(parameter = "s1_b0", spec = "true", bias = f1truemod1_bias),
                                         cbind(parameter = "s2_b0", spec = "true", bias = f2truemod1_bias),
                                         cbind(parameter = "s1_b1", spec = "true", bias = s1b1truemod1_bias),
                                         cbind(parameter = "s2_b1", spec = "true", bias = s2b1truemod1_bias),
                                         cbind(parameter = "f12", spec = "true", bias = f12truemod1_bias),
                                         # cbind(parameter = "s1_p", spec = "true", bias = s1ptruemod1_bias),
                                         # cbind(parameter = "s2_p", spec = "true", bias = s2ptruemod1_bias),
                                         cbind(parameter = "s1_b0", spec = "misspecified", bias = f1mismod1_bias),
                                         cbind(parameter = "s2_b0", spec = "misspecified", bias = f2mismod1_bias),
                                         cbind(parameter = "s1_b1", spec = "misspecified", bias = NA),
                                         cbind(parameter = "s2_b1", spec = "misspecified", bias = NA),
                                         cbind(parameter = "f12", spec = "misspecified", bias = f12mismod1_bias)))
                                         # cbind(parameter = "s1_p", spec = "misspecified", bias = s1pmismod1_bias),
                                         # cbind(parameter = "s2_p", spec = "misspecified", bias = s2pmismod1_bias)))

all_biasandcov_df$bias <- as.numeric(all_biasandcov_df$bias)

# replace parameter names with parse format style for mathematical symbols in figures
all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_b0"),]$parameter  <- "beta[0]^~s1"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_b0"),]$parameter  <- "beta[0]^~s2"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_b1"),]$parameter  <- "beta[1]^~s1"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_b1"),]$parameter  <- "beta[1]^~s2"
all_biasandcov_df[which(all_biasandcov_df$parameter == "f12"),]$parameter  <- "f[12]"
# all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_p"),]$parameter  <- "p^s1"
# all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_p"),]$parameter  <- "p^s2"

# set as factor and set levels
all_biasandcov_df$parameter <- factor(all_biasandcov_df$parameter, levels = c('f[12]', 'beta[0]^~s1', 'beta[0]^~s2', 'beta[1]^~s1', 'beta[1]^~s2'))
all_biasandcov_df$spec <- factor(all_biasandcov_df$spec, levels = c('true', 'misspecified'))

# data summary function for violin plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# custom colour palette
cbPalette <- c( "#E69F00", "#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")


# plot the bias figure using ggplot and save out
tiff("rota_sim_missinghabitatcov_positive_errorplot.tiff", units="in", width=6, height=4, res=300)

ggplot(all_biasandcov_df, aes(x = parameter, y = bias, fill = spec)) + geom_violin() + stat_summary(fun.data=data_summary, size = 0.1, position=position_dodge(width=0.9)) +scale_fill_manual(values=cbPalette, name = 'model specification') +
  theme_classic() + xlab('model parameter')+ geom_hline(data=all_biasandcov_df, aes(yintercept=0))+   scale_x_discrete(labels = scales::parse_format()) +
  theme(legend.position="bottom") +ggtitle("a)")


dev.off()

###### SIMULATION STUDY 2.2 ###### 
############ Shared negative response to habitat covariate ###########
# True parameter values on logit-scale
f1 <- 0 # mean occupancy of species 1
f2 <- 0 # mean occupancy of species 2 
f12 <- 0 # second order interaction term s1 - s2
s1b1 <- -1 # beta_1 (slope term) for s1
s2b1 <- -1 # beta_1 (slope term) for s2
s1a0 <- 0.4075  # mean detection probability of species 1 
s2a0 <- 0.4075 # mean detection probability of species 2 

# Calculate psi to show results of parameter values
psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

# Organize effect sizes
ef <- list(state = c(f1,s1b1, f2, s2b1, f12), det = c(0.4075, 0.4075))

# Blank unmarkedFrame with study design
J <- 1000
K <- 5   # occasions
y <- matrix(NA, J, K)


# create some covariates
site_covs <- data.frame(forest = rnorm(J))

# create unmarked multi dataframe
temp <- unmarkedFrameOccuMulti(list(y,y), siteCovs = site_covs)
str(temp)
head(temp)

# Formulas
sf <- c("~1+forest", "~1+forest", "~1") 
df <- c("~1", "~1")

# Simulate datasets
umf <- simulate(temp, stateformulas=sf, detformulas=df, coefs=ef, nsim=250
str(umf)

# Fit models
truemod2_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=sf, detformulas=df, data=x)
})

# Get CIs for true model (state model)
truemod2_cis <- lapply(1:length(truemod2_fits), function(i){
  data.frame(confint(truemod2_fits[[i]], type='state'), pars = names(coef(truemod2_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get CIs for true model (observation model)
truemod2_det_cis <- lapply(1:length(truemod2_fits), function(i){
  data.frame(confint(truemod2_fits[[i]], type='det'), pars = names(coef(truemod2_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format data
true_out <- do.call(rbind, truemod2_cis)
true_det_out <- do.call(rbind, truemod2_det_cis)
rownames(true_out) <- NULL
rownames(true_det_out) <- NULL
true_out

##### Process results########
# s1_psi/f1
truemod2_f1_est <- sapply(truemod2_fits, function(x) coef(x)[1])
truemod2_f1_CIs <- subset(true_out, pars == "psi([sp1] (Intercept))")
truemod2_s1psi_est <- plogis(truemod2_f1_est)
# bias
f1truemod2_bias = (truemod2_f1_est-f1)/f1
s1psitruemod2_bias = (truemod2_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1truemod2_bias)
mean(s1psitruemod2_bias)
hist(truemod2_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_truemod2_coverage = 1*(truemod2_f1_CIs[2]>=f1&truemod2_f1_CIs[1]<=f1)
f1_truemod2_coverage_table <- table(f1_truemod2_coverage)


# s1b1
truemod2_s1b1_est <- sapply(truemod2_fits, function(x) coef(x)[2])
truemod2_s1b1_CIs <- subset(true_out, pars == "psi([sp1] forest)")
# bias
s1b1truemod2_bias = (truemod2_s1b1_est-s1b1)/s1b1
mean(s1b1truemod2_bias)
hist(truemod2_s1b1_est, main = "s1b1", xlab = "s1b1 est")
abline(v=s1b1, col='red', lwd=2)
# coverage
s1b1_truemod2_coverage = 1*(truemod2_s1b1_CIs[2]>=s1b1&truemod2_s1b1_CIs[1]<=s1b1)
s1b1_truemod2_coverage_table <- table(s1b1_truemod2_coverage)

# f2/s2_psi
truemod2_f2_est <- sapply(truemod2_fits, function(x) coef(x)[3])
truemod2_f2_CIs <- subset(true_out, pars == "psi([sp2] (Intercept))")
truemod2__s2psi_est <- plogis(truemod2_f2_est)
# bias
f2truemod2_bias = (truemod2_f2_est-f2)/f2
# s2_psitruemod2_bias = (truemod2__s2psi_est-(plogis(f2)))/plogis(f2)
mean(f2truemod2_bias)
mean(s2_psitruemod2_bias)
hist(truemod2_f2_est, main="f2", xlab="f2 estimate")
abline(v=f2, col='red', lwd=2)
# coverage
f2_truemod2_coverage = 1*(truemod2_f2_CIs[2]>=f2&truemod2_f2_CIs[1]<=f2)
f2_truemod2_coverage_table <- table(f2_truemod2_coverage)

# s2b1
truemod2_s2b1_est <- sapply(truemod2_fits, function(x) coef(x)[4])
truemod2_s2b1_CIs <- subset(true_out, pars == "psi([sp2] forest)")
# bias
s2b1truemod2_bias = (truemod2_s2b1_est-s2b1)/s2b1
mean(s2b1truemod2_bias)
hist(truemod2_s2b1_est, main="s2b1", xlab="s2b1 est")
abline(v=s2b1, col='red', lwd=2)
# coverage
s2b1_truemod2_coverage = 1*(truemod2_s2b1_CIs[2]>=s2b1&truemod2_s2b1_CIs[1]<=s2b1)
s2b1_truemod2_coverage_table <- table(s2b1_truemod2_coverage)

# interaction term f12
truemod2_f12_est <- sapply(truemod2_fits, function(x) coef(x)[5])
truemod2_f12_CIs <- subset(true_out, pars == "psi([sp1:sp2] (Intercept))")
# bias
f12truemod2_bias = (truemod2_f12_est-f12)
mean(f12truemod2_bias)
hist(truemod2_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_truemod2_coverage = 1*(truemod2_f12_CIs[2]>=f12&truemod2_f12_CIs[1]<=f12)
f12_truemod2_coverage_table <- table(f12_truemod2_coverage)

# p_s1/a0
truemod2_a0_s1_est <- sapply(truemod2_fits, function(x) coef(x)[6])
truemod2_a0_s1_CIs <- subset(true_det_out, pars == "p([sp1] (Intercept))")
truemod2_p_s1_est <- plogis(truemod2_a0_s1_est)
# bias
s1a0truemod2_bias = (truemod2_a0_s1_est-s1a0)/s1a0
s1ptruemod2_bias = (truemod2_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0truemod2_bias)
mean(s1ptruemod2_bias)
hist(truemod2_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=s1a0, col='red', lwd=2)
# coverage
s1a0truemod2_coverage = 1*(truemod2_a0_s1_CIs[2]>=s1a0&truemod2_a0_s1_CIs[1]<s1a0)
s1a0_truemod2_coverage_table <- table(s1a0truemod2_coverage)

# p_s2/a0
truemod2_a0_s2_est <- sapply(truemod2_fits, function(x) coef(x)[7])
truemod2_a0_s2_CIs <- subset(true_det_out, pars == "p([sp2] (Intercept))")
truemod2_p_s2_est <- plogis(truemod2_a0_s2_est)
# bias
s2a0truemod2_bias = (truemod2_a0_s2_est-s2a0)/s2a0
s2ptruemod2_bias = (truemod2_p_s2_est-(plogis(s2a0)))/plogis(s2a0)
mean(s2a0truemod2_bias)
mean(s2ptruemod2_bias)
hist(truemod2_p_s2_est, main="p_s2", xlab="p_s2 estimate")
abline(v=plogis(s2a0), col='red', lwd=2)
# coverage
a0s2truemod2_coverage = 1*(truemod2_a0_s2_CIs[2]>=s2a0&truemod2_a0_s2_CIs[1]<s2a0)
a0s2_truemod2_coverage_table <- table(a0s2truemod2_coverage)

# misspecify model (missing habitat covariate)
modelstateform <- c("~1", "~1", "~1")
str(umf)

# Fit models 
mismod2_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=modelstateform, detformulas=df, data=x)
})

# Get misspecified model CIs (state model)
mismod2_cis <- lapply(1:length(mismod2_fits), function(i){
  data.frame(confint(mismod2_fits[[i]], type='state'), pars = names(coef(mismod2_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get misspecified model CIs (observation model)
mismod2_det_cis <- lapply(1:length(mismod2_fits), function(i){
  data.frame(confint(mismod2_fits[[i]], type='det'), pars = names(coef(mismod2_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format data
mismod2_out <- do.call(rbind, mismod2_cis)
mismod2_det_out <- do.call(rbind, mismod2_det_cis)
rownames(mismod2_out) <- NULL
rownames(mismod2_det_out) <- NULL

#### process results ####
# bias and coverage for each parameter
# f1/s1_psi
mismod2_f1_est <- sapply(mismod2_fits, function(x) coef(x)[1])
mismod2_f1_CIs <- subset(mismod2_out, pars == "psi([sp1] (Intercept))")
mismod2_s1psi_est <- plogis(mismod2_f1_est)
# bias
f1mismod2_bias = (mismod2_f1_est-f1)/f1
s1psimismod2_bias = (mismod2_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1mismod2_bias)
mean(s1psimismod2_bias)
hist(mismod2_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_mismod2_coverage = 1*(mismod2_f1_CIs[2]>=f1&mismod2_f1_CIs[1]<=f1)
f1_mismod2_coverage_table <- table(f1_mismod2_coverage)


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

# interaction term f12
mismod2_f12_est <- sapply(mismod2_fits, function(x) coef(x)[3])
mismod2_f12_CIs <- subset(mismod2_out, pars == "psi([sp1:sp2] (Intercept))")
# bias
f12mismod2_bias = (mismod2_f12_est-f12)
mean(f12mismod2_bias)
hist(mismod2_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_coverage = 1*(mismod2_f12_CIs[2]>=f12&mismod2_f12_CIs[1]<=f12)
mismod2_f12_coverage <- table(f12_coverage)
(mismod2_f12_coverage[2]/250)*100

#p_s1/a0 
mismod2_a0_s1_est <- sapply(mismod2_fits, function(x) coef(x)[4])
mismod2_a0_s1_CIs <- subset(mismod2_det_out, pars == "p([sp1] (Intercept))")
mismod2_p_s1_est <- plogis(mismod2_a0_s1_est)
# bias
s1a0mismod2_bias = (mismod2_a0_s1_est-s1a0)/s1a0
s1pmismod2_bias = (mismod2_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0mismod2_bias)
mean(s1pmismod2_bias)
hist(mismod2_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=plogis(s1a0), col='red', lwd=2)
# coverage
s1a0mismod2_coverage = 1*(mismod2_a0_s1_CIs[2]>=s1a0&mismod2_a0_s1_CIs[1]<s1a0)
s1a0_mismod2_coverage_table <- table(s1a0mismod2_coverage)

# p_s2/a0
mismod2_a0_s2_est <- sapply(mismod2_fits, function(x) coef(x)[5])
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


# combine all bias from all simulations into a df
all_biasandcov_df <- as.data.frame(rbind(cbind(parameter = "s1_b0", spec = "true", bias = f1truemod2_bias),
                                         cbind(parameter = "s2_b0", spec = "true", bias = f2truemod2_bias),
                                         cbind(parameter = "s1_b1", spec = "true", bias = s1b1truemod2_bias),
                                         cbind(parameter = "s2_b1", spec = "true", bias = s2b1truemod2_bias),
                                         cbind(parameter = "f12", spec = "true", bias = f12truemod2_bias),
                                         # cbind(parameter = "s1_p", spec = "true", bias = s1ptruemod2_bias),
                                         # cbind(parameter = "s2_p", spec = "true", bias = s2ptruemod2_bias),
                                         cbind(parameter = "s1_b0", spec = "misspecified", bias = f1mismod2_bias),
                                         cbind(parameter = "s2_b0", spec = "misspecified", bias = f2mismod2_bias),
                                         cbind(parameter = "s1_b1", spec = "misspecified", bias = NA),
                                         cbind(parameter = "s2_b1", spec = "misspecified", bias = NA),
                                         cbind(parameter = "f12", spec = "misspecified", bias = f12mismod2_bias)))

all_biasandcov_df$bias <- as.numeric(all_biasandcov_df$bias)

# replace parameter names with parse format style for mathematical symbols in figures
all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_b0"),]$parameter  <- "beta[0]^~s1"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_b0"),]$parameter  <- "beta[0]^~s2"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_b1"),]$parameter  <- "beta[1]^~s1"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_b1"),]$parameter  <- "beta[1]^~s2"
all_biasandcov_df[which(all_biasandcov_df$parameter == "f12"),]$parameter  <- "f[12]"
# all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_p"),]$parameter  <- "p^s1"
# all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_p"),]$parameter  <- "p^s2"

# set as factor and set levels
all_biasandcov_df$parameter <- factor(all_biasandcov_df$parameter, levels = c('f[12]', 'beta[0]^~s1', 'beta[0]^~s2', 'beta[1]^~s1', 'beta[1]^~s2'))
all_biasandcov_df$spec <- factor(all_biasandcov_df$spec, levels = c('true', 'misspecified'))


# data summary function for violin plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# custom colour palette
cbPalette <- c( "#E69F00", "#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plot figure in ggplot and save out
tiff("rota_sim_missinghabitatcov_bothnegative_errorplot.tiff", units="in", width=6, height=4, res=300)

ggplot(all_biasandcov_df, aes(x = parameter, y = bias, fill = spec)) + geom_violin() + stat_summary(fun.data=data_summary, size = 0.1, position=position_dodge(width=0.9)) +scale_fill_manual(values=cbPalette, name = 'model specification') +
  theme_classic() + xlab('model parameter')+ geom_hline(data=all_biasandcov_df, aes(yintercept=0))+   scale_x_discrete(labels = scales::parse_format()) +
  theme(legend.position="bottom") +ggtitle("b)")


dev.off()


###### SIMULATION STUDY 2.3 ###### 
############ Opposite responses to a shared habitat covariate ###########

# True parameter values on logit-scale
f1 <- 0 # mean occupancy of species 1
f2 <- 0 # mean occupancy of species 2 
f12 <- 0 # second order interaction term s1 - s2
s1b1 <- 1 # beta_1 (slope term) for s1
s2b1 <- -1 # beta_1 (slope term) for s2
s1a0 <- 0.4075  # mean detection probability of species 1 
s2a0 <- 0.4075 # mean detection probability of species 2 

# Calculate psi to show results of parameter values
psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

# Organize effect sizes
ef <- list(state = c(f1,s1b1, f2, s2b1, f12), det = c(0.4075, 0.4075))

# Blank unmarkedFrame with study design
J <- 1000
K <- 5   # occasions
y <- matrix(NA, J, K)

# create some covariates
site_covs <- data.frame(forest = rnorm(J))

# create unmarked multi dataframe
temp <- unmarkedFrameOccuMulti(list(y,y), siteCovs = site_covs)
str(temp)
head(temp)

# Formulas
sf <- c("~1+forest", "~1+forest", "~1") 
df <- c("~1", "~1")

# Simulate datasets
umf <- simulate(temp, stateformulas=sf, detformulas=df, coefs=ef, nsim=250)
str(umf)

# Fit models 
truemod3_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=sf, detformulas=df, data=x)
})

# Get CIs for true model (state model)
truemod3_cis <- lapply(1:length(truemod3_fits), function(i){
  data.frame(confint(truemod3_fits[[i]], type='state'), pars = names(coef(truemod3_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get CIs for true model (observation model)
truemod3_det_cis <- lapply(1:length(truemod3_fits), function(i){
  data.frame(confint(truemod3_fits[[i]], type='det'), pars = names(coef(truemod3_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format the data
true_out <- do.call(rbind, truemod3_cis)
true_det_out <- do.call(rbind, truemod3_det_cis)
rownames(true_out) <- NULL
rownames(true_det_out) <- NULL
true_out

##### Process results########
# s1_psi/f1
truemod3_f1_est <- sapply(truemod3_fits, function(x) coef(x)[1])
truemod3_f1_CIs <- subset(true_out, pars == "psi([sp1] (Intercept))")
truemod3_s1psi_est <- plogis(truemod3_f1_est)
# bias
f1truemod3_bias = (truemod3_f1_est-f1)/f1
s1psitruemod3_bias = (truemod3_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1truemod3_bias)
mean(s1psitruemod3_bias)
hist(truemod3_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_truemod3_coverage = 1*(truemod3_f1_CIs[2]>=f1&truemod3_f1_CIs[1]<=f1)
f1_truemod3_coverage_table <- table(f1_truemod3_coverage)

# s1b1
truemod3_s1b1_est <- sapply(truemod3_fits, function(x) coef(x)[2])
truemod3_s1b1_CIs <- subset(true_out, pars == "psi([sp1] forest)")
# bias
s1b1truemod3_bias = (truemod3_s1b1_est-s1b1)/s1b1
mean(s1b1truemod3_bias)
hist(truemod3_s1b1_est, main = "s1b1", xlab = "s1b1 est")
abline(v=s1b1, col='red', lwd=2)
# coverage
s1b1_truemod3_coverage = 1*(truemod3_s1b1_CIs[2]>=s1b1&truemod3_s1b1_CIs[1]<=s1b1)
s1b1_truemod3_coverage_table <- table(s1b1_truemod3_coverage)

# f2/s2_psi
truemod3_f2_est <- sapply(truemod3_fits, function(x) coef(x)[3])
truemod3_f2_CIs <- subset(true_out, pars == "psi([sp2] (Intercept))")
truemod3__s2psi_est <- plogis(truemod3_f2_est)
# bias
f2truemod3_bias = (truemod3_f2_est-f2)/f2
s2_psitruemod3_bias = (truemod3__s2psi_est-(plogis(f2)))/plogis(f2)
mean(f2truemod3_bias)
mean(s2_psitruemod3_bias)
hist(truemod3_f2_est, main="f2", xlab="f2 estimate")
abline(v=f2, col='red', lwd=2)
# coverage
f2_truemod3_coverage = 1*(truemod3_f2_CIs[2]>=f2&truemod3_f2_CIs[1]<=f2)
f2_truemod3_coverage_table <- table(f2_truemod3_coverage)

# s2b1
truemod3_s2b1_est <- sapply(truemod3_fits, function(x) coef(x)[4])
truemod3_s2b1_CIs <- subset(true_out, pars == "psi([sp2] forest)")
# bias
s2b1truemod3_bias = (truemod3_s2b1_est-s2b1)/s2b1
mean(s2b1truemod3_bias)
hist(truemod3_s2b1_est, main="s2b1", xlab="s2b1 est")
abline(v=s2b1, col='red', lwd=2)
# coverage
s2b1_truemod3_coverage = 1*(truemod3_s2b1_CIs[2]>=s2b1&truemod3_s2b1_CIs[1]<=s2b1)
s2b1_truemod3_coverage_table <- table(s2b1_truemod3_coverage)

# interaction term f12
truemod3_f12_est <- sapply(truemod3_fits, function(x) coef(x)[5])
truemod3_f12_CIs <- subset(true_out, pars == "psi([sp1:sp2] (Intercept))")
# bias
f12truemod3_bias = (truemod3_f12_est-f12)
mean(f12truemod3_bias)
hist(truemod3_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_truemod3_coverage = 1*(truemod3_f12_CIs[2]>=f12&truemod3_f12_CIs[1]<=f12)
f12_truemod3_coverage_table <- table(f12_truemod3_coverage)

# p_s1/a0
truemod3_a0_s1_est <- sapply(truemod3_fits, function(x) coef(x)[6])
truemod3_a0_s1_CIs <- subset(true_det_out, pars == "p([sp1] (Intercept))")
truemod3_p_s1_est <- plogis(truemod3_a0_s1_est)
# bias
s1a0truemod3_bias = (truemod3_a0_s1_est-s1a0)/s1a0
s1ptruemod3_bias = (truemod3_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0truemod3_bias)
mean(s1ptruemod3_bias)
hist(truemod3_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=s1a0, col='red', lwd=2)
# coverage
s1a0truemod3_coverage = 1*(truemod3_a0_s1_CIs[2]>=s1a0&truemod3_a0_s1_CIs[1]<s1a0)
s1a0_truemod3_coverage_table <- table(s1a0truemod3_coverage)

# p_s2/a0
truemod3_a0_s2_est <- sapply(truemod3_fits, function(x) coef(x)[7])
truemod3_a0_s2_CIs <- subset(true_det_out, pars == "p([sp2] (Intercept))")
truemod3_p_s2_est <- plogis(truemod3_a0_s2_est)
# bias
s2a0truemod3_bias = (truemod3_a0_s2_est-s2a0)/s2a0
s2ptruemod3_bias = (truemod3_p_s2_est-(plogis(s2a0)))/plogis(s2a0)
mean(s2a0truemod3_bias)
mean(s2ptruemod3_bias)
hist(truemod3_p_s2_est, main="p_s2", xlab="p_s2 estimate")
abline(v=plogis(s2a0), col='red', lwd=2)
# coverage
a0s2truemod3_coverage = 1*(truemod3_a0_s2_CIs[2]>=s2a0&truemod3_a0_s2_CIs[1]<s2a0)
a0s2_truemod3_coverage_table <- table(a0s2truemod3_coverage)

# misspecify model (missing habitat covariate)
modelstateform <- c("~1", "~1", "~1")

str(umf)
# Fit models to check
mismod3_fits <- pbapply::pblapply(umf, function(x){
  occuMulti(stateformulas=modelstateform, detformulas=df, data=x)
})

# Get misspecified model CIs (state model)
mismod3_cis <- lapply(1:length(mismod3_fits), function(i){
  data.frame(confint(mismod3_fits[[i]], type='state'), pars = names(coef(mismod3_fits[[i]], type='state')),
             rep=i, check.names=FALSE)
})

# Get misspecified model CIs (observation model)
mismod3_det_cis <- lapply(1:length(mismod3_fits), function(i){
  data.frame(confint(mismod3_fits[[i]], type='det'), pars = names(coef(mismod3_fits[[i]], type='det')),
             rep=i, check.names=FALSE)
})

# format data
mismod3_out <- do.call(rbind, mismod3_cis)
mismod3_det_out <- do.call(rbind, mismod3_det_cis)
rownames(mismod3_out) <- NULL
rownames(mismod3_det_out) <- NULL

###### process results #######
# bias and coverage for each parameter
# f1/s1_psi
mismod3_f1_est <- sapply(mismod3_fits, function(x) coef(x)[1])
mismod3_f1_CIs <- subset(mismod3_out, pars == "psi([sp1] (Intercept))")
mismod3_s1psi_est <- plogis(mismod3_f1_est)
# bias
f1mismod3_bias = (mismod3_f1_est-f1)/f1
s1psimismod3_bias = (mismod3_s1psi_est-(plogis(f1)))/plogis(f1)
mean(f1mismod3_bias)
mean(s1psimismod3_bias)
hist(mismod3_f1_est, main="f1", xlab="f1 estimate")
abline(v=f1, col='red', lwd=2)
# coverage
f1_mismod3_coverage = 1*(mismod3_f1_CIs[2]>=f1&mismod3_f1_CIs[1]<=f1)
f1_mismod3_coverage_table <- table(f1_mismod3_coverage)

# f2
mismod3_f2_est <- sapply(mismod3_fits, function(x) coef(x)[3])
mismod3_f2_CIs <- subset(mismod3_out, pars == "psi([sp2] (Intercept))")
mismod3_s2psi_est <- plogis(mismod3_f2_est)
# bias
f2mismod3_bias = (mismod3_f2_est-f2)/f2
s2psimismod3_bias = (mismod3_s2psi_est-(plogis(f2)))/plogis(f2)
# bias
mean(f2mismod3_bias)
mean(s2psimismod3_bias)
hist(mismod3_f2_est, main="f2", xlab="f2 estimate")
abline(v=f2, col='red', lwd=2)
# coverage
f2_mismod3_coverage = 1*(mismod3_f2_CIs[2]>=f2&mismod3_f2_CIs[1]<=f2)
f2_mismod3_coverage_table <- table(f2_mismod3_coverage)

# interaction term f12
mismod3_f12_est <- sapply(mismod3_fits, function(x) coef(x)[3])
mismod3_f12_CIs <- subset(mismod3_out, pars == "psi([sp1:sp2] (Intercept))")
# bias
f12mismod3_bias = (mismod3_f12_est-f12)
mean(f12mismod3_bias)
hist(mismod3_f12_est, main="f12", xlab="f12 estimate")
abline(v=f12, col='red', lwd=2)
# coverage
f12_coverage = 1*(mismod3_f12_CIs[2]>=f12&mismod3_f12_CIs[1]<=f12)
mismod3_f12_coverage <- table(f12_coverage)
(mismod3_f12_coverage[2]/250)*100

# p_s1/a0
mismod3_a0_s1_est <- sapply(mismod3_fits, function(x) coef(x)[4])
mismod3_a0_s1_CIs <- subset(mismod3_det_out, pars == "p([sp1] (Intercept))")
mismod3_p_s1_est <- plogis(mismod3_a0_s1_est)
# bias
s1a0mismod3_bias = (mismod3_a0_s1_est-s1a0)/s1a0
s1pmismod3_bias = (mismod3_p_s1_est-(plogis(s1a0)))/plogis(s1a0)
mean(s1a0mismod3_bias)
mean(s1pmismod3_bias)
hist(mismod3_p_s1_est, main="p_s1", xlab="p_s1 estimate")
abline(v=plogis(s1a0), col='red', lwd=2)
# coverage
s1a0mismod3_coverage = 1*(mismod3_a0_s1_CIs[2]>=s1a0&mismod3_a0_s1_CIs[1]<s1a0)
s1a0_mismod3_coverage_table <- table(s1a0mismod3_coverage)

# p_s2/a0
mismod3_a0_s2_est <- sapply(mismod3_fits, function(x) coef(x)[5])
mismod3_a0_s2_CIs <- subset(mismod3_det_out, pars == "p([sp2] (Intercept))")
mismod3_p_s2_est <- plogis(mismod3_a0_s2_est)
# bias
s2a0mismod3_bias = (mismod3_a0_s2_est-s2a0)/s2a0
s2pmismod3_bias = (mismod3_p_s2_est-(plogis(s2a0)))/plogis(s2a0)
mean(s2a0mismod3_bias)
mean(s2pmismod3_bias)

hist(mismod3_p_s2_est, main="p_s2", xlab="p_s2 estimate")
abline(v=plogis(s2a0), col='red', lwd=2)
# coverage
a0s2mismod3_coverage = 1*(mismod3_a0_s2_CIs[2]>=s2a0&mismod3_a0_s2_CIs[1]<s2a0)
a0s2_mismod3_coverage_table <- table(a0s2mismod3_coverage)


# combine simulation results into a df
all_biasandcov_df <- as.data.frame(rbind(cbind(parameter = "s1_psi", spec = "true", bias = f1truemod3_bias),
                                         cbind(parameter = "s2_psi", spec = "true", bias = f2truemod3_bias),
                                         cbind(parameter = "s1_b1", spec = "true", bias = s1b1truemod3_bias),
                                         cbind(parameter = "s2_b1", spec = "true", bias = s2b1truemod3_bias),
                                         cbind(parameter = "f12", spec = "true", bias = f12truemod3_bias),
                                         # cbind(parameter = "s1_p", spec = "true", bias = s1ptruemod3_bias),
                                         # cbind(parameter = "s2_p", spec = "true", bias = s2ptruemod3_bias),
                                         cbind(parameter = "s1_psi", spec = "misspecified", bias = f1mismod3_bias),
                                         cbind(parameter = "s2_psi", spec = "misspecified", bias = f2mismod3_bias),
                                         cbind(parameter = "s1_b1", spec = "misspecified", bias = NA),
                                         cbind(parameter = "s2_b1", spec = "misspecified", bias = NA),
                                         cbind(parameter = "f12", spec = "misspecified", bias = f12mismod3_bias)))
# cbind(parameter = "s1_p", spec = "misspecified", bias = s1pmismod3_bias), 
# cbind(parameter = "s2_p", spec = "misspecified", bias = s2pmismod3_bias)))

all_biasandcov_df$bias <- as.numeric(all_biasandcov_df$bias)

# replace parameter names with parse format style for mathematical symbols in figures
all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_psi"),]$parameter  <- "beta[0]^~s1"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_psi"),]$parameter  <- "beta[0]^~s2"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_b1"),]$parameter  <- "beta[1]^~s1"
all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_b1"),]$parameter  <- "beta[1]^~s2"
all_biasandcov_df[which(all_biasandcov_df$parameter == "f12"),]$parameter  <- "f[12]"
# all_biasandcov_df[which(all_biasandcov_df$parameter == "s1_p"),]$parameter  <- "p^s1"
# all_biasandcov_df[which(all_biasandcov_df$parameter == "s2_p"),]$parameter  <- "p^s2"

# set as factor and set levels
all_biasandcov_df$parameter <- factor(all_biasandcov_df$parameter, levels = c('f[12]', 'beta[0]^~s1', 'beta[0]^~s2', 'beta[1]^~s1', 'beta[1]^~s2'))
all_biasandcov_df$spec <- factor(all_biasandcov_df$spec, levels = c('true', 'misspecified'))

# data summary function for violin plots
data_summary <- function(x) {
  m <- mean(x)
  ymin <- m-sd(x)
  ymax <- m+sd(x)
  return(c(y=m,ymin=ymin,ymax=ymax))
}

# custom colour palette
cbPalette <- c( "#E69F00", "#999999","#56B4E9", "#009E73", "#F0E442", "#0072B2", "#D55E00", "#CC79A7")

# plot in ggplot and save out
tiff("rota_sim_missinghabitatcov_oppositesigns_errorplot.tiff", units="in", width=6, height=4, res=300)

ggplot(all_biasandcov_df, aes(x = parameter, y = bias, fill = spec)) + geom_violin() + stat_summary(fun.data=data_summary, size = 0.1, position=position_dodge(width=0.9)) +scale_fill_manual(values=cbPalette, name = 'model specification') +
  theme_classic() + xlab('model parameter')+ geom_hline(data=all_biasandcov_df, aes(yintercept=0))+   scale_x_discrete(labels = scales::parse_format()) +
  theme(legend.position="bottom") +ggtitle("c)")


dev.off()


# save.image("occuMulti_missinghabitatcovs_allsimsrun.RData")

load("occuMulti_missinghabitatcovs_allsimsrun.RData")
