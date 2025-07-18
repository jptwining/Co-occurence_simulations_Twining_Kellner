library(unmarked)

# Based on literature review, we decided to use the following settings
# in our simulations
# species: 2, 3, 5
# sites  : 50, 100, 200, 500, 1000
# occasions: 5, 10, 30
# psi: 0.2, 0.6
# p  : 0.2, 0.6
# covariates yes / no
# big and small effect size (=absolute value of interaction term)

# 2 species effect sizes-------------------------------------------------------

# High occupancy---------------------------------------------------------------
# ~0.6 marginal occupancy
# Big effect size of interaction
f1 <- -0.2 
f2 <- -0.2
f12 <- 1

psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

psi[1] + psi[2] # 59
psi[1] + psi[3] # 59

psi[1] / (psi[1] + psi[3]) # 68
psi[2] / (psi[2] + psi[4]) # 45

highocc_2sp_bigeffect <- c(-0.2, -0.2, 1)

# ~0.6 marginal occupancy
# Small effect size of interaction
f1 <- 0.2
f2 <- 0.2
f12 <- 0.4 # <- EFFECT SIZE OF INTERACTION

psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

psi[1] + psi[2] # 61
psi[1] + psi[3] # 61

psi[1] / (psi[1] + psi[3]) # 65
psi[2] / (psi[2] + psi[4]) # 55


highocc_2sp_smalleffect <- c(0.2, 0.2, 0.4)

# No effect size of interaction
f1 <- 0.4
f2 <- 0.4
f12 <- 0 # <- EFFECT SIZE OF INTERACTION

psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

psi[1] + psi[2] # 61
psi[1] + psi[3] # 61

psi[1] / (psi[1] + psi[3]) # 65
psi[2] / (psi[2] + psi[4]) # 55

highocc_2sp_noeffect <- c(0.4, 0.4, 0)


# Low occupancy----------------------------------------------------------------

# ~0.2 marginal occupancy
# Big effect size of interaction
f1 <- -1.7
f2 <- -1.7
f12 <- 1.2

psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

psi[1] + psi[2] # 20
psi[1] + psi[3] # 20

psi[1] / (psi[1] + psi[3]) # 37
psi[2] / (psi[2] + psi[4]) # 15

lowocc_2sp_bigeffect <- c(-1.7, -1.7, 1.2)

# ~0.2 marginal occupancy
# Small effect size of interaction
f1 <- -1.5
f2 <- -1.5
f12 <- 0.6

psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

psi[1] + psi[2] # 20
psi[1] + psi[3] # 20

psi[1] / (psi[1] + psi[3]) # 28
psi[2] / (psi[2] + psi[4]) # 18

lowocc_2sp_smalleffect <- c(-1.5,-1.5,0.6)

# ~0.2 marginal occupancy
# No effect size of interaction
f1 <- -1.4
f2 <- -1.4
f12 <- 0

psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

psi[1] + psi[2] # 20
psi[1] + psi[3] # 20

psi[1] / (psi[1] + psi[3]) # 20
psi[2] / (psi[2] + psi[4]) # 20

lowocc_2sp_noeffect <- c(-1.4,-1.4,0)


# 3 species effect sizes-------------------------------------------------------

# High occupancy---------------------------------------------------------------
# ~0.6 marginal occupancy
# Big effect size of interaction
f1 <- -0.6
f2 <- -0.6
f3 <- -0.6
f12 <- 0.9
f13 <- 0.9
f23 <- 0.9
# no 3-way interaction

psi <- numeric(8)
psi[1] <- exp(f1 + f2 + f3 + f12 + f13 + f23)
psi[2] <- exp(f1 + f2 + f12)
psi[3] <- exp(f1 + f3 + f13)
psi[4] <- exp(f1)
psi[5] <- exp(f2 + f3 + f23)
psi[6] <- exp(f2)
psi[7] <- exp(f3)
psi[8] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[111]","[110]","[101]","[100]","[011]","[010]","[001]","[000]")
round(psi, 3)

psi[1] + psi[2] + psi[3] + psi[4] # 61
psi[1] + psi[2] + psi[5] + psi[6] # 61
psi[1] + psi[3] + psi[5] + psi[7] # 61

(psi[1] + psi[2]) / (psi[1] + psi[2] + psi[5] + psi[6]) # 71
(psi[3] + psi[4]) / (psi[3] + psi[4] + psi[7] + psi[8]) # 45

highocc_3sp_bigeffect <- c(-0.6,-0.6,-0.6,0.9,0.9,0.9)


# Small effect size of interaction
f1 <- -0.1
f2 <- -0.1
f3 <- -0.1
f12 <- 0.4
f13 <- 0.4
f23 <- 0.4
# no 3-way interaction

psi <- numeric(8)
psi[1] <- exp(f1 + f2 + f3 + f12 + f13 + f23)
psi[2] <- exp(f1 + f2 + f12)
psi[3] <- exp(f1 + f3 + f13)
psi[4] <- exp(f1)
psi[5] <- exp(f2 + f3 + f23)
psi[6] <- exp(f2)
psi[7] <- exp(f3)
psi[8] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[111]","[110]","[101]","[100]","[011]","[010]","[001]","[000]")
round(psi, 3)

psi[1] + psi[2] + psi[3] + psi[4] # 59
psi[1] + psi[2] + psi[5] + psi[6] # 59
psi[1] + psi[3] + psi[5] + psi[7] # 59

(psi[1] + psi[2]) / (psi[1] + psi[2] + psi[5] + psi[6]) # 63
(psi[3] + psi[4]) / (psi[3] + psi[4] + psi[7] + psi[8]) # 53

highocc_3sp_smalleffect <- c(-0.1,-0.1,-0.1,0.4,0.4,0.4)


# No interaction
f1 <- 0.4
f2 <- 0.4
f3 <- 0.4
f12 <- 0
f13 <- 0
f23 <- 0
# no 3-way interaction

psi <- numeric(8)
psi[1] <- exp(f1 + f2 + f3 + f12 + f13 + f23)
psi[2] <- exp(f1 + f2 + f12)
psi[3] <- exp(f1 + f3 + f13)
psi[4] <- exp(f1)
psi[5] <- exp(f2 + f3 + f23)
psi[6] <- exp(f2)
psi[7] <- exp(f3)
psi[8] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[111]","[110]","[101]","[100]","[011]","[010]","[001]","[000]")
round(psi, 3)

psi[1] + psi[2] + psi[3] + psi[4] # 60
psi[1] + psi[2] + psi[5] + psi[6] # 60
psi[1] + psi[3] + psi[5] + psi[7] # 60

(psi[1] + psi[2]) / (psi[1] + psi[2] + psi[5] + psi[6]) # 65
(psi[3] + psi[4]) / (psi[3] + psi[4] + psi[7] + psi[8]) # 53

highocc_3sp_noeffect <- c(0.4,0.4,0.4,0,0,0)


# Low occupancy---------------------------------------------------------------
# ~0.2 marginal occupancy
# Big effect size of interaction
f1 <- -1.9
f2 <- -1.9
f3 <- -1.9
f12 <- 1 
f13 <- 1
f23 <- 1
# no 3-way interaction

psi <- numeric(8)
psi[1] <- exp(f1 + f2 + f3 + f12 + f13 + f23)
psi[2] <- exp(f1 + f2 + f12)
psi[3] <- exp(f1 + f3 + f13)
psi[4] <- exp(f1)
psi[5] <- exp(f2 + f3 + f23)
psi[6] <- exp(f2)
psi[7] <- exp(f3)
psi[8] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[111]","[110]","[101]","[100]","[011]","[010]","[001]","[000]")
round(psi, 3)

psi[1] + psi[2] + psi[3] + psi[4] # 20
psi[1] + psi[2] + psi[5] + psi[6] # 20
psi[1] + psi[3] + psi[5] + psi[7] # 20

(psi[1] + psi[2]) / (psi[1] + psi[2] + psi[5] + psi[6]) # 38
(psi[3] + psi[4]) / (psi[3] + psi[4] + psi[7] + psi[8]) # 15

lowocc_3sp_bigeffect <- c(-1.9,-1.9,-1.9,1,1,1)

# Small effect size of interaction
f1 <- -1.6
f2 <- -1.6
f3 <- -1.6
f12 <- 0.5
f13 <- 0.5
f23 <- 0.5
# no 3-way interaction

psi <- numeric(8)
psi[1] <- exp(f1 + f2 + f3 + f12 + f13 + f23)
psi[2] <- exp(f1 + f2 + f12)
psi[3] <- exp(f1 + f3 + f13)
psi[4] <- exp(f1)
psi[5] <- exp(f2 + f3 + f23)
psi[6] <- exp(f2)
psi[7] <- exp(f3)
psi[8] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[111]","[110]","[101]","[100]","[011]","[010]","[001]","[000]")
round(psi, 3)

psi[1] + psi[2] + psi[3] + psi[4] # 20
psi[1] + psi[2] + psi[5] + psi[6] # 20
psi[1] + psi[3] + psi[5] + psi[7] # 20

(psi[1] + psi[2]) / (psi[1] + psi[2] + psi[5] + psi[6]) # 28
(psi[3] + psi[4]) / (psi[3] + psi[4] + psi[7] + psi[8]) # 18

lowocc_3sp_smalleffect <- c(-1.6,-1.6,-1.6,0.5,0.5,0.5)

# No interaction
f1 <- -1.4
f2 <- -1.4
f3 <- -1.4
f12 <- 0
f13 <- 0
f23 <- 0
# no 3-way interaction

psi <- numeric(8)
psi[1] <- exp(f1 + f2 + f3 + f12 + f13 + f23)
psi[2] <- exp(f1 + f2 + f12)
psi[3] <- exp(f1 + f3 + f13)
psi[4] <- exp(f1)
psi[5] <- exp(f2 + f3 + f23)
psi[6] <- exp(f2)
psi[7] <- exp(f3)
psi[8] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[111]","[110]","[101]","[100]","[011]","[010]","[001]","[000]")
round(psi, 3)

psi[1] + psi[2] + psi[3] + psi[4] # 20
psi[1] + psi[2] + psi[5] + psi[6] # 20
psi[1] + psi[3] + psi[5] + psi[7] # 20

(psi[1] + psi[2]) / (psi[1] + psi[2] + psi[5] + psi[6]) # 28
(psi[3] + psi[4]) / (psi[3] + psi[4] + psi[7] + psi[8]) # 18

lowocc_3sp_noeffect <- c(-1.4,-1.4,-1.4,0,0,0)


# 5 species effect sizes-------------------------------------------------------

# High occupancy---------------------------------------------------------------
# ~0.6 marginal occupancy
# Big effect size of interaction
# no 3/4/5-way interaction
f <- c(rep(-1, 5), rep(0.6, 10))

ymat <- matrix(0, 5, 3)
ylist <- replicate(5, ymat, simplify=FALSE)
umf <- unmarkedFrameOccuMulti(ylist)
fdesign <- umf@fDesign[,1:15]

psi <- numeric(32)
for (i in 1:length(psi)){
  psi[i] <- exp(sum(f[fdesign[i,]==1]))
}
psi <- psi / sum(psi)

names(psi) <- gsub("psi","",rownames(fdesign), fixed=TRUE)


# species 1 present
sum(psi[substr(names(psi), 2, 2) == "1"]) # 60
sum(psi[substr(names(psi), 3, 3) == "1"]) # 60

# species 1 present conditional on 2 present
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "1"]) /
             sum(psi[substr(names(psi), 3, 3) == "1"]) #69
# species 1 present conditional on 2 absent
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "0"]) /
             sum(psi[substr(names(psi), 3, 3) == "0"]) #46

highocc_5sp_bigeffect <- f


# Small effect size of interaction
# no 3/4/5-way interaction
f <- c(rep(-0.3, 5), rep(0.3, 10))

ymat <- matrix(0, 5, 3)
ylist <- replicate(5, ymat, simplify=FALSE)
umf <- unmarkedFrameOccuMulti(ylist)
fdesign <- umf@fDesign[,1:15]

psi <- numeric(32)
for (i in 1:length(psi)){
  psi[i] <- exp(sum(f[fdesign[i,]==1]))
}
psi <- psi / sum(psi)

names(psi) <- gsub("psi","",rownames(fdesign), fixed=TRUE)


# species 1 present
sum(psi[substr(names(psi), 2, 2) == "1"]) # 60
sum(psi[substr(names(psi), 3, 3) == "1"]) # 60

# species 1 present conditional on 2 present
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "1"]) /
             sum(psi[substr(names(psi), 3, 3) == "1"]) #64
# species 1 present conditional on 2 absent
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "0"]) /
             sum(psi[substr(names(psi), 3, 3) == "0"]) #55

highocc_5sp_smalleffect <- f

# No interaction
# no 3/4/5-way interaction
f <- c(rep(0.4, 5), rep(0, 10))

ymat <- matrix(0, 5, 3)
ylist <- replicate(5, ymat, simplify=FALSE)
umf <- unmarkedFrameOccuMulti(ylist)
fdesign <- umf@fDesign[,1:15]

psi <- numeric(32)
for (i in 1:length(psi)){
  psi[i] <- exp(sum(f[fdesign[i,]==1]))
}
psi <- psi / sum(psi)

names(psi) <- gsub("psi","",rownames(fdesign), fixed=TRUE)


# species 1 present
sum(psi[substr(names(psi), 2, 2) == "1"]) # 60
sum(psi[substr(names(psi), 3, 3) == "1"]) # 60

# species 1 present conditional on 2 present
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "1"]) /
             sum(psi[substr(names(psi), 3, 3) == "1"]) #60
# species 1 present conditional on 2 absent
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "0"]) /
             sum(psi[substr(names(psi), 3, 3) == "0"]) #60

highocc_5sp_noeffect <- f


# Low occupancy----------------------------------------------------------------
# ~0.2 marginal occupancy
# Big effect size of interaction
# no 3/4/5-way interaction
f <- c(rep(-2.2, 5), rep(0.8, 10))

ymat <- matrix(0, 5, 3)
ylist <- replicate(5, ymat, simplify=FALSE)
umf <- unmarkedFrameOccuMulti(ylist)
fdesign <- umf@fDesign[,1:15]

psi <- numeric(32)
for (i in 1:length(psi)){
  psi[i] <- exp(sum(f[fdesign[i,]==1]))
}
psi <- psi / sum(psi)

names(psi) <- gsub("psi","",rownames(fdesign), fixed=TRUE)


# species 1 present
sum(psi[substr(names(psi), 2, 2) == "1"]) # 20
sum(psi[substr(names(psi), 3, 3) == "1"]) # 20

# species 1 present conditional on 2 present
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "1"]) /
             sum(psi[substr(names(psi), 3, 3) == "1"]) #41
# species 1 present conditional on 2 absent
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "0"]) /
             sum(psi[substr(names(psi), 3, 3) == "0"]) #15

lowocc_5sp_bigeffect <- f


# Small effect size of interaction
# no 3/4/5-way interaction
f <- c(rep(-1.7, 5), rep(0.4, 10))

ymat <- matrix(0, 5, 3)
ylist <- replicate(5, ymat, simplify=FALSE)
umf <- unmarkedFrameOccuMulti(ylist)
fdesign <- umf@fDesign[,1:15]

psi <- numeric(32)
for (i in 1:length(psi)){
  psi[i] <- exp(sum(f[fdesign[i,]==1]))
}
psi <- psi / sum(psi)

names(psi) <- gsub("psi","",rownames(fdesign), fixed=TRUE)


# species 1 present
sum(psi[substr(names(psi), 2, 2) == "1"]) # 21
sum(psi[substr(names(psi), 3, 3) == "1"]) # 21

# species 1 present conditional on 2 present
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "1"]) /
             sum(psi[substr(names(psi), 3, 3) == "1"]) #28
# species 1 present conditional on 2 absent
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "0"]) /
             sum(psi[substr(names(psi), 3, 3) == "0"]) #19

lowocc_5sp_smalleffect <- f

# No effect size of interaction
# no 3/4/5-way interaction
f <- c(rep(-1.4, 5), rep(0, 10))

ymat <- matrix(0, 5, 3)
ylist <- replicate(5, ymat, simplify=FALSE)
umf <- unmarkedFrameOccuMulti(ylist)
fdesign <- umf@fDesign[,1:15]

psi <- numeric(32)
for (i in 1:length(psi)){
  psi[i] <- exp(sum(f[fdesign[i,]==1]))
}
psi <- psi / sum(psi)

names(psi) <- gsub("psi","",rownames(fdesign), fixed=TRUE)


# species 1 present
sum(psi[substr(names(psi), 2, 2) == "1"]) # 20
sum(psi[substr(names(psi), 3, 3) == "1"]) # 20

# species 1 present conditional on 2 present
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "1"]) /
             sum(psi[substr(names(psi), 3, 3) == "1"]) #20
# species 1 present conditional on 2 absent
sum(psi[substr(names(psi), 2, 2) == "1" & substr(names(psi), 3, 3) == "0"]) /
             sum(psi[substr(names(psi), 3, 3) == "0"]) #20

lowocc_5sp_noeffect <- f


# Scenarios--------------------------------------------------------------------

state_scenarios <- list(
  highocc_2sp_bigeffect = highocc_2sp_bigeffect,
  highocc_2sp_smalleffect = highocc_2sp_smalleffect,
  highocc_2sp_noeffect = highocc_2sp_noeffect,

  lowocc_2sp_bigeffect = lowocc_2sp_bigeffect,
  lowocc_2sp_smalleffect = lowocc_2sp_smalleffect,
  lowocc_2sp_noeffect = lowocc_2sp_noeffect,

  highocc_3sp_bigeffect = highocc_3sp_bigeffect,
  highocc_3sp_smalleffect = highocc_3sp_smalleffect,
  highocc_3sp_noeffect = highocc_3sp_noeffect,

  lowocc_3sp_bigeffect = lowocc_3sp_bigeffect,
  lowocc_3sp_smalleffect = lowocc_3sp_smalleffect,
  lowocc_3sp_noeffect = lowocc_3sp_noeffect,

  highocc_5sp_bigeffect = highocc_5sp_bigeffect,
  highocc_5sp_smalleffect = highocc_5sp_smalleffect,
  highocc_5sp_noeffect = highocc_5sp_noeffect,

  lowocc_5sp_bigeffect = lowocc_5sp_bigeffect,
  lowocc_5sp_smalleffect = lowocc_5sp_smalleffect,
  lowocc_5sp_noeffect = lowocc_5sp_noeffect
)

det_scenarios <- list(highdet = qlogis(0.6), lowdet = qlogis(0.2), tinydet=qlogis(0.05))

# Save scenarios for use in power simulations
saveRDS(state_scenarios, "state_scenarios.Rds")
saveRDS(det_scenarios, "det_scenarios.Rds")
