# Read in information on included papers that use Rota et al. 2016
incl <- read.csv("included_papers.csv")

nrow(incl)

# Cameras
sum(incl$cameras, na.rm=TRUE)
mean(incl$cameras, na.rm=TRUE)

# Number of species
median(incl$max_species, na.rm=TRUE)
range(incl$max_species, na.rm=TRUE)
mean(incl$max_species < 3, na.rm=TRUE)

# Distribution of number of sites
hist(incl$max_sites)
hist(incl$max_sites[incl$max_sites < 1000])

# Number of sites
median(incl$max_sites, na.rm=TRUE)
range(incl$max_sites, na.rm=TRUE)
mean(incl$max_sites <= 1000, na.rm=TRUE)

# Number of occasions
median(incl$max_occasions, na.rm=TRUE)
range(incl$max_occasions, na.rm=TRUE)
mean(incl$max_occasions < 30, na.rm=TRUE)

# Occupancy estimates
mean(incl$psi_min, na.rm=TRUE)
mean(incl$psi_max, na.rm=TRUE)

# Detection probability estimates
mean(incl$p_min, na.rm=TRUE)
mean(incl$p_max, na.rm=TRUE)

# Interaction order
mean(incl$max_order < 3, na.rm=TRUE)

# Papers that have covariates on interactions
mean(incl$int_covs, na.rm=TRUE)

# Significant interactions
sum(!is.na(incl$n_int_sig))
mean(incl$n_int_sig > 0, na.rm=TRUE)

# Size of interactions
sum(!is.na(incl$int_min))
max_abs_int <- cbind(abs(incl$int_min), abs(incl$int_max))
max_abs_int <- apply(max_abs_int, 1, max)
median(max_abs_int, na.rm=TRUE)
range(max_abs_int, na.rm=TRUE)

# Estimate size of mean interaction
# on occupancy scale
f1 <- -0.55
f2 <- -0.55
f12 <- 1.7 # <- EFFECT SIZE OF INTERACTION

psi <- numeric(4)
psi[1] <- exp(f1 + f2 + f12)
psi[2] <- exp(f1)
psi[3] <- exp(f2)
psi[4] <- 1
psi <- psi / sum(psi)
names(psi) <- c("[11]","[10]","[01]","[00]")
round(psi, 3)

# Marginal occupancies of species 1 and 2
psi[1] + psi[2] # 60
psi[1] + psi[3] # 60

# Occupancy of species 1 | species 2 present
psi[1] / (psi[1] + psi[3]) # 76
# Occupancy of species 1 | species 2 absent
psi[2] / (psi[2] + psi[4]) # 37
