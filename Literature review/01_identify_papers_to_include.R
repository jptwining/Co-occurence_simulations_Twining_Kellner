# Read in list of all papers citing Rota et al. 2016
papers <- read.csv("All_Citing_Papers.csv")

nrow(papers)
sum(papers$Uses_Model, na.rm=TRUE)

# Identify papers to keep because they actually use the model
keep <- !is.na(papers$Uses_Model) & papers$Uses_Model == 1
out <- papers[keep, c("DOI", "Source.Title", "Publication.Year", "Notes")]

# Format rows/columns
names(out) <- c("DOI", "Journal", "Year", "Notes")
rownames(out) <- NULL

# New data to collect
out$nsites <- NA
out$nspecies <- NA
out$noccasions <- NA
out$psi_min <- NA
out$psi_max <- NA
out$p_min <- NA
out$p_max <- NA
out$int_min <- NA
out$int_max <- NA
out$n_int_sig <- NA
out$n_int_covs <- NA
out$n_sites_detected <- NA
out$n_overlap_sites <- NA

# Write output
write.csv(out, "included_papers.csv", row.names=FALSE)
