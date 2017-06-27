# micro_ecol
Scripts used in the analysis of amplicon microbiota data in R

## Pairwise_adonis_all.r
Function to run pairwise PERMANOVAs of all levels of a variable of interest across beta-diversity measures. Adapted from Pedro Martinez Arbizu (https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R). Allows use of Bray-Curtis, Jaccard and UniFrac. Allows `strata` function within adonis
