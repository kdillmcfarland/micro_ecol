# micro_ecol
Scripts used in the analysis of amplicon microbiota data

## clean_repFasta_FAST.py

Renames fasta sequences output from mothur from automatically generated name to OTU #

## Pairwise_adonis

Functions to run pairwise PERMANOVAs with multiple comparison (fdr) correction across all levels variables of interest across beta-diversity measures. Adapted from Pedro Martinez Arbizu (https://www.researchgate.net/post/How_can_I_do_PerMANOVA_pairwise_contrasts_in_R). 

* Multiple variables of interest allowed. 
* Bray-Curtis, Jaccard, or UniFrac allowed. 
* Stratification of model allowed with strata option.
* Output saved to object in R environment. 

## pairwise_permdisp.R

Permutational tests of beta-dispersion (PERMDISP) across multiple variables of interest and pairwise with TukeyHSD. Multiple variables of interest allowed. Bray-Curtis or Jaccard only.

## R_krusk.R

Automates Kruskal-Wallis test execution on OTUs of interest as determined by simper_pretty.R. See https://github.com/asteinberger9 for up to date version

## simper_pretty.R

Automates SIMAPERA execution for comparisons OTUs across variables of interest. See https://github.com/asteinberger9 for up to date version