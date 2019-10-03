"Pairwise ADONIS: Pairwise permutational ANOVA (PERMANOVA, adonis)

Options for 1 or more variables of interest, stratification for repeated measures, and beta measure in Bray-Curtis, Jaccard, and UniFrac

DOES NOT SUPPORT interaction terms within models (e.g. a*b). If interaction terms are desired, create a new variable in the metadata such as `paste(a,b, sep='_')`

#################

Kim Dill-Mcfarland
University of Washington, kadm@uw.edu

Copyright (C) 2019 Kim Dill-Mcfarland

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Input parameters:
REQUIRED
  OTU.table = OTU count table with sample names as ROWNAMES and OTUs as additional
              columns.
  metadata = metadata with sample names as ROWNAMES and any factors of interest as 
              additional columns
  factors = factors from metadata table that you want to test, input as a vector/list
  sim.method = Beta-diversity measure. Options are Bray-Curtis 'bray', Jaccard 'jaccard',
                weighted UniFrac 'wunifrac' and unweighted UniFrac 'uwunifrac'
  p.adjust.m = multiple comparison correction of p-values across all pairwise comparisons.
                Methods include all those in `p.adjust()`
  
OPTIONAL
  tree = taxonomic tree of OTUs for UniFrac. Default is NULL
  taxonomy = taxonomic identities of OTUs for UniFrac. OTU names as ROWNAMES and columns as
              taxonomy levels (e.g. Phylum, Class, etc.) Default is NULL
  strata.factor = The variable from metadata that you want to stratify by. Default is NULL

Example
  pairwise.adonis.all(OTU.table=OTU, metadata=meta, 
                      factors=c('group1','group2'), 
                      sim.method='bray',
                      p.adjust.m='bonferroni', 
                      tree=NJ.tree, taxonomy=tax, 
                      strata.factor='Subject')
"
####

pairwise.adonis <- function(OTU.table, metadata, factors, sim.method, p.adjust.m,
                            tree=NULL, taxonomy=NULL, strata.factor=NULL)
{require(vegan)
 require(phyloseq)
require(dplyr)

####### Blank table for results #######
results <- data.frame(matrix(ncol = 5, nrow = 0))
colnames(results) <- c("factor","comparison","F.stat","Rsq","pval")

####### CHECK DATA #######
  #Data size  
  if (nrow(OTU.table) != nrow(metadata)){
    stop("ERROR: OTU table does not have the same number of samples as metadata.
       Please correct.")
  } 
  #Data type
  else if (!is.matrix(OTU.table) | !is.matrix(metadata)){
    stop("ERROR: OTU table and metadata must be matrices.")
  }
  else if (!is.numeric(OTU.table)){
    stop("ERROR: OTU table must be numeric,")
  }
  else {
    
    #Sort samples
    OTU.table <- OTU.table[order(rownames(OTU.table)), , drop=FALSE]
    metadata <- metadata[order(rownames(metadata)), , drop=FALSE]
    #Check if orders match
    if(!identical(rownames(OTU.table), rownames(metadata))){ 
      stop("ERROR: Sample names in your OTU table and metadata do not match.
         Please make sure sample names are in the first column of each.")
    }
    else {
      print("Data and samples are formatted correctly. Continuing with adonis.")
    }
  }
  
####### Calculate beta-diversity measures #######
  if(sim.method == "bray")
  {dist.mat = as.matrix(vegdist(OTU.table, distance="bray"))}
  
  else if(sim.method == "jaccard")
  {dist.mat = as.matrix(vegdist(OTU.table, distance="jaccard"))}
  
  else if(sim.method %in% c("wunifrac","uwunifrac"))
  {#Create phyloseq object
    phy.object <- phyloseq(otu_table(OTU.table, taxa_are_rows=FALSE), 
                           tax_table(taxonomy),
                           phy_tree(tree))
    
    #Then calculate distances
    if(sim.method == "wunifrac"){
      dist.mat = as.matrix(UniFrac(phy.object, 
                                   weighted=TRUE, normalized=TRUE))
    }else{
      dist.mat = as.matrix(UniFrac(phy.object, 
                                   weighted=FALSE, normalized=TRUE))}}else{print("ERROR: Please define valid beta-diversity measure.")}
  
####### Create pairwise comparisons #######
  for (factor1 in factors){
    #List all levels of factor interest
    levels <- unique(metadata[,factor1])
        ##Remove NAs if exist
        levels <- levels[!is.na(levels)]
    #Create all pairwise comparisons
    co <- as.matrix(combn(unique(levels), 2))
    
    #Subset data to samples for 1 pairwise comparison for 1 factor of interest
    for(elem in 1:ncol(co)){
      #Define two levels of factor
      level1 <- as.character(co[1, elem])
      level2 <- as.character(co[2, elem])
      
      #Determine which samples are in two-level subset
      samples.lv1.2 <- metadata[,factor1] %in% c(level1,level2)
      
####### Run adonis #######
      #For stratified models (like repeated measure)
      if(!is.null(strata.factor)){
        
        ad <- adonis(dist.mat[samples.lv1.2, samples.lv1.2]
                    ~ metadata[samples.lv1.2, factor1], 
                    strata=metadata[samples.lv1.2, strata.factor])
      }
      #Unstratified models
      else{
        ad <- adonis(dist.mat[samples.lv1.2, samples.lv1.2]
                    ~ metadata[samples.lv1.2, factor])
      }
      #Extract results 
      result.elem <- data.frame(
        factor=factor1,
        comparison=paste(co[1,elem], co[2,elem], sep="_"),
        F.stat=ad$aov.tab[1,4],
        Rsq=ad$aov.tab[1,5],
        Pval=ad$aov.tab[1,6])
      #Append to results table
      results <- rbind(results, result.elem)
    }#Complete elem loop
  }#complete factor loop

####### Calculate FDR corrected P-values within each factor #######
results.FDR <- results %>% 
  group_by(factor) %>% 
  mutate(FDR = p.adjust(pval, method=p.adjust.m))

####### Save results to global environment #######
assign("pairwise.adonis.result", results.FDR, env=.GlobalEnv)

}#End fxn

