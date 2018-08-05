"Pairwise ADONIS: Pairwise permutational ANOVA (PERMANOVA, adonis) within ONE variable of interest (ANY beta measure)

Kim Dill-Mcfarland
University of British Columbia

Copyright (C) 2018 Kim Dill-Mcfarland

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

Input files:

OTU.table = OTU count table with OTUs as column names and samples as rownames. No data other than OTU counts may be present within the dataframe
metadata = metadata with samples as rownames and any variables of interest as columns
factors = variable from metadata table that you want to test if beta-diversity varies by, input as meta$variable
sim.method = Beta-diversity measure. Options are Bray-Curtis 'bray', Jaccard 'jaccard', weighted UniFrac 'wunifrac' and unweighted UniFrac 'uwunifrac'
tree = taxonomic tree of OTUs to be used in creating phyloseq object for UniFrac calculation
stratify = Do you want to run a PERMANOVA model using the 'strata' option? TRUE or FALSE
strata.factor = If strafify = TRUE, the variable from the metadata that you want to stratify by, input as meta$variable
p.adjust.m = multiple comparison correction of p-values across all pairwise comparisons. Methods include all those in `p.adjust`

Example
  pairwise.adonis.all(OTU, meta, meta$group, 'bray', NJ.tree, TRUE, meta$subject, 'bonferroni')
"

pairwise.adonis.all <- function(OTU.table, metadata, factors, sim.method, tree, stratify, strata.factor, p.adjust.m)
{require(vegan)
 require(phyloseq)
  
  #Subset data to remove NAs
  factors.complete=factors[complete.cases(factors)]
  meta.complete = metadata[complete.cases(factors),]
  OTU.complete = OTU.table[row.names(meta.complete),]
  
  #Create phyloseq object
  OTU.UF = otu_table(OTU.complete, taxa_are_rows=FALSE)
  NJ.tree.UF = tree
  phy.object = phyloseq(OTU.UF, Tax.UF, NJ.tree.UF)
  
  #Create objects to hold output
  co = combn(unique(as.character(factors.complete)),2)
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()

 #Calculate beta-diversity measures
 if(sim.method == "bray"){dist.mat = as.matrix(vegdist(OTU.complete, distance="bray"))}
 if(sim.method == "jaccard"){dist.mat = as.matrix(vegdist(OTU.complete, distance="jaccard"))}
 if(sim.method == "wunifrac"){dist.mat = as.matrix(UniFrac(phy.object, weighted=TRUE, normalized=TRUE))}
 if(sim.method == "uwunifrac"){dist.mat = as.matrix(UniFrac(phy.object, weighted=FALSE, normalized=TRUE))}
  
  #Run all pairwise permanova comparisons, subsampling the data to include the two levels of interest
  for(elem in 1:ncol(co)){
    #For stratified models (like repeated measure)
    if(stratify == TRUE){
      ad = adonis(dist.mat[factors.complete %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                   factors.complete %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
          ~ factors.complete[factors.complete %in% c(as.character(co[1,elem]),as.character(co[2,elem]))], 
          strata=strata.factor[factors.complete %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])
    }
    #Unstratified models
    else{ad = adonis(dist.mat[factors.complete %in% c(as.character(co[1,elem]),as.character(co[2,elem])),
                     factors.complete %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]
          ~ factors.complete[factors.complete %in% c(as.character(co[1,elem]),as.character(co[2,elem]))]);
      pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
      F.Model = c(F.Model, ad$aov.tab[1,4]);
      R2 = c(R2,ad$aov.tab[1,5]);
      p.value = c(p.value, ad$aov.tab[1,6])
  }}
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)
}