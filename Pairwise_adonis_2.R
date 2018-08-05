"Pairwise ADONIS: Pairwise permutational ANOVA (PERMANOVA, adonis) within MULTIPLE variables of interest (Bray-Curtis or Jaccard only)

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

x: OTU table. Counts with rows as samples and columns as specied. OTU names must start with 'Otu'. Sample names must be under 'Sample' variable
y: metadata table. Sample names must be under 'Sample' variable and must match the OTU table
factors: List of variables of interest in metadata (y) for testing in PERMDISP
sim.method: Similarity calculation method such as 'bray' or 'jaccard'
p.adjust.m: p-value adjustment method like 'fdr'
name: Beginning of output file name. Can include filepathsim.method: Similarity calculation method such as 'bray' or 'jaccard'
"

pairwise.adonis <- function(x, y, factors, sim.method, p.adjust.m, name)
{require(vegan)
 require(tidyverse)
  
  # Create matrix of all unique pairwise comparisons within a given variable in factors
  for(factor in 1:length(factors)){
    levels = y[[factors[factor]]]
    co = as.matrix(combn(unique(levels), 2))
    #Create lists to hold results
    pairs = c()
    R2 = c()
    p.value = c()
  
  for(elem in 1:ncol(co)){
    #Subset data to samples for 1 pairwise comparison for 1 variable in factors
      level1 = as.character(co[1, elem])
      level2 = as.character(co[2, elem])
        
      x.subset = x %>% 
        dplyr::filter(x[factors[factor]] == level1 |
                      x[factors[factor]] == level2) %>% 
        arrange(Sample) %>% 
        select(starts_with('Otu'))
      
      y.subset = y %>% 
        dplyr::filter(y[factors[factor]] == level1 |
                      y[factors[factor]] == level2) %>% 
        arrange(Sample) %>% 
        select(factors[factor])
      
      #Calculate beta distances
      dist.mat = vegdist(x.subset, method=sim.method) 
      
      #Run adonis on subset data
      ad = adonis(dist.mat ~ get(factors[factor], y.subset))
      
      #Create names for labeling results
      pairs = c(pairs, paste(co[1,elem], co[2,elem], sep="_"));
      #Pull out R-squared and p-values into list
      R2 = c(R2, ad$aov.tab[1,5]);
      p.value = c(p.value, ad$aov.tab[1,6])
      }
    
    #Print 1 variable in factors results to a table
    results = data.frame(
              comparison = pairs,
              Rsq = R2,
              pval = p.value,
              p.adj = p.adjust(p.value, method=p.adjust.m)
              )
    
    write.table(results, file=paste(name,'adonis.txt',sep="_"), append=TRUE, sep="\t", col.names = FALSE,  row.names=FALSE)
    }

  #After 1 variable in factors, read in results table. This will be added to the th4e append=TRUE option with the next factor's loop
  y=read.table(paste(name,'adonis.txt',sep="_"), header=FALSE, sep="\t", fill = TRUE)
  #Remove the file
  file.remove(paste(name,'adonis.txt',sep = "_"))
  #Relabel columns
  colnames(y) = c("Comparison", "Rsq", "pval", "p.adj")
  #Write final results file
  write.table(y, file=paste(name,'adonis.txt', sep="_"), row.names=FALSE, sep="\t")
}

#How to use
#pairwise.adonis(OTU, meta, c("variable1", "variable2"), "bray", "fdr", "filename")