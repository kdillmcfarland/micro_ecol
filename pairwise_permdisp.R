"Pairwise PERMDISP: Permutational tests of beta-dispersion (PERMDISP) across multiple variables of interest and pairwise with TukeyHSD

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
name: Beginning of output file name. Can include filepath
"

pairwise.permdisp <- function(x, y, factors, sim.method, name)
{require(vegan)
 require(tidyverse)
 # Vectors to hold results
 p.value = c()
 
  for(factor in 1:length(factors)){
      #Subset data to unlabeled OTU table and arrange
      x.subset = x %>% 
        arrange(Sample) %>% 
        select(starts_with('Otu'))
      y.subset = y %>% 
        arrange(Sample) %>% 
        select(factors[factor])
      
      #Calculate beta distances
      dist.mat = vegdist(x.subset, method=sim.method) 
      
      if(factors[factor] == "Random"){
      #Run PERMDSIP only on subset data
      perm = permutest(betadisper(dist.mat, get(factors[factor], y.subset)), permutations=1000)
      #Pull out p-values into list
      p.value = c(perm$tab$`Pr(>F)`[1])
      #Print results to a table
      results = data.frame(
        pval = p.value,
        label = factors[factor])
      
      write.table(results, file=paste(name,'permdisp.txt',sep="_"), append=TRUE, sep="\t", col.names = FALSE,  row.names=FALSE)  
    }
    
    else{
      #Run PERMDSIP and Tukey on subset data
      perm = permutest(betadisper(dist.mat, get(factors[factor], y.subset)), permutations=1000)
      tukey = TukeyHSD(betadisper(dist.mat, get(factors[factor], y.subset)))
      #Pull out R-squared and p-values into list
      p.value = c(perm$tab$`Pr(>F)`[1], tukey$group[,4])
      #Print results to a table
      results = data.frame(
              pval = p.value,
              label = c(factors[factor], labels(tukey$group[,4])))
    
    write.table(results, file=paste(name,'permdisp.txt',sep="_"), append=TRUE, sep="\t", col.names = FALSE,  row.names=FALSE)
    }}

  #After 1 variable in factors, read in results table. This will be added to the the append=TRUE option with the next factor's loop
  y=read.table(paste(name,'permdisp.txt',sep="_"), header=FALSE, sep="\t", fill = TRUE)
  #Remove the file
  file.remove(paste(name,'permdisp.txt',sep = "_"))
  #Relabel columns
  colnames(y) = c("pval", "label")
  #Write final results file
  write.table(y, file=paste(name,'permdisp.txt', sep="_"), row.names=FALSE, sep="\t")
}

#How to use
#pairwise.permdisp(OTU, meta, c("variable1", "variable2"), "bray", "filename")