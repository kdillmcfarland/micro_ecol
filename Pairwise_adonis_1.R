"Pairwise ADONIS: Pairwise permutational ANOVA (PERMANOVA, adonis) within ONE variable of interest (Bray-Curtis or Jaccard only)

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

x: OTU table. Counts with rows as samples and columns as specied. OTU names must start with 'Otu'. No other data may be present within
factors: Variable of interest in metadata. Must be in the form metadata$variable
sim.method: Similarity calculation method such as 'bray' or 'jaccard'
p.adjust.m: p-value adjustment method like 'fdr'
"

pairwise.adonis <- function(x, factors, sim.method, p.adjust.m)
{require(vegan)
  co = as.matrix(combn(unique(factors), 2))
  pairs = c()
  F.Model =c()
  R2 = c()
  p.value = c()
  
  for(elem in 1:ncol(co)){
    ad = adonis(x[factors %in% c(as.character(co[1,elem]),
                                 as.character(co[2,elem])),] ~ 
                  factors[factors %in% c(as.character(co[1,elem]),
                                         as.character(co[2,elem]))],
                method =sim.method);
    pairs = c(pairs,paste(co[1,elem],'vs',co[2,elem]));
    F.Model =c(F.Model,ad$aov.tab[1,4]);
    R2 = c(R2,ad$aov.tab[1,5]);
    p.value = c(p.value,ad$aov.tab[1,6])}
  
  p.adjusted = p.adjust(p.value,method=p.adjust.m)
  pairw.res = data.frame(pairs,F.Model,R2,p.value,p.adjusted)
  return(pairw.res)}

#How to use
#pairwise.adonis(OTU, meta$variable, "bray", "fdr")
#pairwise.adonis(x, factors, "jaccard", "fdr")