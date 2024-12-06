library(ggplot2)
library(facetscales)

setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov")

dfFilename <- 'AllKolmogorov.df'
ris <- readRDS(file = dfFilename)

pat1 <- 'MotifRepl-U'
pat2 <- 'PatTransf-U'
pat3 <- 'Uniform'
pat4 <- 'Uniform-T1'
pat5 <- "-20{5,6}"
for(row in 1:nrow(ris)) {
	ris[row,1] <- gsub( pat1, "MR", ris[row, 1])
	ris[row,1] <- gsub( pat2, "PT", ris[row, 1])
	ris[row,1] <- gsub( pat3, "G-NM", ris[row, 1])
	ris[row,1] <- gsub( pat4, "C-NM", ris[row, 1])
	ris[row,1] <- gsub( pat5, "", ris[row, 1])
}

saveRDS(ris, file = dfFilename)

