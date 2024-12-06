

setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov")

cnt <- 0

for( k in c(4, 6, 8, 10)) {

	for( len in c(200000, 2000000)) {
		dfFilename = sprintf("Kolmogorov-k=%d-len=%d.df", k, len)
		ris <- readRDS(file = dfFilename)
		ris['K'] <- k
		ris['len'] <- len
	
		if (cnt == 0) {
			 final <- data.frame(ris)
		}	
		else {
			final <- rbind(ris, final)
		}	
		cnt <- cnt + 1
		cat( sprintf("%d - %d - %d\n", k, len, cnt))
	}
	
}

# sostituisce le labels per il modello
pat1 <- 'MotifRepl-U'
pat2 <- 'PatTransf-U'
pat3 <- 'Uniform-T1'
pat4 <- 'Uniform'
pat5 <- "-20{5,6}"
for(row in 1:nrow(final)) {
	final[row,1] <- gsub( pat1, "MR", final[row, 1])
	final[row,1] <- gsub( pat2, "PT", final[row, 1])
	final[row,1] <- gsub( pat3, "T1-NM", final[row, 1]) # attenzione all'ordine
	final[row,1] <- gsub( pat4, "NM", final[row, 1])
	final[row,1] <- gsub( pat5, "", final[row, 1])
 	final[row,5] <- if (final[row,5] == 200000) '200,000' else if (final[row, 5] == 2000000) '2,000,000' else final[row, 5]
}

saveRDS(final, file = 'AllKolmogorov.df')


	