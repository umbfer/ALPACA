library(DescTools)
library(ggplot2)
library(dplyr)


setwd("/home/cattaneo/results/Escherichiacoli")
# setwd('/Users/pipp8/Universita/Progetti/BioInformatica/Present-Absent/Src')

dfFilename = "HammingDistanceEC.df"

gammas <- c(0.005, 0.010, 0.050, 0.100, 0.200, 0.300, 0.500)

if (file.exists(dfFilename)) {
	cat( sprintf("Data file exists. Loading %s\n", dfFilename))			
	# carica il dataframe dal file
	df <- readRDS(file = dfFilename)
} else {

	df <- data.frame( Name = character(), Gamma = double(), Dist = double(), len = numeric(), stringsAsFactors=FALSE)

	model = 'EscherichiaColi'

	for( g in gammas)	{
	  f1 <- sprintf("%s.fasta", model)
	  con1 = file(f1, "r")
	  f2 <- sprintf("%s-G=%.3f.fasta", model, g)
		con2 = file(f2, "r")
    distance <- 0
    totLen <- 0

    cat( sprintf("Starting for gamma = %.3f\n", g))
		while( TRUE) {

			s1Name = readLines(con1, n = 1)
			s1 = readLines(con1, n = 1)
			s2Name = readLines(con2, n = 1)
			s2 = readLines(con2, n = 1)

			if (nchar(s1Name) == 0 | nchar(s2Name) == 0 | 
			    nchar(s1) == 0 | nchar(s2) == 0) {
			    break
			}
			if (startsWith(s1, ">") | startsWith(s2, ">")) {
			  stop("Malformed input file")
			}
			
			# compute hamming distance
			d <- StrDist(s1, s2, method = 'hamming')
			distance <- distance + d[1]
			totLen <- totLen + nchar(s1)
		}
		close(con1)
		close(con2)
		df[nrow(df) + 1,] <- list( model, g, distance/totLen, totLen)
		cat( sprintf("%s. G=%f d = %d. done.\n", model, g, distance/totLen))
	}
	saveRDS( df, file = dfFilename)
}


df$Gamma = as.factor(df$Gamma)

sp <- ggplot( dfNM, aes(x = Gamma, y = Dist)) +
  geom_line() +
  geom_point(aes(shape = 2, size = 3)) +
  labs(y = "Hamming Distance")


# dev.new(width = 6, height = 6)
outfname <- sprintf("JaccardEscherichiacoli-AllK.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 6, dpi = 300)
# print(sp)
# dev.off() #only 129kb in size

