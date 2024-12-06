library(ggplot2)

setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov")
# setwd("/Volumes/Catalina/PowerStatistics/Kolmogorov")

ris <- data.frame( Name = character(), D = double(), PV = double(), stringsAsFactors=FALSE)

EvalKS <- function( model, k, len, gValues)
{
	for( g in gValues) {
		gVal <- if (g == "no g") "" else sprintf(".G=%.3f", g)

		writeLines(  sprintf("model: %s, k: %d, len: %d, %s", model, k, len, gVal))
		k1 <- sprintf("%s-%d%s", model, len, gVal)

		start <- nrow(ris) + 1
		f1 <- sprintf("k=%d/%s/dist-k=%d_%s-%d.%d%s-All.distbin",
					  k, seqDistDir, k, model, nPairs, len, gVal)
		
		inFile <- file( f1, "rb")
		nKeys = readBin( inFile, integer(), n = 1)
		
		# df1 <- read.table(header=FALSE, f1)
		
		for(pair in 1:nPairs) {
			# f2 <- sprintf("k=%d/%s/%s/%d/dist-k=%d_%s-%d.%d%s-%05d-%s.dist",
			#			  k, seqDistDir, model, len, k, model, nPairs, len, gVal, pair, seqB)	
			# xv <- df1[, pair*2]
			# yv <- df1[, pair*2 + 1]
			xv <- readBin( inFile, double(), n = nKeys)
			yv <- readBin( inFile, double(), n = nKeys)
			ksRis <- ks.test(xv, yv)

			ris[nrow(ris) + 1,] <<- list( k1, as.numeric(ksRis[1]), as.numeric(ksRis[2]))
			# print( ksRis)
			cat( sprintf("%5d /%5d\r", pair, nPairs))
		}
		end <- nrow(ris)
		writeLines(' ')
		writeLines( sprintf("processed: %d pairs from: %s", nPairs, f1))
		writeLines( sprintf("Min(D):%f, Max(D):%f, Avg(D):%f, Var(D):%f",
							min(ris$D[c(start:end)]), max(ris$D[c(start:end)]), mean(ris$D[c(start:end)]), var(ris$D[c(start:end)])))
		writeLines(sprintf("Min(P):%f, Max(P):%f, Avg(P):%f, Var(P):%f",
						   min(ris$PV[c(start:end)]), max(ris$PV[c(start:end)]), mean(ris$PV[c(start:end)]), var(ris$PV[c(start:end)])))
	}
}


seqDistDir <- "seqDists"
seqA <- "A"
seqB <- "B"
nPairs <- 1000
# k <- 4
# len <- 100000

#for( k in c(4, 6, 8, 10)) {
for( k in c(4, 6, 8, 10)) {

	for( len in c(200000, 2000000)) {

		dfFilename = sprintf("Kolmogorov-k=%d-len=%d.df", k, len)		
		if (file.exists(dfFilename)) {
			cat( sprintf("Data file exists. Skipping k = %d, len = %d", k, len))	
			
			# carica il dataframe dal file
			ris <- readRDS(file = dfFilename)
		}
		else {
			# resetta il dataframe con i risultati
			ris <- data.frame( Name = character(), D = double(), PV = double(), stringsAsFactors=FALSE)
	
			model <- "MotifRepl-U"
			gValues <- c(0.010, 0.050, 0.100)
			EvalKS( model, k, len, gValues)
	
			model <- "PatTransf-U"
			gValues <- c(0.010, 0.050, 0.100)
			EvalKS( model, k, len, gValues)
	
			model <- "Uniform"
			gValues <- c("no g")
			EvalKS( model, k, len, gValues)
	
			model <- "Uniform-T1"
			gValues <- c("no g")
			EvalKS( model, k, len, gValues)
	
	#		stop("break")	
			colnames(ris) <- c("Name", "D", "PV")
			
			saveRDS(ris, file = dfFilename)
		}
		
		# pat1 <- sprintf("%d\\.+|-%d$", len, len)
		pat1 <- 'MotifRepl-U'
		pat2 <- 'PatTransf-U'
		pat3 <- 'Uniform'
		pat4 <- 'Uniform-T1'
		pat5 <- "-20{5,6}"
		for(row in 1:nrow(ris)) {
			ris[row,1] <- gsub( pat1, "MR", ris[row, 1])
			ris[row,1] <- gsub( pat2, "PT", ris[row, 1])
			ris[row,1] <- gsub( pat3, "GeneratingNM", ris[row, 1])
			ris[row,1] <- gsub( pat4, "ControlNM", ris[row, 1])
			ris[row,1] <- gsub( pat5, "", ris[row, 1])
		}

		# title = sprintf("Kolmogorov-Smirnov test for len=%d, k=%d", len, k)
		title = sprintf("k = %d", k)
		mv = max(ris[,2])
		if (mv >= 0.1)
			sf = 10
		else if (mv >= 0.01)
			sf = 100
		else
			sf = 1000

		fs = ceiling( mv * sf) / sf
		cat(sprintf("fs: %f <- %f", fs, mv))
		sp2 <-	ggplot(ris, aes(x = Name, y = D, color = Name)) +
			ggtitle( title) + labs( x = "") +
			labs(color = "Model") +
			# theme(axis.text.x = element_blank()) + # axis.ticks.x = element_blank()) +
			theme_light(base_size = 10) + 
			# theme(legend.position = "bottom") +
			theme(legend.position = "none") + 
			theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.95, hjust=1)) +
			scale_colour_brewer(palette="Dark2") +
			scale_y_continuous(name = "KS", limits = c(0, fs)) +
			geom_boxplot(lwd = 0.2, alpha = 0.9) # width=15, outlier.shape=20

		# dev.new(width = 6, height = 9)
		outfname <- sprintf( "Kolmogorov-k=%d-len=%d.png", k, len)
		ggsave( outfname, device = png(), width = 9, height = 6, dpi = 300)
		# ggsave( outfname)
		# print( sp2)
		# stop("break")
		# readline(prompt="Press [enter] to continue")
		# dev.off() #only 129kb in size
		cat( sprintf(" done.\n"))
	} # for each len
} # for each k
