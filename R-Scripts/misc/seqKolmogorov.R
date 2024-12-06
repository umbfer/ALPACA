library(ggplot2)

# setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov")
setwd("/Volumes/Catalina/PowerStatistics/Kolmogorov")

ris <- data.frame( Name = character(), D = double(), PV = double(), stringsAsFactors=FALSE)

EvalKS <- function( model, k, len, gValues)
{
	for( g in gValues) {
		gVal <- if (g == "no g") "" else sprintf(".G=%.3f", g)

		writeLines(  sprintf("model: %s, k: %d, len: %d, %s", model, k, len, gVal))
		k1 <- sprintf("%s-%d%s", model, len, gVal)

		start <- nrow(ris) + 1
		for(pair in 1:nPairs) {
			f1 <- sprintf("k=%d/%s/%s/%d/dist-k=%d_%s-%d.%d%s-%05d-%s.dist",
						  k, seqDistDir, model, len, k, model, nPairs, len, gVal, pair, seqA)
			f2 <- sprintf("k=%d/%s/%s/%d/dist-k=%d_%s-%d.%d%s-%05d-%s.dist",
						  k, seqDistDir, model, len, k, model, nPairs, len, gVal, pair, seqB)

			df1 <- read.table(header=FALSE, f1)
			xv <- df1[,2]
			df2 <- read.table(header=FALSE, f2)
			yv <- df2[,2]
			ksRis <- ks.test(xv, yv)

			ris[nrow(ris) + 1,] <<- list( k1, as.numeric(ksRis[1]), as.numeric(ksRis[2]))
			# print( ksRis)
			cat( sprintf("%5d /%5d\r", pair, nPairs))
		}
		end <- nrow(ris)
		writeLines(' ')
		inFName = sub(pattern = ".*_(.*)-[0-9]*-[AB]\\..*$", replacement = "\\1", f1)
		processed = sub(pattern = ".*-([0-9]*)-[AB]\\..*$", replacement = "\\1", f1)
		writeLines( sprintf("%s, processed: %s pairs", inFName, processed))
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
for( k in c(4, 6, 8)) {

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
		
		title = sprintf("D value Test di Kolmogorov/Smirnov (len: %d, k=%d)", len, k)

		sp2 <-	ggplot(ris, aes(x = Name, y = D, color = Name)) +
			ggtitle( title) + labs( x = "") +
			# theme(axis.text.x = element_blank()) + # axis.ticks.x = element_blank()) +
			theme_light(base_size = 10) + theme(legend.position = "bottom") +
			theme(axis.text.x = element_text(size = 6, angle = 45, vjust = 0.95, hjust=1)) +
			scale_colour_brewer(palette="Dark2") +
			scale_y_continuous(name = "D", limits = c(0, 0.05)) +
			geom_boxplot(lwd = 0.2, alpha = 0.9) # width=15, outlier.shape=20

		# dev.new(width = 6, height = 9)
		outfname <- sprintf( "Kolmogorov-k=%d-len=%d.png", k, len)
		# ggsave( outfname, device = png(), width = 9, height = 6, units = "cm", dpi = 300)
		ggsave( outfname)
		# print( sp2)
		# stop("break")
		# readline(prompt="Press [enter] to continue")
		# dev.off() #only 129kb in size
	} # for each len
} # for each k
