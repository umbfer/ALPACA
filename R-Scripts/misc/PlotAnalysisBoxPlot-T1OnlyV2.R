library(rjson)
library(ggplot2)
library(dplyr)
library(RColorBrewer)

# produce 3 grafici (T1Box-A=0.10.png) con i valori del type 1 error uno per ciascun Alpha (obsoleto)

setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset5-1000")

alphaValues = c("010", "050", "100")
# alphaValues = c("050")
gammaValues = c(0.01, 0.05, 0.10)

kValues = c(4, 6, 8, 10)

for(alpha in alphaValues) {

    jsonDir <- sprintf( "json-%s", alpha)

    files <- list.files(jsonDir, "*.json")

    print( sprintf("Processing %d files from: %s/%s", length(files), getwd(), jsonDir))

    dirname <- "T1BoxPlot"
    if (!dir.exists(dirname)) {
    	dir.create(dirname)
    }

    seqLengths <- 50 # 1000, 10000, 100000, 1.000.000 = # numero di risultati per ciascun plot
    nAlpha <- 1
    nGamma <- length(gammaValues) # cambiare in 3
    nK <- length(kValues)
    nPlot <- nAlpha * nK + 1 # nGamma * nAlpha * nK + 1
    nCol <- 6
    nCol2 <- 6
	colors <- rainbow(nK + 1)

	dfAll <- data.frame()
	tblMeasures <- c()
	ndx = 1
	
    for (file in files) {

        if (grepl( "jaccard|mash", file, ignore.case = TRUE)) {
            cat(sprintf("skipping file: %s\n", file))
            next
        }
        cat( sprintf("Processing: %s, Alpha: 0.%s ... ", file, alpha))

        exp <- fromJSON(file = paste(jsonDir, file, sep="/"))
        values <-exp[['values']]
		
        df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=nCol, byrow=T))})
        df <- do.call(rbind, df)
        colnames(df)[1:nCol] <- c("len", "alpha", "k", "power", "T1", "gamma")

#		if (exp$header$alternateModel == "PatternTransfer") {
#			print("skipped.")
#			next
#		}

		measureName <- gsub( "[dD]istance", "", exp$header$distanceName)
		tblMeasures[ndx] = measureName

		# aggiunge una colonna con il nome della misura in ogni riga
		# df["measure"] = ndx
		df["measure"] = measureName
				
		# aggiunge al dataframe finale i valori per un solo gamma N.B. T1 check
		gv <- 0.010
		dfAll <- rbind(dfAll, filter( df, gamma == gv))
		
		ndx <- ndx + 1
		cat("done.\n")
	}
	
	dfAll$measure <- as.factor(dfAll$measure)
	dfAll$k <- factor( dfAll$k, levels = c( "var", 4, 6, 8, 10))
	# dfAll$k <- as.factor(dfAll$k)
		
    title = sprintf("T1 Error Check (alpha = %d%%)", dfAll$alpha[1]*100)

	sp2 <-	ggplot(dfAll, aes(x = measure, y = T1, fill = k)) +
			ggtitle( title) + labs( x = "") +
			theme_light() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) +
			# theme(axis.text.x = element_blank()) + # axis.ticks.x = element_blank()) +
			# scale_colour_brewer(palette="Dark2") +
			scale_fill_manual(values = c("#1C8F63", "#CE4907", "#6159A3", "#DD0077", "#535353")) +
			scale_y_continuous(name = "T1 Error", limits = c(0, 0.30)) +
			geom_boxplot(lwd = 0.2, alpha = 0.9) #outlier.shape = 20, width=15

	# dev.new(width = 10, height = 5)
    outfname <- sprintf( "%s/T1Box-A=%.2f.png", dirname, dfAll$alpha[1])
    ggsave( outfname, device = png(), width = 9, height = 6, units = "cm", dpi = 300)
    # print( sp2)
    # readline(prompt="Press [enter] to continue")
    dev.off() #only 129kb in size
} # for each alpha