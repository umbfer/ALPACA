library(rjson)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
 
 
setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/syntheticAllK")

varKDataPath <- "/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset4-1000"


alphaValues = c("010", "050", "100")
# alphaValues = c("050")
gammaValues = c(0.01, 0.05, 0.10)

kValues = c(4, 6, 8, 10)

for(alpha in alphaValues) {

    jsonDir <- sprintf( "json-%s", alpha)

    dirname <- "PowerBoxPlot"
    if (!dir.exists(dirname)) {
    	dir.create(dirname)
    }

    seqLengths <- 17 # 1000, 10000, 100000, 1.000.000 = # numero di risultati per ciascun plot
    nAlpha <- 1
    nGamma <- length(gammaValues) # cambiare in 3
    nK <- length(kValues)
    nPlot <- nAlpha * nK + 1 # nGamma * nAlpha * nK + 1
    nCol <- 6
    nCol2 <- 6
	colors <- rainbow(nK + 1)
	measureList <- c()
			
		files <- list.files(jsonDir, "^d2*")
	    print( sprintf("Processing %d files from: %s/%s", length(files), getwd(), jsonDir))

		dfAll <- data.frame()
		tblMeasures <- c()
		ndx = 1
	    
	    for (file in files) {
	
	        cat( sprintf("Processing: %s, Alpha: 0.%s ... ", file, alpha))
	
	        exp <- fromJSON(file = paste(jsonDir, file, sep="/"))
	        values <-exp[['values']]
			
	        df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=nCol, byrow=T))})
	        df <- do.call(rbind, df)
	        colnames(df)[1:nCol] <- c("len", "alpha", "k", "power", "T1", "gamma")
	
			measureName <- gsub( "[dD]istance", "", exp$header$distanceName)
			tblMeasures[ndx] = measureName
			if ( ! measureName %in% measureList) measureList <- c(measureList, measureName)

	
			# aggiunge una colonna con il nome della misura in ogni riga
			# df["measure"] = ndx
			df["measure"] = measureName
			df["model"] = exp$header$alternateModel

	
			varKFile <- sprintf("%s/json-%s/%s", varKDataPath, alpha, file)
	        expVarK <- fromJSON(file = varKFile)
	        valuesVarK <-expVarK[['values']]
	        dfVarK <- lapply( valuesVarK, function( p) { data.frame(matrix(unlist(p), ncol=nCol2, byrow=T))})
	        dfVarK <- do.call(rbind, dfVarK)
	        colnames(dfVarK)[1:nCol2] <- c("len", "alpha", "k", "power", "T1", "gamma")
	
			# aggiunge una colonna con il nome della misura in ogni riga
			# dfVarK[ "measure"] = ndx
			dfVarK["measure"] = measureName
			dfVarK["model"] = exp$header$alternateModel

			# sostituisce ai valori di k variabile una label "var"
			dfVarK["k"] = "var"
			
			dfAll <- rbind(dfAll, df)
			dfAll <- rbind(dfAll, dfVarK)
			
			ndx <- ndx + 1
			print("done")
		}
		
		dfAll$measure <- as.factor(dfAll$measure)
		# dfAll$k <- as.factor(dfAll$k)
		dfAll$k <- factor( dfAll$k, levels = c( "var", 4, 6, 8, 10))
		dfAll$model <- as.factor(dfAll$model)
		
		for( measureName in measureList) {
		# stop("break")

			for (gv in gammaValues) {	
			    title = sprintf("%s (A = %d%%, G = %d%%)", measureName, dfAll$alpha[1]*100, gv*100)
			
				sp2 <-	ggplot(filter(dfAll, measure == measureName, gamma == gv), aes(x = k, y = power, fill = model)) +
						ggtitle( title) + theme_light() + labs( x = "K len") + 
						theme(legend.position = "bottom") + # theme(axis.text.x = element_blank()) + # axis.ticks.x = element_blank()) +
						# theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) +
						# scale_colour_brewer(palette="Dark2") +
						scale_fill_manual(values = c("#1C8F63", "#CE4907", "#6159A3", "#DD0077", "#535353")) +
						scale_y_continuous(name = "Power", limits = c(0, 1)) +
						geom_boxplot(lwd = 0.5, alpha = 0.9) # outlier.shape = 20, 
			
				# dev.new(width = 10, height = 5)
			    outfname <- sprintf( "%s/Power%s-A=%.2f-G=%.2f.png", dirname, measureName, dfAll$alpha[1], gv)
			    ggsave( outfname, device = png(), width = 15, height = 10, units = "cm", dpi = 300)
			    # print( sp2)
			    # stop("break")
			    # readline(prompt="Press [enter] to continue")
			    dev.off() #only 129kb in size
			} # for each gamma
		} # for each measure
} # for each alpha