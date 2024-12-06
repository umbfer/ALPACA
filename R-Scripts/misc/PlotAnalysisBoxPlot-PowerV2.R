library(rjson)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
 
 
setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset5-1000")

# produce 18 grafici (i.e.  PowerBox-PatternTransfer-A=0.10-G=0.10.png) per ogni AM, per ogni alpha e per ogni gamma
# ogni grafico contiene i boxplot della power statistics per ciascun misura e per ogni valore di k (obsoleto)

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

    seqLengths <- 50 # 200.000 - 10.000.0000 step 200.000 # numero di risultati per ciascun plot
    nAlpha <- 1
    nGamma <- length(gammaValues) # cambiare in 3
    nK <- length(kValues)
    nPlot <- nAlpha * nK  # nGamma * nAlpha * nK 
    nCol <- 6
    nCol2 <- 6
	colors <- rainbow(nK + 1)
	
	for (am in c("MotifReplace", "PatternTransfer")) {
		
		files <- list.files(jsonDir, sprintf("*%s*", am))
	    print( sprintf("Processing %d files from: %s/%s", length(files), getwd(), jsonDir))

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
	
			measureName <- gsub( "[dD]istance", "", exp$header$distanceName)
			tblMeasures[ndx] = measureName
	
			# aggiunge una colonna con il nome della misura in ogni riga
			# df["measure"] = ndx
			df["measure"] = measureName
						
			dfAll <- rbind(dfAll, df)
			
			ndx <- ndx + 1
			cat("done.\n")
		}
		
		# stop("break")
		
		dfAll$measure <- as.factor(dfAll$measure)
		dfAll$k <- factor( dfAll$k, levels = c( 4, 6, 8, 10))
		# dfAll$k <- as.factor(dfAll$k)
		# saveRDS( dfAll, file = sprintf("AllData-%s.RDS", am))
		for (gv in gammaValues) {	
		    title = sprintf("Power Statistic (Alpha =%d%%, G = %d%%, A.M. = %s)", dfAll$alpha[1]*100, gv*100, exp$header$alternateModel)
		
			sp2 <-	ggplot(filter( dfAll, gamma == gv), aes(x = measure, y = power, fill = k)) +
					ggtitle( title) + labs( x = "") +
					# theme(axis.text.x = element_blank()) + # axis.ticks.x = element_blank()) +
					theme_light() + theme(legend.position = "bottom") + theme(axis.text.x = element_text(angle = 45, vjust = 0.95, hjust=1)) +
					# scale_colour_brewer(palette="Dark2") +
					scale_fill_manual(values = c("#1C8F63", "#CE4907", "#6159A3", "#DD0077", "#535353")) +
					scale_y_continuous(name = "Power", limits = c(0, 1)) +
					geom_boxplot(lwd = 0.2, alpha = 0.9) # width=15, outlier.shape=20
		
			# dev.new(width = 10, height = 5)
		    outfname <- sprintf( "%s/PowerBox-%s-A=%.2f-G=%.2f.png", dirname, am, dfAll$alpha[1], gv)
		    ggsave( outfname, device = png(), width = 15, height = 10, units = "cm", dpi = 300)
		    # print( sp2)
		    # stop("break")
		    # readline(prompt="Press [enter] to continue")
		    dev.off() #only 129kb in size
		} # for each gamma
	} # for each alternate model
} # for each alpha