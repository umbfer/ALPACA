library(rjson)
library(ggplot2)
library(dplyr)

# produce un grafico a linea con la power statistics (kulczynski2-PatternTransfer-Un-G=0.01.png) per ciascuna misura per ciascun valore di gamma
# sulle ascisse le lunghezze delle sequenze, sulle ordinate la power 4 linee una per ogni valore di k
# (obsoleto)

setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset5-1000")

alphaValues = c("010", "050", "100")
# alphaValues = c("050")
gammaValues = c(0.01, 0.05, 0.10)

kValues = c(4, 6, 8, 10)

for(alpha in alphaValues) {

    jsonDir <- sprintf( "json-%s", alpha)

    files <- list.files(jsonDir, "*.json")

    print( sprintf("Processing %d files from: %s/%s", length(files), getwd(), jsonDir))

    seqLengths <- 50 # da 200.000 a 10.000.000 step 200.000 = # numero di risultati per ciascun plot
    nAlpha <- 1
    nGamma <- length(gammaValues) # cambiare in 3
    nK <- length(kValues)
    nPlot <- nAlpha * nK + 1 # nGamma * nAlpha * nK + 1
    nCol <- 6
	nCol2 <- 6
	colors <- rainbow(nK + 1)

    for (file in files) {

        cat( sprintf("Processing: %s, Alpha: 0.%s ... ", file, alpha))

        exp <- fromJSON(file= paste(jsonDir, file, sep="/"))
        values <-exp[['values']]

        df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=nCol, byrow=T))})
        df <- do.call(rbind, df)
        colnames(df)[1:nCol] <- c("len", "alpha", "k", "power", "T1", "gamma")

        alphastr <- sprintf( "%.3f", df$alpha[1])
        dirname <- sprintf( "imgs-a=%s", substr(alphastr, 3, 5))
        if (!dir.exists(dirname)) {
            dir.create(dirname)
        }

		cat( " Gamma = ")
		for( gv in gammaValues) {
			cat( sprintf(" %.3f, ", gv))
				
	        title = sprintf("%s-%s-%s G=%.2f", gsub( "[dD]istance", "", exp$header$distanceName),
	        			exp$header$alternateModel, substr(exp$header$nullModel, start=1, stop=2), gv)

			k4 <- filter(df, k == 4, gamma == gv)

			# un grafico per ogni valore di gamma con gli n valori di k
			# sp2 <- ggplot( subset(df, gamma == gv)) +
			#	geom_line( aes( x = len, y = power, colour = k)) 
			sp2 <- ggplot( k4, aes(x=len)) + ggtitle( title) +  
				scale_y_continuous(name = "Power Statistics", limits = c(0, 1)) +
				scale_x_continuous("Sequence Len", limits = c(10000, 10000000))

			sp2 <- sp2 + geom_line( aes( y = power, colour = "k=4")) # + geom_point(aes(y = power))
			sp2 <- sp2 + geom_line( aes( y = filter(df, k == 6, gamma == gv)$power, colour = "k=6")) # + geom_point(aes(y = filter(df, k == 6, gamma == gv)$power, colour = "k=6"))
			sp2 <- sp2 + geom_line( aes( y = filter(df, k == 8, gamma == gv)$power, colour = "k=8"))
			sp2 <- sp2 + geom_line( aes( y = filter(df, k == 10, gamma == gv)$power, colour = "k=10"))

			lbls <- c()
			for( kv in kValues) {				
				lbl = sprintf("k=%d", kv)
				lbls <- append(lbls, lbl)
#				x[[lbls[i]]] <- filter(df, k == kv, gamma == gv)
#				dfk <- filter(df, k == kv, gamma == gv)
#				sp2 <- sp2 + geom_line( aes( y = filter(df, k==kv, gamma==gv)$power, colour = lbls[i]))
			}
			sp2 <- sp2 + scale_colour_manual(breaks = lbls, values = colors, name = "k-len")
	
			# dev.new(width = 6, height = 4)
	        outfname <- sprintf( "%s/%s-G=%.2f.png", dirname, tools::file_path_sans_ext(file), gv)
	        ggsave( outfname, device = png(), width = 15, height = 10, units = "cm", dpi = 300)
	        # readline(prompt="Press [enter] to continue")
	        dev.off() #only 129kb in size
	    }
        cat( sprintf(" done.\n"))
    }
}