library(rjson)
library(ggplot2)
library(dplyr)

setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset5-1000")

# varKDataPath <- "/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset4-1000"

alphaValues = c("010", "050", "100")
# alphaValues = c("050")
gammaValues = c(0.01, 0.05, 0.10)

kValues = c(4, 6, 8, 10)

threshold = 0.7


for(alpha in alphaValues) {

    jsonDir <- sprintf( "json-%s", alpha)

    files <- list.files(jsonDir, "*.json")

    print( sprintf("Processing %d files from: %s/%s", length(files), getwd(), jsonDir))

    seqLengths <- 50 # da 200.000 a 10.000.000 step 200.000 = # numero di risultati per ciascun plot
    brks <- seq(1000000, 10000000, 1000000)
    ticks <- c()
    for( v in brks) {
    	ticks <- append(ticks, sprintf("%s", format(as.numeric(v), nsmall=0, big.mark=",", scientific=FALSE)))
    }
    nAlpha <- 1
    nGamma <- length(gammaValues) # cambiare in 3
    nK <- length(kValues)
    nPlot <- nAlpha * nK + 1 # nGamma * nAlpha * nK + 1
    nCol <- 6
	nCol2 <- 6
	colors <- rainbow(nK)
	
    for (file in files) {

        cat( sprintf("Processing: %s, Alpha: 0.%s ... ", file, alpha))

        exp <- fromJSON(file= paste(jsonDir, file, sep="/"))
        values <-exp[['values']]

        df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=nCol, byrow=T))})
        df <- do.call(rbind, df)
        colnames(df)[1:nCol] <- c("len", "alpha", "k", "power", "T1", "gamma")

        alphastr <- sprintf( "%.3f", df$alpha[1])
        dirname <- sprintf( "PointLenK-a=%s", substr(alphastr, 3, 5))
        if (!dir.exists(dirname)) {
            dir.create(dirname)
        }

		cat( " Gamma = ")
		for( gv in gammaValues) {
			cat( sprintf(" %.3f, ", gv))
				
	        title = sprintf("%s-%s-%s G=%.2f", gsub( "[dD]istance", "", exp$header$distanceName),
	        			exp$header$alternateModel, substr(exp$header$nullModel, start=1, stop=2), gv)

			data <- filter(df, power >= threshold, gamma == gv)
#			stop("break")
			if (nrow( data) == 0 ) {
				cat("no data")
				data <- data.frame(len = c(1000000), power = c(0.8), k = c(4))
				title <- sprintf("%s-NO Data!!!", title)
			}
			# un grafico per ogni valore di gamma con gli n valori di k
			sp2 <- ggplot( data, aes(x = len, y = k, color=factor(k))) +
				# ggtitle( title) + 
				scale_y_continuous(name = "Length(k)", limits = c(2,10), breaks = c(4, 6, 8, 10)) +
				scale_x_continuous("Sequence Len", limits = c(100000, 10000000), breaks = brks, labels = ticks) + # labels = function(x) format(x, scientific = FALSE)) +
				theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.95, hjust=1)) +

				geom_point(size = 1) +
				labs(title = title) # , color = "lenght(k)")

			lbls <- c()
			clrs <- c()
			levs <- levels(factor(data$k))
			for( kv in levs) {				
				lbls <- append(lbls, sprintf("k=%s", kv))
				clrs <- append(clrs,  colors[match(kv, kValues)]) # per usare gli stessi colori per ciascun k
			}
			# sp2 <- sp2 + scale_color_manual(breaks = lbls, values = colors, name = "k-len")
			sp2 <- sp2 + scale_color_manual(labels = lbls, breaks = levs, values = clrs, name = "Length") +
			guides(shape = FALSE) # remove the legend for shapes


			dev.new(width = 6, height = 4)
			print(sp2)
	        outfname <- sprintf( "%s/%s-G=%.2f.png", dirname, tools::file_path_sans_ext(file), gv)
	        # ggsave( outfname, device = png(), width = 15, height = 10, units = "cm", dpi = 300)
	        # readline(prompt="Press [enter] to continue")
	        Sys.sleep(5)
	        dev.off() #only 129kb in size
	    }
        cat( sprintf(" done.\n"))
    }
}