library(rjson)
library(ggplot2)
library(dplyr)

setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/syntheticAllK")

alphaValues = c("010", "050", "100")
# alphaValues = c("050")
gammaValues = c(0.01, 0.05, 0.10)

kValues = c(4, 6, 8, 10)

for(alpha in alphaValues) {

    jsonDir <- sprintf( "json-%s", alpha)

    files <- list.files(jsonDir, "*.json")

    print( sprintf("Processing %d files from: %s/%s", length(files), getwd(), jsonDir))

    seqLengths <- 17 # 1000, 10000, 100000, 1.000.000 = # numero di risultati per ciascun plot
    nAlpha <- 1
    nGamma <- length(gammaValues) # cambiare in 3
    nK <- length(kValues)
    nPlot <- nAlpha * nK + 1 # nGamma * nAlpha * nK + 1
    nCol <- 6
	colors <- rainbow(nK + 1)

    for (file in files) {

        cat( sprintf("Processing: %s, Alpha: 0.%s ... ", file, alpha))

        exp <- fromJSON(file= paste(jsonDir, file, sep="/"))
        values <-exp[['values']]

        df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=nCol, byrow=T))})
        df <- do.call(rbind, df)
        colnames(df)[1:nCol] <- c("len", "alpha", "k", "power", "T1", "gamma")

        alphastr <- sprintf( "%.3f", df$alpha[1])
        dirname <- sprintf( "img-a=%s", substr(alphastr, 3, 5))
        if (!dir.exists(dirname)) {
            dir.create(dirname)
        }
		
		t <- switch(alpha, "010" = 20, "050" = 10, "100" = 5) # fattore di amplificazione del valore di T1
		print( sprintf("t = %d\n", t))

		gndx <- 0
		for( gv in gammaValues) {
				
	        title = sprintf("%s-%s-%s G=%.2f", gsub( "[dD]istance", "", exp$header$distanceName),
	        			exp$header$alternateModel, substr(exp$header$nullModel, start=1, stop=2), gv)

			gndx <- gndx + 1
			T1kv <- if (exp$header$alternateModel == "MotifReplace") kValues[gndx] else if (gndx < 3) kValues[gndx] else kValues[4] 
			T1lbl <- sprintf("T1 k=%d", T1kv)
			x <- c()
			k4 <- filter(df, k == 4, gamma == gv)
			shp <- c(1)
			
			# un grafico per ogni valore di gamma con gli n valori di k
			# sp2 <- ggplot( subset(df, gamma == gv)) +
			#	geom_line( aes( x = len, y = power, colour = k)) 
			sp2 <- ggplot( k4, aes(x=len)) + ggtitle( title) +  
				scale_y_continuous(name = "Power Statistics", limits = c(0, 1), 
					sec.axis = sec_axis(~. * t, name = "T1 Error Check", # 40 => 0.025 fondo scala
		 			labels = function(b) { paste0( round( b * 100 / (t*t), 1), "%") })) + 
				scale_x_continuous("Sequence Len", limits = c(1000, 10000000), trans='log10') +
#				scale_shape_manual(values=1:nK+1) +

				geom_hline(yintercept = df$alpha[1] * t, linetype="dashed", color = "gray25") +
				annotate(geom="text", x=3000000, y= df$alpha[1] * t + 0.500 / t,
					label=paste0("Alpha = ", df$alpha[1]), color="gray25") + 
				geom_line( aes( x = len, y = filter(df, k == T1kv, gamma == gv)$T1*t, colour = T1lbl))
				# geom_point(aes(x = len, y = T1*t, colour = T1lbl))

			sp2 <- sp2 + geom_line( aes( y = power, colour = "k=4")) # + geom_point(aes(y = power))
			sp2 <- sp2 + geom_line( aes( y = filter(df, k == 6, gamma == gv)$power, colour = "k=6")) # + geom_point(aes(y = filter(df, k == 6, gamma == gv)$power, colour = "k=6"))
			sp2 <- sp2 + geom_line( aes( y = filter(df, k == 8, gamma == gv)$power, colour = "k=8"))
			sp2 <- sp2 + geom_line( aes( y = filter(df, k == 10, gamma == gv)$power, colour = "k=10"))

			lbls <- c(T1lbl)
			for( kv in kValues) {				
				lbl = sprintf("k=%d", kv)
				lbls <- append(lbls, lbl)
#				x[[lbls[i]]] <- filter(df, k == kv, gamma == gv)
#				dfk <- filter(df, k == kv, gamma == gv)
#				sp2 <- sp2 + geom_line( aes( y = filter(df, k==kv, gamma==gv)$power, colour = lbls[i]))
			}
			sp2 <- sp2 + scale_colour_manual(breaks = lbls, values = colors, name = "k-len")
	
			# dev.new(width = 6, height = 4)
			cat( sprintf(" Gamma = %.3f, ", gv))
	        outfname <- sprintf( "%s/%s-G=%.2f.png", dirname, tools::file_path_sans_ext(file), gv)
	        ggsave( outfname, device = png(), width = 15, height = 10, units = "cm", dpi = 300)
	        # readline(prompt="Press [enter] to continue")
	        dev.off() #only 129kb in size
	    }
        print( " done")
    }
}