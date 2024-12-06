library(rjson)
library(ggplot2)
library(dplyr)

# produce un grafico a linea con la T1 error check (harmonicmean-MotifReplace-UnT1.png) per ciascuna misura
# sulle ascisse le lunghezze delle sequenze, sulle ordinate il risultato T1 error check linee una per ogni valore di k
# (obsoleto)


setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset5-1000")


alphaValues = c("010", "050", "100")
# alphaValues = c("050")
gammaValues = c(0.01, 0.05, 0.10)

kValues = c(4, 6, 8, 10)

brks <- seq(1000000, 10000000, 1000000)
ticks <- c()
for( v in brks) {
	ticks <- append(ticks, sprintf("%s", format(as.numeric(v), nsmall=0, big.mark=",", scientific=FALSE)))
}

for(alpha in alphaValues) {

    jsonDir <- sprintf( "json-%s", alpha)

    files <- list.files(jsonDir, "*.json")

    print( sprintf("Processing %d files from: %s/%s", length(files), getwd(), jsonDir))

    dirname <- sprintf( "T1img-a=%s", alpha)
    if (!dir.exists(dirname)) {
    	dir.create(dirname)
    }

    seqLengths <- 50 # 1000, 10000, 100000, 1.000.000 = # numero di risultati per ciascun plot
    nAlpha <- 1
    nGamma <- length(gammaValues) # cambiare in 3
    nK <- length(kValues)
    nCol <- 6
	colors <- rainbow(nK)

    for (file in files) {

        cat( sprintf("Processing: %s, Alpha: 0.%s ... ", file, alpha))

        exp <- fromJSON(file = paste(jsonDir, file, sep="/"))
        values <-exp[['values']]
		
        df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=nCol, byrow=T))})
        df <- do.call(rbind, df)
        colnames(df)[1:nCol] <- c("len", "alpha", "k", "power", "T1", "gamma")

		if (exp$header$alternateModel == "PatternTransfer")
				next
		
		MaxT1 <- switch(alpha, "010" = 0.030, "050" = 0.150, "100" = 0.3) # fattore di amplificazione del valore di T1

        title = sprintf("T1 Error Check for: %s-%s", gsub( "[dD]istance", "", exp$header$distanceName), exp$header$nullModel)

		x <- c()
		gv <- 0.010
		
		k4 <- filter(df, k == 4, gamma == gv)
		shp <- c(1)
		
		# un solo grafico con i valori di T1-check per ogni valore di k
		#	geom_line( aes( x = len, y = power, colour = k)) 
		sp2 <- ggplot( k4, aes(x=len)) + ggtitle( title) + theme_bw() + # theme_classic() +
			scale_y_continuous(name = "T1 Error Check", limits = c(0, MaxT1)) +
			scale_x_continuous("Sequence Len", limits = c(100000, 10000000), breaks = brks, labels = ticks) + # labels = function(x) format(x, scientific = FALSE)) +
			theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.95, hjust=1)) +
			geom_hline(yintercept = df$alpha[1], linetype="dashed", color = "gray25") +
			annotate(geom="text", x = 9000000, y = 0.002, label=paste0("Alpha = ", df$alpha[1]), color="gray25") # df$alpha[1] + 0.002,

		
		sp2 <- sp2 + geom_line( aes( y = T1, colour = "k = 4")) # + geom_point(aes(y = power))
		sp2 <- sp2 + geom_line( aes( y = filter(df, k == 6, gamma == gv)$T1, colour = "k = 6"))
		# + geom_point(aes(y = filter(df, k == 6, gamma == gv)$power, colour = "k=6"))
		sp2 <- sp2 + geom_line( aes( y = filter(df, k == 8, gamma == gv)$T1, colour = "k = 8"))
		sp2 <- sp2 + geom_line( aes( y = filter(df, k == 10, gamma == gv)$T1, colour = "k = 10"))

		lbls <- c()
		for( kv in kValues) {				
			lbl = sprintf("k = %d", kv)
			lbls <- append(lbls, lbl)
		}
		sp2 <- sp2 + scale_colour_manual(breaks = lbls, values = colors, name = "")

		# dev.new(width = 6, height = 4)
		# print(sp2)
        outfname <- sprintf( "%s/%sT1.png", dirname, tools::file_path_sans_ext(file), gv)
        ggsave( outfname, device = png(), width = 15, height = 10, units = "cm", dpi = 300)
        # readline(prompt="Press [enter] to continue")
        dev.off() #only 129kb in size
	    
        print( " done")
    } # for each file
} # for each alpha