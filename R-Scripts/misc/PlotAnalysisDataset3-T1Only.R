library(rjson)
setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset2-100")

files <- list.files("json", "*.json")

print( paste("Processing", length(files), "files from:", getwd()))

seqLengths <- 12 # 1000, 10000, 100000, 1.000.000 = # numero di risultati per ciascun plot
nAlpha <- 1
nGamma <- 0 # cambiare in 3
nPlot <- (nGamma + 1) * nAlpha
nCol <- 6

for (file in files) {

	cat( "Processing: ", file)
	
    exp <- fromJSON(file= paste("json", file, sep="/"))
    values <-exp[['values']]

    df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=nCol, byrow=T))})
    df <- do.call(rbind, df)
	colnames(df)[1:nCol] <- c("len", "alpha", "k", "power", "T1", "gamma")

	outfname <- paste( "img/", tools::file_path_sans_ext(file), ".png", sep = "")
    png( outfname, width=15, height=10, units="cm", res=300)
    # par(mar=c(4,4,1,1))

	# title = paste(gsub( "[dD]istance", "", exp$header$distanceName), exp$header$alternateModel, exp$header$nullModel, sep = "-")
	title = paste(gsub( "[dD]istance", "", exp$header$distanceName), "T1 Check", sep = "-")
	xrange <- range(900, 15000000)	
	yrange <- range(0, 1)
    with( df, plot( xrange, yrange, type='n', xlab = "Sequence Length",
    				ylab = "Power Value", log="x", main = title )) # sub = "k = 6"
  	grid( 8, 5,col = "lightgray", lty = "dotted")
  	
  	abline(h = 0.05, lty = 4, col = "gray25")
	text(2500000,0.06, "alpha = 5%", col = "gray25", adj = c(0, -.1))

	leg <- seq(1,nPlot)
    cl <- seq(1,nPlot)
    plotchar <- seq(1,nPlot)
	linetype <- c(1:nAlpha)
    colors <- rainbow(nPlot)

	t <- 1
	st <- 1
	en <- seqLengths
	for( l in 1:nAlpha) {
		
		al <- ifelse( (nAlpha == 1), "", paste("a=", df$alpha[st], ""))

		# disegna prima il plot relativo al test T1
		with(df, lines(len[st:en], T1[st:en],  type='b', lty=linetype[l], col=colors[t], pch=plotchar[t]))
		
		leg[t] <- paste(al, "T1 Check", "")
		cl[t] <- colors[t]

		# for( i in 1:nGamma) {
		#   	t <- ((nGamma + 1) * (l - 1)) + i + 1
		#
		# 	with(df, lines(len[st:en], power[st:en],  type='b', lty=linetype[l], col=colors[t], pch=plotchar[t]))
		#
		#   	leg[t] <- paste(al, "g=", df$gamma[st], "")
		# 	cl[t] <- colors[t]
		#
		# 	st <- st + seqLengths
		# 	en <- seqLengths * (i+1)
		# }
	}
		
	legend( 'topleft', # 200, 1.0, #
  		legend = leg, col = cl, pch = plotchar,  bty = "y",  pt.cex = 0.8, cex = 0.5, text.col = "black", #
		horiz = F, inset = c(0.01, 0.01))
	
	# readline(prompt="Press [enter] to continue")
	dev.off() #only 129kb in size
	print( " ... done")
}