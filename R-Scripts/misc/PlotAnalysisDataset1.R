library(rjson)
setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset1")
getwd()


files <- list.files(".", "*.json")

for (file in files) {
	
	cat( "Processing: ", file)
	
    exp <- fromJSON(file=file)
    values <-exp[['values']]
    
    df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=4, byrow=T))})
    df <- do.call(rbind, df)
    
	outfname <- paste( "imgs/", tools::file_path_sans_ext(file), ".png", sep = "")
    # png( outfname, width=15, height=10, units="cm", res=300)
    # par(mar=c(4,4,1,1))

	title = paste(gsub( "[dD]istance", "", exp$header$distanceName), exp$header$alternateModel, exp$header$nullModel, sep = "-")
    with( df, plot( X1[1:8], X3[1:8], type='b', pch=1,
  		xlab = "Sequence Length", xlim = c(200, 30000), ylab = "Power Value", ylim = c(0,1), log="x",
  		main = title, sub = "k = 5"))
  	grid( 6, 5,col = "lightgray", lty = "dotted")
	with(df, lines(X1[9:16],  X3[9:16],  type='b', col='blue', pch=2))
	with(df, lines(X1[17:24], X3[17:24], type='b', col='green', pch=3))
	with(df, lines(X1[25:32], X3[25:32], type='b', col='red', pch=4))
	with(df, lines(X1[33:40], X3[33:40], type='b', col='orange', pch=5))
	
	legend( 'topleft', # 200, 1.0, #
  		legend = c("g = 0.001", "g = 0.005", "g = 0.010", "g = 0.050", "g = 0.100"), #
		col = c('black', 'blue', 'green', 'red', 'orange'), #
		pch = c(1,2,3,4,5),  bty = "y",  pt.cex = 0.8,   cex = 0.6, text.col = "black", #
		horiz = F, inset = c(0.01, 0.01))
	
	# dev.off() #only 129kb in size

	readline(prompt="Press [enter] to continue")
}