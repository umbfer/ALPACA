library(rjson)
setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset2-500")

files <- list.files("json", "*.json")

print( paste("Processing", length(files), "files from:", getwd()))

for (file in files) {
	
	cat( "Processing: ", file)
	
    exp <- fromJSON(file= paste("json", file, sep="/"))
    values <-exp[['values']]

    df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=5, byrow=T))})
    df <- do.call(rbind, df)

	outfname <- paste( "img/", tools::file_path_sans_ext(file), ".png", sep = "")
    png( outfname, width=15, height=10, units="cm", res=300)
    #par(mar=c(4,4,1,1))

	title = paste(gsub( "[dD]istance", "", exp$header$distanceName), exp$header$alternateModel, exp$header$nullModel, sep = "-")
    with( df, plot( X1[1:15], X4[1:15], type='b', pch=1,
  		xlab = "Sequence Length", xlim = c(900, 15000000), ylab = "Power Value", ylim = c(0,1), log="x",
  		main = title )) # sub = "k = 6"
  	grid( 8, 5,col = "lightgray", lty = "dotted")
	with(df, lines(X1[16:30],  X4[16:30],  type='b', col='blue', pch=2))
	with(df, lines(X1[31:45], X4[31:45], type='b', col='green', pch=3))
	
	legend( 'topleft', # 200, 1.0, #
  		legend = c("g = 1%", "g = 5%", "g = 10%"), #
		col = c('black', 'blue', 'green'), #
		pch = c(1,2,3),  bty = "y",  pt.cex = 0.8, cex = 0.6, text.col = "black", #
		horiz = F, inset = c(0.01, 0.01))
	
	# readline(prompt="Press [enter] to continue")
	dev.off() #only 129kb in size
	print( " ... done")
}