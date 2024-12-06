library(ggplot2)

setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov/AnalisiAlternativeModel/")

# produce il plot a barre di due sequenze per ciascun AM


k=4
model='MotifRepl-U'
model='PatTransf-U'
nPairs=1000
len=2000000

g=0.10

pair <- 1

for( g in c(0.01, 0.05, 0.10)) {
	
	for(model in c('MotifRepl-U', 'PatTransf-U', 'Uniform')) {
		
		gVal <- if (model == "Uniform") "" else sprintf(".G=%.3f", g)
		
		f1 <- sprintf("k=%d/dist-k=%d_%s-%d.%d%s-All.distbin", k, k, model, nPairs, len, gVal)
		inFile <- file( f1, "rb")
		nKeys = readBin( inFile, integer(), n = 1)
	
		xv <- readBin( inFile, double(), n = nKeys)
		df1 <- as.data.frame( xv)
		colnames(df1) <- c('values')
		df1['seq'] <- 's1'
		df1['pos'] <- 1:nKeys
	
		yv <- readBin( inFile, double(), n = nKeys)
		df2 <- as.data.frame( yv)
		colnames(df2) <- c('values')
		df2['seq'] <- 's2'
		df2['pos'] <- 1:nKeys
		df1 <- rbind(df1, df2)
	#}	
	
		title = sprintf("model:%s, k = %d, g = %.2f", model, k, g)
		sp2 <-	ggplot( df1, aes( x = pos, y = values, fill = seq)) +
			ggtitle( title) + labs( x = "k-mers", y = "k-mer probability") +
		    geom_bar(stat="identity", position=position_dodge())
			
	#			labs(color = "Model") +
				# theme(axis.text.x = element_blank()) + # axis.ticks.x = element_blank()) +
	#			theme_light(base_size = 10) + 
				# theme(legend.position = "bottom") +
	#			theme(legend.position = "none") + 
	#			theme(axis.text.x = element_text(size = 7, angle = 45, vjust = 0.95, hjust=1)) +
	#			scale_colour_brewer(palette="Dark2") +
	#			scale_y_continuous(name = "k-mes probability", limits = c(0, fs)) +
	#			geom_boxplot(lwd = 0.2, alpha = 0.9) # width=15, outlier.shape=20
	
			# dev.new(width = 9, height = 6)
			outfname <- sprintf( "%s-Histogram-k=%d-len=%d%s.png", model, k, len, gVal)
			ggsave( outfname, device = png(), width = 9, height = 6, dpi = 300)
			# ggsave( outfname)
			# print( sp2)
			# stop("break")
			# readline(prompt="Press [enter] to continue")
			# dev.off() #only 129kb in size
			cat( sprintf("%s%s.  done.\n", model, gVal))
	
	} # per ogni modello
} # per ogni gamma