library(ggplot2)

setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov/AnalisiAlternativeModel/seqs/subSeqs")
#setwd("/Volumes/Catalina/PowerStatistics/Kolmogorov/")

k=4
model='MotifRepl-U'
model='PatTransf-U'
nPairs=1000
len=2000000

g=0.10

pair <- 1

model <- 'PatTransf-U'

for( g in c(0.01, 0.05, 0.10)) {
	
	for(s in 1:10) {
		
		gVal <- if (model == "Uniform") "" else sprintf(".G=%.3f", g)
		
		# f1 <- sprintf("distk=%d_%s-%d.%d%s-s1.1-%d.dist", k, 'Uniform', nPairs, len, '', s)
		f2 <- sprintf("distk=%d_%s-%d.%d%s-s1.1-%d.dist", k, model, nPairs, len, gVal, s)
		f3 <- sprintf("distk=%d_%s-%d.%d%s-s1.2-%d.dist", k, 'Uniform', nPairs, len, '', s)
		f4 <- sprintf("distk=%d_%s-%d.%d%s-s1.2-%d.dist", k, model, nPairs, len, gVal, s)

		df1 <- data.frame()
		# df1 = read.delim( f1, header = FALSE, stringsAsFactors = FALSE, quote = "", sep = "\t")
		# colnames(df1) <- c('kmer', 'values')
		# df1['seq'] <- 's1'
	
		df2 = read.delim( f2, header = FALSE, stringsAsFactors = FALSE, quote = "", sep = "\t")
		colnames(df2) <- c('kmer', 'values')
		df2['seq'] <- 's1PT'
		df1 <- rbind(df1, df2)

		df3 = read.delim( f3, header = FALSE, stringsAsFactors = FALSE, quote = "", sep = "\t")
		colnames(df3) <- c('kmer', 'values')
		df3['seq'] <- 's2'
		df1 <- rbind(df1, df3)

		df4 = read.delim( f4, header = FALSE, stringsAsFactors = FALSE, quote = "", sep = "\t")
		colnames(df4) <- c('kmer', 'values')
		df4['seq'] <- 's2PT'
		df1 <- rbind(df1, df4)
		
		df1$kmer <- as.factor(df1$kmer)
	
		title = sprintf("model:%s, k = %d, g = %.2f, subseq = %s", model, k, g, s)
		sp2 <-	ggplot( df1, aes( x = kmer, y = values, fill = seq)) +
			ggtitle( title) + labs( x = "k-mers", y = "k-mer probability") +
			theme(axis.text.x = element_text(size = 6, angle = 45, vjust = 0.95, hjust=1)) +
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
			outfname <- sprintf( "%s.png",  tools::file_path_sans_ext(f4))
			ggsave( outfname, device = png(), width = 9, height = 6, dpi = 300)
			# ggsave( outfname)
			print( sp2)
			# stop("break")
			# readline(prompt="Press [enter] to continue")
			dev.off() #only 129kb in size
			cat( sprintf("%d ok ", s))	
	} # per ogni sottosequenza
	cat( sprintf("%s%s.  done.\n", model, gVal))
} # per ogni gamma