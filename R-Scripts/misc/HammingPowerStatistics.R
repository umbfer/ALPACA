library(DescTools)
library(ggplot2)
library(dplyr)


setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov/AnalisiAlternativeModel/seqs")
#setwd("/Volumes/Catalina/PowerStatistics/Kolmogorov/")

nPairs = 1000
lenghts = c(200000, 5000000)

dfFilename = "HammingDistance-All.df"			

if (file.exists(dfFilename)) {
	cat( sprintf("Data file exists. Loading %s\n", dfFilename))			
	# carica il dataframe dal file
	df <- readRDS(file = dfFilename)
} else {

	df <- data.frame( Name = character(), Dist = double(), len = numeric(), stringsAsFactors=FALSE)

	for( len in lenghts)	{
		model = 'Uniform'
		f1 <- sprintf("%s-%d.%d%s.fasta", model, nPairs, len, '')
		con1 = file(f1, "r")
		
		for( pair in 1:nPairs) {
			
			s1Name = readLines(con1, n = 1)
			s1 = readLines(con1, n = 1)
			s2Name = readLines(con1, n = 1)
			s2 = readLines(con1, n = 1)
		
			d <- StrDist(s1, s2, method = 'hamming')
			df[nrow(df) + 1,] <- list( 'NM', d[1]/len, len)
			
			cat( sprintf("%5d /%5d\r", pair, nPairs))
		}
		close(con1)
		cat( sprintf("%s.  done.\n", 'Uniform'))
	
		for( model in c('MotifRepl-U', 'PatTransf-U')) {
			
			for( g in c(0.01, 0.05, 0.10)) {
		
				gVal <- sprintf(".G=%.3f", g)
				k1 <- sprintf("%s%s", if (model == 'MotifRepl-U') 'MR' else 'PT', gVal)
						
				f2 <- sprintf("%s-%d.%d%s.fasta", model, nPairs, len, gVal)
				con2 = file(f2, "r")
			
				for( pair in 1:nPairs) {
		
					s1AMName = readLines(con2, n = 1)
					s1AM = readLines(con2, n = 1)
					s2AMName = readLines(con2, n = 1)
					s2AM = readLines(con2, n = 1)
				
					d <- StrDist(s1AM, s2AM, method = 'hamming')
					df[ nrow(df) + 1,] <- list( k1, d[1]/len, len)
					
					cat( sprintf("%5d /%5d\r", pair, nPairs))
				} # for all pairs
	
				close(con2)
				cat( sprintf("\n%s ok\n", f2))	
			} # for all gamma
			cat( sprintf("%s.  done.\n", model))
		} # for all models
	} # for all lenngth
	saveRDS( df, file = dfFilename)
}

df$Name = factor(df$Name, levels = c('NM', 'MR.G=0.010','MR.G=0.050','MR.G=0.100', 'PT.G=0.010','PT.G=0.050','PT.G=0.100'))


HammingPower <- data.frame( Name = character(), len = numeric(), Power = double(), alpha = double(), threshold = double(), stringsAsFactors=FALSE)

aValues <- c( 0.01, 0.05, 0.10)

for( len in lenghts)	{
	df2 <- filter(df, df$Name == 'NM' & df$len == len)
	NM <- df2[order(df2$Dist),2]
	
	threshold <- NM[length(NM) * aValues[1]]
	threshold <- c(threshold, NM[length(NM) * aValues[2]])
	threshold <- c(threshold, NM[length(NM) * aValues[3]])

	for(i in 1:3) {
		for( model in c('MR', 'PT')) {
			for( g in c(0.01, 0.05, 0.10)) {
				k1 <- sprintf("%s.G=%.3f", model, g)
			
				AM <- filter(df, df$Name == k1 & df$len == len)[,2]
				pwr <- 0
				for(dist in AM) {
					if (dist <= threshold[i])
						pwr <- pwr + 1
				}
				HammingPower[nrow( HammingPower) + 1,] <- list( k1, len, pwr / length(AM), aValues[i], threshold[i])
			}
		}
	}
}

HammingPower$Name <- factor(HammingPower$Name)
HammingPower$len <- factor(HammingPower$len)
levels(HammingPower$len) <- c('n = 200 000', 'n = 5 000 000')
	
sp <- ggplot( HammingPower, aes(x = Name, y = Power)) +
# sp <- ggplot( HammingPower, aes(x = Name, y = Power, fill = Name)) + 
 	# geom_boxplot( aes(color = Name), outlier.size = 0.3) +
 	geom_point( aes(color = Name), shape = 15, size = 4) +
 	facet_grid(cols = vars( len)) +
 	scale_y_continuous(name = "HD Power Statistics", limits = c(0, 1)) +
 	theme_bw() + theme( axis.text.x = element_text(size = 10, angle = 45, hjust =1)) +
 	theme(legend.position = "none") + labs(x = "")
 	# ggtitle("Power Statistics for Hamming Distance") + labs(y = "Power")


# dev.new(width = 6, height = 9)
outfname <- sprintf( "HammingPowerStatistics.png",  tools::file_path_sans_ext(dfFilename))
ggsave( outfname, device = png(), width = 9, height = 4, dpi = 300)
# ggsave( outfname, device = png(), dpi = 300)
print( sp)

# stop("break")
# readline(prompt="Press [enter] to continue")
#			dev.off() #only 129kb in size
