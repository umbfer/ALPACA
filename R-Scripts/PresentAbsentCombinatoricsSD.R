library(ggplot2)
library(dplyr)

###### DESCRIPTION

# compute Power Statistics and T1 error from raw data prduced by PresenAbsent.py script
# in CSV format


###### OPTIONS


###### CODE
plot_labeller <- function(variable,value){
  # cat(sprintf("variable: <%s>, value: <%s>\n", variable, as.character(value)))
  if (variable == 'len') {
    # N.B. len e' un factor
    return(len_names[as.character(value)])
  } else if (variable == 'k') {
    return(sprintf("k = %s", as.character(value)))
  }  else {
    return(as.character(value))
  }
}


columnClasses = c(
  #   model	    gamma	  seqLen     pairId	     k
  "character", "numeric", "integer", "integer", "integer",
  #   A	        B	         C	        D	        N
  "numeric", "numeric", "numeric", "numeric", "numeric",
  # 15 x misure present absent
  # Anderberg	Antidice	 Dice	     Gower	    Hamman	  Hamming	   Jaccard	  Kulczynski
  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
  # Matching	 Ochiai	     Phi	     Russel	   Sneath    	Tanimoto	  Yule
  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
  # mash 4 x 3 (P value, Mash distance, A, N)
  "numeric", "numeric", "numeric", "numeric",
  "numeric", "numeric", "numeric", "numeric",
  "numeric", "numeric", "numeric", "numeric",
  # dati entropia 5 x seq x 2 (A-B)
  # sequence-A
  # NKeysA	 2*totalCntA  deltaA	    HkA	       errorA
  "numeric", "numeric", "numeric", "numeric", "numeric",
  # sequence-B
  # NKeysB	 2*totalCntB  deltaB	    HkB	      errorB
  "numeric", "numeric", "numeric", "numeric", "numeric")

unstd = c( 1, 16, 32)

columnClasses1 = c(
  #   model	    gamma	  seqLen     pairId	     k
  "character", "numeric", "integer", "integer", "integer",
  #   A	        B	         C	        D	        N
  "numeric", "numeric", "numeric", "numeric", "numeric",
  # 15 x misure present absent
  # Anderberg	Antidice	 Dice	     Gower	    Hamman	  Hamming	   Jaccard	  Kulczynski
  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
  # Matching	 Ochiai	     Phi	     Russel	   Sneath    	Tanimoto	  Yule
  "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
  # mash 4 x 3 (P value, Mash distance, A, N)
  "numeric", "numeric", "numeric", "numeric",
  "numeric", "numeric", "numeric", "numeric",
  "numeric", "numeric", "numeric", "numeric",
  # solo in 1,32 ????
  # D2,     Euclidean, Euclid_norm
  "numeric", "numeric", "numeric",
  # dati entropia 5 x seq x 2 (A-B)
  # sequence-A
  # NKeysA	 2*totalCntA  deltaA	    HkA	       errorA
  "numeric", "numeric", "numeric", "numeric", "numeric",
  # sequence-B
  # NKeysB	 2*totalCntB  deltaB	    HkB	      errorB
  "numeric", "numeric", "numeric", "numeric", "numeric")


geneSize = 1
patternLen = 32

# Sets the path of the directory containing the output of FADE
setwd( sprintf("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent/%d,%d",
 geneSize, patternLen))

# Defines the name of the file containing a copy of the dataframe created by this script
dfFilename <- sprintf("PresentAbsentEC-%d-%d.RDS", geneSize, patternLen)
csvFilename <- sprintf("PresentAbsentECData-%d-%d.csv", geneSize, patternLen)
nullModel <- 'ShuffledEColi'
T1Model <- paste( sep='', nullModel, '-T1')


if (file.exists(dfFilename)) {
  cat( sprintf("Loading existing data file: %s\n", dfFilename))
  outDF <- readRDS(file = dfFilename)
} else {
  cat( sprintf("Loading data from experimental results from CSV file: %s\n", csvFilename))

  # read experiments raw results
  df <-read.csv( file = csvFilename, colClasses = if (geneSize %in% unstd) columnClasses1 else columnClasses)
  df$model = factor(df$model)
  # df$gamma = factor(df$gamma)
  # df$k = factor(df$k)
  # df$seqLen = factor(df$seqLen)

  ll = levels(factor(df$gamma))
  gValues = as.double(ll[2:length(ll)])
  gValue = 0.05
  kValues = as.integer(levels(factor(df$k)))
  lengths = as.integer(levels(factor(df$seqLen)))
  col <- colnames(df)
  measures <- c(col[11:25], col[27], col[31], col[35])
  altModels = levels(df$model)[1:2]

  outDF <- data.frame(Model = character(), len = numeric(), k = numeric(), A = numeric(), N = numeric(),
                      sqm = numeric(), jaccard = numeric(), normEuclidean = numeric(), stringsAsFactors=FALSE)

  for( modello in levels(df$model)[1:3]) {
    for( len in lengths) {
      for(kv in kValues)  {
        cat(sprintf("Modello = %s, len = %d, k = %d, ", modello, len, kv))
        if (modello == "ShuffledEColi") {
          nm <- filter( df, df$model == modello & df$seqLen == len & df$k == kv)
          modKey = modello
        } else {
          nm <- filter( df, df$model == modello & df$seqLen == len & df$k == kv & df$gamma == gValue)
          modKey = sprintf("%s-g=%.2f", modello, gValue)
        }
        if (nrow(nm) != 1000) { # 3 valori di gamma
          stop("errore estrazione dati")
        }
        nr = nrow(outDF) + 1
        outDF[nr,] <- list(modKey, len, kv, mean(nm$A), mean(nm$N), sd(nm$A/nm$N),
                                      mean(nm$Jaccard), mean(nm$Euclid_norm))
        cat(sprintf("A/N = %f, stdev(A/N) = %f, Jaccard = %f, NormEuclidean = %f\n",
                    outDF[nr,]$A/outDF[nr,]$N, outDF[nr,]$sqm, outDF[nr,]$jaccard, outDF[nr,]$normEuclidean))
      }
    }
  }
  outDF$k = factor(outDF$k)
  outDF$Model = factor(outDF$Model)
  saveRDS(outDF, file = dfFilename)
  cat(sprintf("Dataset %s %d rows saved.\n", dfFilename, nrow(outDF)))
}

title <- sprintf("A/N values for model %s gs=%s pl = %d", levels(outDF$model), geneSize, patternLen)
# dev.new(width = 6, height = 9)
# print(sp1)
# dev.off() #only 129kb in size

# Plot valori A/N
sp1 <- ggplot( outDF, aes(x = len, y = A/N, label = sprintf("%.3f", A/N))) +
  geom_line(aes(color=k)) +
  geom_point() +
  geom_text( hjust=0, vjust=0, size=2) +
  facet_grid( rows = vars( k), cols = vars( Model), labeller = labeller( k = label_both)) +
  theme(legend.position = "none") +
  labs(y = "A/N", x = "len") + scale_x_continuous(trans='log10')

outfname <- sprintf("../Combinatorics/Combinatorics-EC-AN.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)

# Plot valori Deviazione Standard
sp1 <- ggplot( outDF, aes(x = len, y = sqm, label = sprintf("%.5f", sqm))) +
  geom_line(aes(color=k)) +
  geom_point() +
  geom_text( hjust=0, vjust=0, size=1.8) +
  facet_grid( rows = vars( k), cols = vars( Model), labeller = labeller( k = label_both)) +
  theme(legend.position = "none") +
  labs(y = "Std Dev(A/N)", x = "len") + scale_x_continuous(trans='log10')

outfname <- sprintf("../Combinatorics/Combinatorics-EC-SD.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)

# Plot valori distanza euclidea normalizzata
sp1 <- ggplot( outDF, aes(x = len, y = normEuclidean, label = sprintf("%.3f", normEuclidean))) +
  geom_line(aes(color=k)) +
  geom_point() +
  geom_text( hjust=0, vjust=0, size=2) +
  facet_grid( rows = vars( k), cols = vars( Model), labeller = labeller( k = label_both)) +
  theme(legend.position = "none") +
  labs(y = "Normalized Euclidean Distance", x = "len") + scale_x_continuous(trans='log10')

outfname <- sprintf("../Combinatorics/Combinatorics-EC-NormEuclid.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)

# Plot valori distanza di Jaccard
sp1 <- ggplot( outDF, aes(x = len, y = jaccard, label = sprintf("%.3f", jaccard))) +
  geom_line(aes(color=k)) +
  geom_point() +
  geom_text( hjust=0, vjust=0, size=2) +
  facet_grid( rows = vars( k), cols = vars( Model), labeller = labeller( k = label_both)) +
  theme(legend.position = "none") +
  labs(y = "Jaccard Distance", x = "len") + scale_x_continuous(trans='log10')

outfname <- sprintf("../Combinatorics/Combinatorics-EC-Jaccard.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)

# Plot valori distanze NED e Jaccard
sp1 <- ggplot( outDF, aes(x = len)) +
  geom_line(aes(y = normEuclidean, colour="normEuclid")) +
  geom_point(aes(y = normEuclidean)) +
  # geom_text( aes(y = normEuclidean, label = sprintf("%.3f", normEuclidean)), hjust=0, vjust=0, size=2) +
  geom_line(aes(y = jaccard, colour="jaccard")) +
  geom_point(aes(y = jaccard)) +
  # geom_text( aes(y = jaccard, label = sprintf("%.3f", jaccard)), hjust=0, vjust=0, size=2) +
  facet_grid( rows = vars( k), cols = vars( Model), labeller = labeller( k = label_both)) +
  # theme(legend.position = "none") +
  labs(y = "Jaccard vs Norm Euclidean Distances", x = "len") + scale_x_continuous(trans='log10')

outfname <- sprintf("../Combinatorics/Combinatorics-EC-Distances.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)
