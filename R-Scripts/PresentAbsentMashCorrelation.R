library(DescTools)
library(dplyr)
library(ggplot2)


###### DESCRIPTION

# compute Power Statistics and T1 error from raw data prduced by PresenAbsent.py script
# in CSV format


###### OPTIONS

# Sets the path of the directory containing the output of FADE
setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")

# Defines the name of the file containing a copy of the dataframe created by this script
dfFilename <- "PresentAbsent-Correlation.RDS"
csvFilename <- 'PresentAbsentData-all.csv'
dirname <- "PlotAN"

###### CODE

# if (file.exists(dfFilename)) {
#   cat( sprintf("Data file %s exists. Do you want to overwrite (Y/N) ?\n", dfFilename))
#   res <- readline()
#   if (res != "Y") {
#     quit(save="ask")
#   }
# }


# if (!file.exists(dfFilename)) {
if ( TRUE) {
  # converte il file CSV in dataframe
  
  columnClasses = c(
                #   model	    gamma	    seqLen	   pairId	       k	
                "character", "numeric", "integer", "integer", "integer",
                #   A	        B	         C	        D	        N
                "numeric", "numeric", "numeric", "numeric", "numeric",
                # 15 x misure present absent
                # Anderberg	Antidice	 Dice	     Gower	    Hamman	  Hamming	   Jaccard	  Kulczynski
                "NULL",   "NULL",      "NULL",  "NULL",    "NULL",   "NULL",     "NULL",    "NULL",
                # Matching	 Ochiai	     Phi	     Russel	   Sneath    	Tanimoto	  Yule
                "NULL",    "NULL",    "NULL",     "NULL",  "NULL",     "NULL",    "NULL",
                # mash 4 x 3
                # Pvalue,   distance,     A,        N
                "numeric", "numeric", "numeric", "numeric", 
                "numeric", "numeric", "numeric", "numeric",
                "numeric", "numeric", "numeric", "numeric",
                # dati entropia 5 x seq x 2 (A-B)
                # sequence-A
                # NKeysA	 2*totalCntA  deltaA	    HkA	       errorA
                "NULL",     "NULL",    "NULL",     "NULL",    "NULL",
                # sequence-B
                # NKeysB	 2*totalCntB  deltaB	    HkB	      errorB
                "NULL",     "NULL",    "NULL",     "NULL",    "NULL")

  dati <-read.csv(file = csvFilename, header = TRUE, sep = ",", colClasses = columnClasses)
  dati$model = factor(dati$model)
  
  ll = levels(factor(dati$gamma))
  gValues = as.double(ll[2:length(ll)])
  kValues = as.integer(levels(factor(dati$k)))
  lengths = as.integer(levels(factor(dati$seqLen)))
  col <- colnames(dati)

  ll = levels(factor(dati$gamma))
  gValues = as.double(ll[2:length(ll)])
  kValues = as.integer(levels(factor(dati$k)))
  lengths = as.integer(levels(factor(dati$seqLen)))
  col <- colnames(dati)
  models = levels(dati$model)
  
  corDF <- data.frame( model = character(),  gamma = double(), seqLen = integer(), k = integer(),
                      sks = integer(), correlation = numeric(), stringsAsFactors=FALSE)
  
  for(mm in models) {
    gammas = if (mm == 'Uniform' || mm == 'Uniform-T1') c(0) else c(0.01, 0.05, 0.10)
    for(g in gammas) {
      for( len in lengths) {
        for( k in kValues) {
          kv <- as.integer(k)
          df <- filter( dati, dati$model == mm & dati$seqLen == len & dati$gamma == g & dati$k == kv)
          # c <- cor.test(df$A/(df$N-df$D), df$A..1000./df$N..1000.,  method = "spearman")
          c <- cor.test(df$A, df$A..1000.,  method = "spearman")
          r1 <- as.numeric(c[4])
          corDF[nrow( corDF)+1,] <- list(mm, g, len, kv, 1000, r1)
          
          # c <- cor.test(df$A/(df$N-df$D), df$A..10000./df$N..10000.,  method = "spearman")
          c <- cor.test(df$A, df$A..10000.,  method = "spearman")
          r2 <- as.numeric(c[4])
          corDF[nrow( corDF)+1,] <- list(mm, g, len, kv, 10000, r2)
          

          # c <- cor.test(df$A/(df$N-df$D), df$A..100000./df$N..100000.,  method = "spearman")
          c <- cor.test(df$A, df$A..100000.,  method = "spearman")
          r3 <- as.numeric(c[4])
          corDF[nrow( corDF)+1,] <- list(mm, g, len, kv, 100000, r3)
          
          cat( sprintf("%s(%.3f), len = %d, k = %d -> %s, %s, %s\n", mm, g, len, kv, r1, r2, r3))
        }
      }
    }
  }
  
  corDF$model = factor(corDF$model)
  saveRDS(corDF, file = dfFilename)
  cat(sprintf("Dataset %s %d rows saved.", dfFilename, nrow(corDF)))

} else {
  # carica il dataframe esistente
  corDF <- readRDS(file = dfFilename)
}

corDF$kf = factor(corDF$k)
corDF$lf = factor(corDF$seqLen)

nm <- filter(corDF, corDF$model == 'Uniform') 
corDF <- filter(corDF, corDF$model != 'Uniform' & corDF$model != 'Uniform-T1') # solo gli AM per ogni gamma

nm$gamma <- 0.01
corDF <- rbind(corDF, nm)

nm$gamma <- 0.05
corDF <- rbind(corDF, nm)

nm$gamma <- 0.10
corDF <- rbind(corDF, nm)

cat(nrow(corDF))

  
for (kvf in levels(factor(corDF$k))) {
  kv = as.integer( kvf)
  dff <- filter(corDF, corDF$k == kv) # seleziona solo i valori k correnti
  
  cat(sprintf("k = %d -> %d rows\n", kv, nrow(dff)))
  
  sp <- ggplot( dff, aes(x = lf, y = correlation, group = model)) +
    geom_line( aes(color = model)) +
    geom_point( aes(color = model)) +
    facet_grid(rows = vars(gamma), cols = vars(sks)) +
    # facet_grid_sc(rows = vars( gamma), scales = 'free') +
    # scale_x_continuous(name = NULL, breaks=c(1000, 10000, 100000, 1000000, 10000000),
    #                   labels=c("", "1e+4", "", "1e+6", ""), limits = c(1000, 10000000), trans='log10') +
    scale_y_continuous(name = "Spearman Correlation A vs A(mash))") +
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 70),
                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                          panel.spacing=unit(0.1, "lines")) +
    guides(colour = guide_legend(override.aes = list(size=3)))
  # ggtitle( am)
  
  # dev.new(width = 9, height = 6)
  # print(sp)
  outfname <- sprintf( "%s/CorrelationAN-k=%d.pdf", dirname, kv)
  ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
  dev.off() #only 129kb in size
}
