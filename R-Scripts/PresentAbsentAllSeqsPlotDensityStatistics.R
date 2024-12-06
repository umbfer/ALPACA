library(DescTools)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(r2r)


###### DESCRIPTION

# Plot experiment results for Macaco Mulatta Genomes artificially modified


###### OPTIONS
###### CODE
plot_labeller <- function(variable, value){
  # cat(sprintf("variable: <%s>, value: <%s>\n", variable, as.character(value)))
  if (variable == 'kv') {
    # N.B. kv e' un factor
    return(sprintf("k = %s", as.character(value)))
  } else if (variable == 'k') {
    return(sprintf("k = %d", value))
  }  else if (variable == 'lenFac') {
    # lenFac è un factor
    return(formatC(as.numeric(as.character(value)), format="f", digits=0, big.mark="."))
  }else {
    return(as.character(value))
  }
}


TranslationTable  <- hashmap(default = 0)
TranslationTable[["Mash.Distance.10000."]] <- "Mash"
TranslationTable[["GCF_0001654452_Mmur_30"]] <- "gray mouse lemur"
TranslationTable[["GCF_0033397651_Mmul_10"]] <- "rhesus monkey"
TranslationTable[["GCF_0125594852_MFA2"]] <- "crab eating macaque"
TranslationTable[["GCF_0009559451_Caty_10"]] <- "sooty mangabey"

TerminologyServer <- function( key) {
  v = TranslationTable[[key]]
  if (v == 0) {
    return( key)
  } else {
    return( v)
  }
}

###### CODE

# Sets the path of the directory containing the output of FADE
setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")

bs <- "uniform"

similarities = c('D2')

# Defines the name of the file containing a copy of the dataframe created by this script

dfFilename <- sprintf( "%s,32/AllRealSeq2.RDS", bs)
csvFilename <- sprintf("%s,32/AllRealSeq2.csv", bs)
dirname <- sprintf("%s,32/ReportMMul", bs)

if (!dir.exists(dirname)) {
  dir.create(dirname)
}

if (!file.exists(dfFilename)) {
  # carica il CSV dell'esperimento
  columnClasses = c(
    #   sequenceA  sequenceB  start.time  real.time    Theta        k
    "character", "character", "integer", "numeric", "integer", "integer",
    #   A	        B	      C	         D	        N         A/N
    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
    # 15 x misure present absent
    # Anderberg	Antidice	 Dice	     Gower	    Hamman	  Hamming	   Jaccard	  Kulczynski
    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
    # Matching	 Ochiai	     Phi	     Russel	   Sneath    	Tanimoto	  Yule
    "numeric", "numeric", "numeric", "numeric", "numeric", "numeric", "numeric",
    # mash 4 x 1 (P value, Mash distance, A, N size = 10.000)
    "numeric", "numeric", "numeric", "numeric",
    # D2     Euclidean    Euclid_norm
    "numeric", "numeric", "numeric",
    # NKeysA   totalCntA   deltaA       HkA      errorA
    "numeric", "numeric", "numeric", "numeric","numeric",
    # NKeysB   totalCntB   deltaB       HkB      errorB
    "numeric", "numeric", "numeric", "numeric","numeric")

  df <-read.csv( file = csvFilename, sep = ";", dec = ",", colClasses = columnClasses)
  saveRDS( df, file = dfFilename)
  cat(sprintf("Dataset %s %d rows saved.", dfFilename, nrow(df)))
} else {
  # dataframe già disponibile
  df <-readRDS( file = dfFilename)
}

# 'data.frame':	48 obs. of  44 variables:
# $ sequenceA           : Factor w/ 1 level "GCF_003339765.1_Mmul_1.0": 1 1 1 1 1 1 1 1 1 1 ...
# $ sequenceB           : Factor w/ 6 levels "GCF_003339765.1_Mmul_1.0-10",..: 3 1 2 4 5 6 3 1 2 4 ...
# $ start.time          : int  1702206937 1702236286 1702268144 1702300312 1702338382 1702370789 1702207380 1702236826 1702268569 1702300989 ...
# $ real.time           : Factor w/ 48 levels "222,4092066",..: 17 23 14 30 38 29 4 12 5 21 ...
# $ Theta               : int  5 10 20 80 90 95 5 10 20 80 ...
# $ k                   : int  4 4 4 4 4 4 8 8 8 8 ...
# $ A                   : Factor w/ 33 levels "0","103","108.466.454",..: 13 13 13 13 13 13 25 25 25 25 ...
# $ B                   : Factor w/ 31 levels "0","1.390.399.430",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ C                   : Factor w/ 33 levels "0","1.039.140.839",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ D                   : Factor w/ 29 levels "0","1.094.310.351.548",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ N                   : Factor w/ 8 levels "1.099.511.627.776",..: 4 4 4 4 4 4 7 7 7 7 ...
# $ A.N                 : Factor w/ 32 levels "0","0,000324535",..: 13 13 13 13 13 13 13 13 13 13 ...
# $ Anderberg           : Factor w/ 30 levels "0,213759018",..: 30 30 30 30 30 30 30 30 30 30 ...
# $ Antidice            : Factor w/ 32 levels "0","0,016613598",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Dice                : Factor w/ 32 levels "0","0,004205805",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Gower               : Factor w/ 31 levels "-0,3132203","-0,962032841",..: 31 31 31 31 31 31 31 31 31 31 ...
# $ Hamman              : Factor w/ 27 levels "0","0,010904372",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Hamming             : Factor w/ 27 levels "0","0,002733565",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Jaccard             : Factor w/ 32 levels "0","0,00837638",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Kulczynski          : Factor w/ 32 levels "0","0,004189791",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Matching            : Factor w/ 27 levels "0","0,002733565",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Ochiai              : Factor w/ 32 levels "0","0,004197798",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Phi                 : Factor w/ 30 levels "-1,03E+33","-1,42E+35",..: 30 30 30 30 30 30 30 30 30 30 ...
# $ Russel              : Factor w/ 22 levels "0","0,008593202",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Sneath              : Factor w/ 27 levels "0","0,001368653",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Tanimoto            : Factor w/ 28 levels "0","0,005452227",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Yule                : Factor w/ 31 levels "0","0,000340283",..: 24 24 24 24 24 24 24 24 24 24 ...
# $ Mash.Pv..10000.     : Factor w/ 25 levels "0","0,0148949",..: 1 1 1 1 1 1 1 1 1 1 ...
# $ Mash.Distance.10000.: Factor w/ 16 levels "0","0,0348505",..: 7 7 7 7 7 7 6 6 6 6 ...
# $ A..10000.           : int  136 136 136 136 136 136 10000 10000 10000 10000 ...
# $ N..10000.           : int  136 136 136 136 136 136 10000 10000 10000 10000 ...
# $ D2                  : Factor w/ 46 levels "0","1,00952E+12",..: 34 33 31 26 25 24 23 19 17 8 ...
# $ Euclidean           : Factor w/ 42 levels "1.832.067,52",..: 7 24 36 4 5 6 30 38 42 8 ...
# $ Euclid_norm         : int  0 0 0 0 0 0 0 0 0 0 ...
# $ NKeysA              : num  256 256 256 256 256 ...
# $ X2.totalCntA        : num  5.87e+09 5.87e+09 5.87e+09 5.87e+09 5.87e+09 ...
# $ deltaA              : Factor w/ 8 levels "0,002831776",..: 8 8 8 8 8 8 7 7 7 7 ...
# $ HkA                 : Factor w/ 8 levels "-15,28220155",..: 8 8 8 8 8 8 1 1 1 1 ...
# $ errorA              : Factor w/ 8 levels "-0,000127184",..: 7 7 7 7 7 7 8 8 8 8 ...
# $ NKeysB              : num  256 256 256 256 256 ...
# $ X2.totalCntB        : num  5.87e+09 5.87e+09 5.87e+09 5.87e+09 5.87e+09 ...
# $ deltaB              : Factor w/ 32 levels "0,002854579",..: 32 32 32 32 32 32 31 31 31 31 ...
# $ HkB                 : Factor w/ 47 levels "-15,47734591",..: 42 43 44 47 46 45 1 2 3 6 ...
# $ errorB              : Factor w/ 42 levels "-0,00011904",..: 37 36 35 34 34 34 42 41 40 38 ...

cat(sprintf("Dataset %s loaded. (%d rows).\n", dfFilename, nrow(df)))

classes <- c("rare", "normal", "saturated")
thresholds = c(0.01, 0.99)


for(i in 1:nrow(df)) {
  t <- df[i, 'A.N']
  if (t < thresholds[1])        { lbl <- classes[1]}
  else if (t > thresholds[2])   { lbl <- classes[3]}
  else                          { lbl <- classes[2]}
  df[i, 'class'] <- lbl
}
# riordina per classi
df$class <- factor(df$class, levels = classes)

# escludiamo le misure: euclidean norm, anderberg, gowel , phi e yule. 15120 -> 10800 observations
# df <- filter( df, Measure != "Anderberg" & Measure != "Gower" & Measure != "Phi" & Measure != "Yule" &
#  Measure != "Euclid_norm" & Measure != "Mash.Distance.1000." & Measure != "Mash.Distance.100000.")
#
# cat(sprintf("Filtered measures: Anderberg, Gower, Phi, Yule, Euclid_norm, Mash.Distance.1000, Mash.Distance.100000. (%d rows).\n",nrow(df)))
# in questo caso non interessa il valore di theta (sempre 0)
tgtDF <- data.frame( Pair = character(), Measure = character(), k = integer(),
                     A = numeric(), B = numeric(), C = numeric(), D = numeric(), N = numeric(), density = numeric(),
                     distance=double(), stringsAsFactors=FALSE)

measures = colnames(df)[13:27]
extras = c("Mash.Distance.10000.", "D2", "Euclidean")
measures = append(measures, extras)

for(i in 1:nrow(df)) {
  r <- df[i,]
  for(m in measures) {
    nr <- c( sprintf( "%s/%s", TerminologyServer(r[[1]]), TerminologyServer(r[[2]])), TerminologyServer(m),
             r[6:12], df[i, m]) #theta non interessa
    tgtDF[nrow( tgtDF)+1,] <- nr
  }
}

# escludiamo le misure: euclidean norm, anderberg, gowel , phi e yule. 864 -> 672 observations
tgtDF <- filter( tgtDF, Measure != "Anderberg" & Measure != "Gower" & Measure != "Phi" & Measure != "Yule" &
  Measure != "Euclid_norm" & Measure != "Mash.Distance.1000." & Measure != "Mash.Distance.100000.")

# nessun filtro 864 -> 864
# tgtDF <- filter( tgtDF, Measure != "Euclid_norm" & Measure != "Mash.Distance.1000." & Measure != "Mash.Distance.100000.")

cat(sprintf("Filtered measures: Anderberg, Gower, Phi, Yule, Euclid_norm, Mash.Distance.1000, Mash.Distance.100000. (%d rows).\n",nrow(df)))

tgtDF$k = factor(tgtDF$k)


tgtDF$Pair <- factor(tgtDF$Pair, levels = c( "rhesus monkey/crab eating macaque",
                                             "gray mouse lemur/crab eating macaque",
                                             "gray mouse lemur/rhesus monkey"))

kValues = levels(factor(tgtDF$k))
measures <- levels(factor(tgtDF$Measure))

for(i in 1:nrow(tgtDF)) {
  t <- tgtDF[i, 'density']
  if (t < thresholds[1])        { lbl <- classes[1]}
  else if (t > thresholds[2])   { lbl <- classes[3]}
  else                          { lbl <- classes[2]}
  tgtDF[i, 'class'] <- lbl
}

# riordina per classi
tgtDF$class <- factor(tgtDF$class, levels = classes)

cat(sprintf("Data Frame converted (%d observations).\n", nrow(tgtDF)))

#
# primo grafico solo della density
df2 = filter( tgtDF, Measure == "Jaccard")

sp1 <- ggplot(data=df2, aes(x=Pair, y=density, label=density)) +
  geom_line(aes(color = k)) +
  geom_point() +
  geom_text(aes(label = round(density, 2)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
  facet_grid( cols = vars(class), rows = vars(k), labeller = labeller( k = label_both)) +
  scale_y_continuous(name = "A/N",
                     breaks=c(0, 0.5, 1),
                     labels=c("0", "0.5", "1")) +
  theme_light() + theme( panel.spacing=unit(0.2, "lines"),
                         legend.position = "none",
                         strip.text.x = element_text( size = 8, angle = 0),
                         axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                         axis.title.x = element_blank())
#                         axis.text.x=element_blank())

outfname <- sprintf( "%s/PanelRealGenomeDensities.pdf", dirname)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off()
totPrinted <- 1


# secondo grafico distanze per ciascuna misura (Theta sull'asse delle x)
df2 = filter( tgtDF, Measure != "D2" & Measure != "Euclidean")

sp1 <- ggplot( df2, aes(x = Pair, y = distance, fill = Pair)) +
    geom_bar( aes(color = Pair), width = 0.7, position = "dodge", stat = "identity") +
    facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) +
    # theme_light() + theme( panel.spacing=unit(0.1, "lines"),
    #                       legend.position = "none",
    #                       axis.text.x=element_blank()) +
    theme_light() + theme( strip.text.x = element_text( size = 8, angle = 0),
                           axis.text.x = element_blank(), # element_text( size = rel( 0.7), angle = 45, hjust=1),
                           panel.spacing=unit(0.1, "lines"),
                           axis.title.x = element_blank()) +
    scale_y_continuous(name = "Distance",
                     breaks=c(0, 0.5, 1),
                     labels=c("0", "0.5", "1")) +
    guides(colour = guide_legend(override.aes = list(size=1)))


# dev.new(width = 6, height = 6)
# print(sp1)
outfname <- sprintf( "%s/PanelRealGenomeDistances.pdf", dirname)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1

# terzo grafico di riferimento per Euclid
df2 = filter( tgtDF, Measure == "Euclidean")

sp1 <- ggplot( df2, aes(x = Pair, y = distance, fill = Pair)) +
  geom_bar( aes(color = Pair), width = 0.7, position = "dodge", stat = "identity") +
  facet_grid( rows = vars(k), scales = "free_y", labeller = labeller( k = label_both)) +
  theme_light() + theme(axis.text.x=element_blank(),
                        axis.title.x = element_blank()) +
  # theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
  #                       axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
  #                       panel.spacing=unit(0.1, "lines")) +
  # scale_y_continuous(name = "Euclidean Distance")
  #                    breaks=c(0, 0.5, 1),
  #                    labels=c("0", "0.5", "1"))
  labs(y = "Euclidean Distance")
# guides(colour = guide_legend(override.aes = list(size=1)))


# dev.new(width = 6, height = 9)
# print(sp1)
outfname <- sprintf( "%s/PanelRealGenomeEuclidDistance.pdf", dirname)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1

# Quarto grafico di riferimento per simmilarità D2
df2 = filter( tgtDF, Measure == "D2")

sp1 <- ggplot( df2, aes(x = Pair, y = distance, fill = Pair)) +
  geom_bar( aes(color = Pair), width = 0.7, position = "dodge", stat = "identity") +
  facet_grid( rows = vars(k), scales = "free_y", labeller = labeller( k = label_both)) +
  theme_light() + theme(axis.text.x=element_blank(),
                        axis.title.x = element_blank()) +
  # theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
  #                       axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
  #                       panel.spacing=unit(0.1, "lines")) +
  # scale_y_continuous(name = "Euclidean Distance")
  #                    breaks=c(0, 0.5, 1),
  #                    labels=c("0", "0.5", "1"))
  labs(y = "D2 Similarity")
# guides(colour = guide_legend(override.aes = list(size=1)))


# dev.new(width = 6, height = 9)
# print(sp1)
outfname <- sprintf( "%s/PanelRealGenomeD2Distance.pdf", dirname)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1


cat(sprintf("Real genome distances Done. %d plot printed", totPrinted))