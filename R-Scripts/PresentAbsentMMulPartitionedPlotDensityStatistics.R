library(DescTools)
library(dplyr)
library(ggplot2)
library(hrbrthemes)



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

###### CODE

# Sets the path of the directory containing the output of FADE
setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")

bs <- "uniform"

similarities = c('D2')

# Defines the name of the file containing a copy of the dataframe created by this script

dfFilename <- sprintf( "%s,32/ReportMMul.RDS", bs)
csvFilename <- sprintf("%s,32/ReportMMul.csv", bs)
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

df <-readRDS( file = dfFilename)
cat(sprintf("Dataset %s loaded. (%d rows).\n", dfFilename, nrow(df)))

classes <- c("rare", "normal", "saturated")
thresholds = c(0.01, 0.99)

df$kf = factor(df$k)
df$tf = factor(df$Theta)

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

tgtDF <- data.frame( Measure = character(), Theta = integer(), k = integer(),
                     A = numeric(), B = numeric(), C = numeric(), D = numeric(), N = numeric(), density = numeric(),
                     distance=double(), stringsAsFactors=FALSE)

measures = colnames(df)[13:27]
extras = c("Mash.Distance.10000.", "D2", "Euclidean")
measures = append(measures, extras)

for(i in 1:nrow(df)) {
  r <- df[i,]
  for(m in measures) {
    MesName <- if (m == "Mash.Distance.10000.") "Mash" else m
    nr <- c( MesName, r[5:12], df[i, m])
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
# tgtDF$Theta = factor(tgtDF$Theta)

kValues = levels(factor(tgtDF$kf))
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
sp <- ggplot(data=df, aes(x=Theta, y=A.N, label=A.N)) +
  geom_line(aes(color = kf)) +
  geom_point() +
  geom_text(aes(label = round(A.N, 2)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
  facet_grid( cols = vars(class), rows = vars(k), labeller = labeller( k = label_both)) +
  scale_y_continuous(name = "A/N",
                     breaks=c(0, 0.5, 1),
                     labels=c("0", "0.5", "1")) +
  theme_light() + theme( panel.spacing=unit(0.2, "lines"),
                         legend.position = "none")

outfname <- sprintf( "%s/PanelMMulDensity.pdf", dirname)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off()
totPrinted <- 1


# secondo grafico distanze per ciascuna misura (Theta sull'asse delle x)
df2 = filter( tgtDF, Measure != "D2" & Measure != "Euclidean")

sp1 <- ggplot( df2, aes(x = Theta, y = distance, fill = k)) +
    geom_line(aes(color = k)) +
    geom_point( size = 0.8) +
    #  geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
    facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) +
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                           axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                           panel.spacing=unit(0.1, "lines")) +
    theme(legend.position = "none") +
    scale_y_continuous(name = "Distance",
                     breaks=c(0, 0.5, 1),
                     labels=c("0", "0.5", "1"))
    # labs(y = "Distance") +
    # guides(colour = guide_legend(override.aes = list(size=1)))


# dev.new(width = 6, height = 6)
# print(sp1)
outfname <- sprintf( "%s/PanelDistanceMMul.pdf", dirname)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1

# terzo grafico di riferimento per Euclid
df2 = filter( tgtDF, Measure == "Euclidean")

sp1 <- ggplot( df2, aes(x = Theta, y = distance)) +
  # geom_bar( width = 0.7, position = "dodge", stat = "identity") +
  geom_line(aes(color = k)) +
  geom_point() +
  # geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
  facet_grid( rows = vars(k), scales = "free_y", labeller = labeller( k = label_both)) +
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
outfname <- sprintf( "%s/PanelEuclidDistanceMMul.pdf", dirname)
ggsave( outfname, device = pdf(), width = 6, height = 9, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1

# Quarto grafico di riferimento per simmilarità D2
df2 = filter( tgtDF, Measure == "D2")

sp1 <- ggplot( df2, aes(x = Theta, y = distance)) +
  # geom_bar( width = 0.7, position = "dodge", stat = "identity") +
  geom_line(aes(color = k)) +
  geom_point() +
  # geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
  facet_grid( rows = vars(k), scales = "free_y", labeller = labeller( k = label_both)) +
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
outfname <- sprintf( "%s/PanelD2DistanceMMul.pdf", dirname)
ggsave( outfname, device = pdf(), width = 6, height = 9, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1


cat(sprintf("MMul Done. %d plot printed", totPrinted))