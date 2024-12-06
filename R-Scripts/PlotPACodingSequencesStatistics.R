library(DescTools)
library(dplyr)
library(ggplot2)
library(hrbrthemes)



###### DESCRIPTION

# Plot experiment results for real genome sequences compared against the same sequence artificially modified


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
  } else {
    return(as.character(value))
  }
}

###### CODE

# Sets the path of the directory containing the output of FADE
setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")

bs <- "uniform"
dirname <- sprintf("%s,32/ReportCodingSeq", bs)
df1Filename <- sprintf("%s/CodingSeqDistancesAll.RDS", dirname )
df2Filename <- sprintf("%s/CodingSeqCVAll.RDS", dirname )

similarities  <- c('D2')
# misure di riferimento
l1 <- c("D2")
# misure analizzate
l2 <- c("Antidice", "Dice", "Jaccard", "Kulczynski", "Ochiai", "Russel")
# misure dominate da A e D (B e C diventano irrilevanti)
l3 <- c("Hamman", "Hamming", "Matching", "Sneath", "Tanimoto")
# misure escluse dal calcolo
l4 <- c("Anderberg", "Gower", "Yule", "Mash.distance.1000.", "Mash.distance.10000.", "Mash.distance.100000.", "Euclidean")

# riordina la lista delle misure
sortedMeasures  <- c(l1, l2, l3)
#measureNanesDF <- data.frame( ref = l1, g1 = l2, alt = l3, no = l4, stringsAsFactors = FALSE)

zoomLevels <- c(95, 90, 80, 70, 60) # soglia sul valore di theta per avere 1, 2, 3, 4, 5 valori sull'asse delle x
zoom <- 5

# Defines the name of the file containing a copy of the dataframe created by this script
#  Yeast, CElegans, HomoSapiens, Schistosoma, Lemur, MacacaMulatta, PiceaAbies
# genomes <- c( "Yeast", "CElegans", "HomoSapiens", "Schistosoma", "Lemur", "MacacaMulatta", "PiceaAbies")
genomes <- c( "Yeast-codingSeq", "CElegans-codingSeq", "HomoSapiens-codingSeq")
sortedGenomes <- c( "Yeast-codingSeq", "CElegans-codingSeq", "HomoSapiens-codingSeq")

types <- c("1:3", "4:6", "7:9", "9:11", "4:11", "all")
elements <- c("1:3", "4:6", "7:9", "9:11", "4:11", "1:11")
cvDF <- data.frame( Genome = character(), Measure = character(), k = integer(),
                    type = character(), cv = double(), stringsAsFactors=FALSE)

tgtDF <- data.frame( Genome = character(), Measure = character(), Theta = integer(), k = integer(),
                     A = numeric(), B = numeric(), C = numeric(), D = numeric(), N = numeric(), density = numeric(),
                     distance=double(), stringsAsFactors=FALSE)

nObs  <- 88
nRowXObs  <- 13
dfSize  <- nObs * nRowXObs


if (!file.exists(df1Filename) || !file.exists(df2Filename) ) {
  # calcola i 2 dataframe
  cnt  <- 0
  for( sequenceName in genomes) {

    dfFilename <- sprintf( "%s,32/Synthetics-DatiEsperimento/Report%s.RDS", bs, sequenceName)
    csvFilename <- sprintf("%s,32/Synthetics-DatiEsperimento/%s.csv", bs, sequenceName)

    if (!dir.exists(dirname)) {
      dir.create(dirname)
    }
    tgtDir  <- sprintf("%s/%s", dirname, sequenceName)
    if (!dir.exists( tgtDir)) {
      dir.create(tgtDir)
    }

    if (!file.exists(dfFilename)) {
      # carica il CSV dell'esperimento
      columnClasses  <- c(
        #   sequenceA  sequenceB  start.time  real.time    Theta        k
        "character", "character", "numeric", "numeric", "integer", "integer",
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

      df <-read.csv( file = csvFilename, sep = ",", dec = ".", colClasses = columnClasses)
      saveRDS( df, file = dfFilename)
      cat(sprintf("Dataset %s %d rows saved.\n", dfFilename, nrow(df)))
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

    df$kf  <- factor(df$k)
    df$tf  <- factor(df$Theta)

    measures  <- colnames(df)[13:27]
    # extras = c("Mash.Distance.10000.", "D2", "Euclidean")
    extras  <- c("D2", "Euclidean")
    measures  <- append(measures, extras)

    # calcola la trasposta ... un rigo per ogni misura
    for(i in 1:nrow(df)) {
      r <- df[i,]
      for(m in measures) {
        MesName <- if (m == "Mash.Distance.10000.") "Mash" else m
        nr <- c( sequenceName, MesName, r[5:12], df[i, m])
        tgtDF[nrow(tgtDF)+1,] <- nr
      }
    }


    # escludiamo le misure: euclidean norm, anderberg, gowel , phi e yule. 864 -> 672 observations
    tgtDF <- filter(tgtDF, Measure != "Anderberg" & Measure != "Gower" & Measure != "Phi" & Measure != "Yule" &
      Measure != "Euclid_norm" &
      Measure != "Mash.Distance.1000." & Measure != "Mash.Distance.10000." & Measure != "Mash.Distance.100000.")


    # nessun filtro 864 -> 864
    # tgtDF <- filter( tgtDF, Measure != "Euclid_norm" & Measure != "Mash.Distance.1000." & Measure != "Mash.Distance.100000.")

    cat(sprintf("Filtered measures: Anderberg, Gower, Phi, Yule, Euclid_norm, Mash.Distance.1000, Mash.Distance.10000, Mash.Distance.100000. (%d rows).\n",nrow(df)))

    tgtDF$k  <- factor(tgtDF$k)
    # tgtDF$Theta  <- factor(tgtDF$Theta)

    kValues  <- levels(tgtDF$k)
    measures <- levels(factor(tgtDF$Measure))

    cat(sprintf("Data Frame %s filtered (%d observations).\n", sequenceName, nrow(tgtDF) - cnt * dfSize))
    # tgtDF$Theta  <- factor(tgtDF$Theta)

    df_total  <- data.frame()
    for( km in kValues) {
      for ( mes in measures) {
        df <- filter(tgtDF, Genome == sequenceName & k == km & Measure == mes)
        # non calcoliamo piu' il rapporto incrementale
        # newcol = 'ri'
        # # cat(sprintf("k = %s, mes = %s, #rows %d\n", km, mes, nrow(df)))
        # for( i in 1:nrow(df) - 1) {
        #   # df[ i, newcol] <- df[i, 'distance'] / df[i+1, 'distance']
        #   # utilizziamo il rapporto incrementale
        #   df[ i, newcol] <- (df[i+1, 'distance'] - df[i, 'distance']) / ( df[i+1, 'Theta'] - df[i, 'Theta'])
        # }
        # solo per i 4 genomi
        if (sequenceName %in% sortedGenomes) {
          b  <- c( 0, 0, 0, 0, 0)
          for( i in 0:3 ) {
            s  <- if (i<3) i * 3 + 1 else i * 3
            cv  <- sd(df$distance[s:(s+2)]) / mean(df$distance[s:(s+2)])
            nr <- list( sequenceName, mes, as.integer(km), types[i+1], cv)
            cvDF[nrow( cvDF)+1,] <- nr
          }
          cv  <- sd(df$distance[4:11]) / mean(df$distance[4:11])
          nr <- list( sequenceName, mes, as.integer(km), types[5], cv)
          cvDF[nrow( cvDF)+1,] <- nr
          nr <- list( sequenceName, mes, as.integer(km), types[6], sd(df$distance)/mean(df$distance))
          cvDF[nrow( cvDF)+1,] <- nr
        }
        df_total <- rbind(df_total, df)
      }
    }

    cnt  <- cnt + 1
    if (nrow(tgtDF) != cnt * dfSize || nrow(df_total) != dfSize)   {
      stop("errore nel calcolo di df_total per il calcolo delle distanze relative")
    }
  } # per ogni genoma
  # salva il df finale

  saveRDS( tgtDF, df1Filename)

  # ordina i genomi per lunghezza crescente
  cvDF$Genome <- factor(cvDF$Genome, levels = sortedGenomes)
  cvDF$Measure <- factor(cvDF$Measure, levels = sortedMeasures)
  cvDF$k <-factor(cvDF$k,levels(factor(cvDF$k)))
  cvDF$type <-factor(cvDF$type)
  # salva il datafrme con i coefficienti di variazione
  saveRDS( cvDF, df2Filename)
  for( mes in measures) {
    df  <- filter( cvDF, Measure == mes)
    cat(sprintf("%s -> min: %f, max: %f\n", mes, min(df$cv, na.rm = TRUE), max(df$cv,na.rm = TRUE)))
  }
}
# i due dataframe già esistono N.B. cancellare per ricolacolare i valori
tgtDF <-readRDS( file = df1Filename)

cvDF <-readRDS( file = df2Filename)
# ordina i genomi per lunghezza crescente
cvDF$Genome <- factor(cvDF$Genome, levels = sortedGenomes)
cvDF$Measure <- factor(cvDF$Measure)
cvDF$k <-factor(cvDF$k,levels(factor(cvDF$k)))
cvDF$type <-factor(cvDF$type)

tgtDF$AD <- (tgtDF$A+tgtDF$D) / tgtDF$N

for( sequenceName in genomes) {

  totPrinted <- 0

  #  grafico distanze per ciascuna misura (Theta sull'asse delle x)
  df  <- filter(tgtDF, Genome == sequenceName & Measure != "D2" & Measure != "Euclidean")
  df$Measure  <- factor( df$Measure, levels = sortedMeasures)

  sp1 <- ggplot(df, aes(x = Theta, y = distance, fill = k)) +
      geom_line(aes(color = k)) +
      geom_point() +
      #  geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
      facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) + #, scales = "free_y"
      theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                             axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                             panel.spacing=unit(0.1, "lines"),
                             legend.position = "none",
                             axis.text.y = element_blank(),
                             axis.title.y = element_blank())
      # scale_y_continuous(name = "Distance")
      #                 breaks=c(0, 0.5, 1),
      #                 labels=c("0", "0.5", "1"))
      # labs(y = "Distance") +
      # guides(colour = guide_legend(override.aes = list(size=1)))


  # dev.new(width = 6, height = 6)
  # print(sp1)
  outfname <- sprintf( "%s/%s/PanelDistances-%s.pdf", dirname, sequenceName, sequenceName)
  ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
  dev.off() # only 129kb in size
  totPrinted <- totPrinted + 1


  # secondo grafico di riferimento per Euclide
  # df  <- filter(tgtDF, Genome == sequenceName & Measure == "Euclidean")
  #
  # sp1 <- ggplot(df, aes(x = Theta, y = distance)) +
  #   # geom_bar( width = 0.7, position = "dodge", stat = "identity") +
  #   geom_line(aes(color = k)) +
  #   geom_point() +
  #   # geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
  #   facet_grid( rows = vars(k), cols = vars(Measure), scales = "free_y", labeller = labeller( k = label_both)) +
  #   theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
  #                         axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
  #                         panel.spacing=unit(0.1, "lines"),
  #                         legend.position = "none",
  #                         axis.title.y = element_blank(),
  #                         axis.text.y = element_blank(),
  #                         strip.text.y = element_blank()) +
  #   labs( x = " ")
  #
  # # labs(y = "Euclidean Distance")
  # # guides(colour = guide_legend(override.aes = list(size=1)))
  #
  # # dev.new(width = 6, height = 9)
  # # print(sp1)
  # outfname <- sprintf( "%s/%s/PanelEuclid-%s.pdf", dirname, sequenceName, sequenceName)
  # ggsave( outfname, device = pdf(), width = 0.9, height = 6, units = "in", dpi = 300)
  # dev.off() # only 129kb in size
  # totPrinted <- totPrinted + 1

  # Terzo grafico di riferimento per simmilarità D2
  df  <- filter(tgtDF, Genome == sequenceName & Measure == "D2")

  sp1 <- ggplot(df, aes(x = Theta, y = distance)) +
    # geom_bar( width = 0.7, position = "dodge", stat = "identity") +
    geom_line(aes(color = k)) +
    geom_point() +
    # geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
    facet_grid( rows = vars(k), cols = vars(Measure), scales = "free_y", labeller = labeller( k = label_both)) +
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                          panel.spacing=unit(0.1, "lines"),
                          legend.position = "none",
                          axis.text.y = element_blank(),
                          strip.text.y = element_blank()) +
    labs( x = " ", y = "Distances")
  # labs(y = "D2 Similarity")
  # guides(colour = guide_legend(override.aes = list(size=1)))

  # dev.new(width = 6, height = 9)
  # print(sp1)
  outfname <- sprintf( "%s/%s/PanelD2-%s.pdf", dirname, sequenceName, sequenceName)
  ggsave( outfname, device = pdf(), width = 1.1, height = 6, units = "in", dpi = 300)
  dev.off() # only 129kb in size
  totPrinted <- totPrinted + 1


  # Pannello della density x k *******************************************************************
  #  grafico distanze per ciascuna misura (Theta sull'asse delle x)
  df  <- filter(tgtDF, Genome == sequenceName & Measure != "D2" & Measure != "Euclidean")
  df$Measure  <- factor( df$Measure, levels = sortedMeasures)

  sp1 <- ggplot(df, aes(x = Theta, y = distance, fill = k)) +
    geom_line(aes(color = k)) +
    geom_point( size = 0.8) +
    #  geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
    facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) + #, scales = "free_y"
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                          panel.spacing=unit(0.1, "lines"),
                          legend.position = "none",
                          axis.text.y = element_blank(),
                          axis.title.y = element_blank())
  # scale_y_continuous(name = "Distance")
  #                 breaks=c(0, 0.5, 1),
  #                 labels=c("0", "0.5", "1"))
  # labs(y = "Distance") +
  # guides(colour = guide_legend(override.aes = list(size=1)))

  # dev.new(width = 6, height = 6)
  outfname <- sprintf( "%s/%s/PanelDistances-%s.pdf", dirname, sequenceName, sequenceName)
  ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
  dev.off() # only 129kb in size
  totPrinted <- totPrinted + 1


  # grafico delle densità A/N
  df  <- filter(tgtDF, Genome == sequenceName & Measure == "Jaccard")

  sp1 <- ggplot(data=df, aes(x=Theta, y=density, label=density)) +
    geom_line(aes(color = k)) +
    geom_point() +
    geom_text(aes(label = round(density, 2)), size = 5, nudge_y = 0.2, show.legend = FALSE) +
    facet_grid( rows = vars(k), labeller = labeller( k = label_both)) +
    scale_y_continuous(name = "A/N",
                       breaks=c(0, 0.5, 1),
                       labels=c("0", "0.5", "1")) +
    theme_light() + theme( panel.spacing=unit(0.2, "lines"),
                           legend.position = "none",
                           strip.text.x = element_text( size = 8, angle = 0),
                           axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                           axis.title.x = element_blank())
  #                         axis.text.x=element_blank())

  outfname <- sprintf( "%s/%s/PanelDensities-%s.pdf", dirname, sequenceName, sequenceName)
  ggsave( outfname, device = pdf(), width = 6, height = 9, units = "in", dpi = 300)
  dev.off()
  totPrinted <- totPrinted + 1


  # Zooming ******************************************************************

  #  grafico distanze per ciascuna misura (Theta sull'asse delle x)
  df  <- filter(tgtDF, Genome == sequenceName & Measure != "D2" & Measure != "Euclidean" & Theta >= zoomLevels[zoom])
  df$Measure  <- factor( df$Measure, levels = sortedMeasures)

  sp1 <- ggplot(df, aes(x = Theta, y = distance, fill = k)) +
    geom_line(aes(color = k)) +
    geom_point(size = 0.8) +
    #  geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
    facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) + #, scales = "free_y"
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                          panel.spacing=unit(0.1, "lines"),
                          legend.position = "none",
                          axis.text.y = element_blank(),
                          axis.title.y = element_blank()) +
    scale_y_continuous(name = "Distance",
                          breaks=c(0, 0.5, 1),
                          labels=c("0", "0.5", "1")) +
    labs(y = "CV")
  # guides(colour = guide_legend(override.aes = list(size=1)))


  # dev.new(width = 6, height = 6)
  # print(sp1)
  outfname <- sprintf( "%s/%s/PanelDistances-%s-zoom%d.pdf", dirname, sequenceName, sequenceName, zoom)
  ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
  dev.off() # only 129kb in size
  totPrinted <- totPrinted + 1


  # df  <- filter(tgtDF, Genome == sequenceName & Measure == "Euclidean" & Theta >= zoomLevels[zoom])
  #
  # sp1 <- ggplot(df, aes(x = Theta, y = distance)) +
  #   # geom_bar( width = 0.7, position = "dodge", stat = "identity") +
  #   geom_line(aes(color = k)) +
  #   geom_point(size = 0.8) +
  #   # geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
  #   facet_grid( rows = vars(k), cols = vars(Measure), scales = "free_y", labeller = labeller( k = label_both)) +
  #   theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
  #                         axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
  #                         panel.spacing=unit(0.1, "lines"),
  #                         legend.position = "none",
  #                         axis.title.y = element_blank(),
  #                         axis.text.y = element_blank(),
  #                         strip.text.y = element_blank()) +
  #   labs( x = " ")
  #
  # # labs(y = "Euclidean Distance")
  # # guides(colour = guide_legend(override.aes = list(size=1)))
  #
  # # dev.new(width = 6, height = 9)
  # # print(sp1)
  # outfname <- sprintf( "%s/%s/PanelEuclid-%s-zoom%d.pdf", dirname, sequenceName, sequenceName, zoom)
  # ggsave( outfname, device = pdf(), width = 0.9, height = 6, units = "in", dpi = 300)
  # dev.off() # only 129kb in size
  # totPrinted <- totPrinted + 1


  df  <- filter(tgtDF, Genome == sequenceName & Measure == "D2" & Theta >= zoomLevels[zoom])

  sp1 <- ggplot(df, aes(x = Theta, y = distance)) +
    geom_line(aes(color = k)) +
    geom_point(size = 0.8) +
    # geom_text(aes(label = round(distance, 5)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
    facet_grid( rows = vars(k), cols = vars(Measure), scales = "free_y", labeller = labeller( k = label_both)) +
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                          panel.spacing=unit(0.1, "lines"),
                          legend.position = "none",
                          axis.text.y = element_blank(),
                          strip.text.y = element_blank()) +
    labs( x = " ", y = "Distances")
  # labs(y = "D2 Similarity")
  # guides(colour = guide_legend(override.aes = list(size=1)))

  # dev.new(width = 6, height = 9)
  # print(sp1)
  outfname <- sprintf( "%s/%s/PanelD2-%s-zoom%d.pdf", dirname, sequenceName, sequenceName, zoom)
  ggsave( outfname, device = pdf(), width = 1.1, height = 6, units = "in", dpi = 300)
  dev.off() # only 129kb in size
  totPrinted <- totPrinted + 1

  cat(sprintf("%s synthetics Done. %d plot printed\n", sequenceName, totPrinted))
}

# per ultimo il grafico con il coefficiente di variazione
# cvDF = 728 observations = 4 genome x 13 misure x 8 k

df  <- filter(cvDF, Measure != "D2" & Measure != "Euclidean" & type == "all")
df$Measure = factor( df$Measure, levels = sortedMeasures)

sp1 <- ggplot(df, aes(x = Genome, y = cv, fill = k)) +
  geom_point(size = 0.8, aes(color = k)) +
  facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) +
  theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                        axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                        panel.spacing=unit(0.1, "lines"),
                        legend.position = "none",
                        axis.text.y = element_blank(),
                        axis.title.y = element_blank(),
                        axis.title.x = element_blank())
  # scale_y_continuous(name  <- "CV")
  #                    breaks=c(0, 0.5, 1),
  #                    labels=c("0", "0.5", "1")) +
  # labs(y = "CV 1:11")

outfname <- sprintf( "%s/PanelCV-all.pdf", dirname)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1



# df  <- filter(cvDF, Measure == "Euclidean" & type == "all")
#
# sp1 <- ggplot(df, aes(x = Genome, y = cv)) +
#   geom_point(size = 0.8, aes(color = k)) +
#   facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) +
#   theme_light() + theme( legend.position = "none",
#                          panel.spacing = unit(0.1, "lines"),
#                          axis.title.x = element_blank(),
#                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
#                          strip.text.x = element_text( size = 8, angle = 0),
#                          axis.title.y = element_blank(),
#                          axis.text.y = element_blank(),
#                          strip.text.y = element_blank())
#
# outfname <- sprintf( "%s/PanelCVEuclid.pdf", dirname)
# ggsave( outfname, device = pdf(), width = 0.9, height = 6, units = "in", dpi = 300)
# dev.off() # only 129kb in size
# totPrinted <- totPrinted + 1

df  <- filter(cvDF, Measure == "D2" & type == "all")

sp1 <- ggplot(df, aes(x = Genome, y = cv)) +
  geom_point(size = 0.8, aes(color = k)) +
  facet_grid( rows = vars(k), cols = vars(Measure), , labeller = labeller( k = label_both)) +
  theme_light() + theme(panel.spacing=unit(0.1, "lines"),
                        legend.position = "none",
                        axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                        axis.title.x = element_blank(),
                        strip.text.x = element_text( size = 8, angle = 0),
                        axis.text.y = element_blank(),
                        strip.text.y = element_blank()) +

  labs( y = "Variability Index")

outfname <- sprintf( "%s/PanelCVD2.pdf", dirname)
ggsave( outfname, device = pdf(), width = 1.1, height = 6, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1

# zooming ---------------------------------------------------------------------------

for( i in 1:5) {
  df  <- filter(cvDF, Measure != "D2" & Measure != "Euclidean" & type == types[i])
  df$Measure = factor( df$Measure, levels = sortedMeasures)

  sp1 <- ggplot(df, aes(x = Genome, y = cv, fill = k)) +
    geom_point(size = 0.8, aes(color = k)) +
    facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) +
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                          panel.spacing=unit(0.1, "lines"),
                          legend.position = "none") +
    #                       axis.text.y = element_blank(),
    #                       axis.title.y = element_blank()) +
    # scale_y_continuous(name = "CV")
    #                    breaks=c(0, 0.5, 1),
    #                    labels=c("0", "0.5", "1")) +
    labs( x = " ", y = sprintf("Variability Index (%s)", elements[i]))

  outfname <- sprintf( "%s/PanelCV-all-zoom%d.pdf", dirname,i)
  ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
  dev.off() # only 129kb in size
  totPrinted <- totPrinted + 1

  # df  <- filter(cvDF, Measure == "Euclidean" & type == types[i])
  #
  # sp1 <- ggplot(df, aes(x = Genome, y = cv)) +
  #   geom_point(size = 0.8, aes(color = k)) +
  #   facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) +
  #   theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
  #                         axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
  #                         panel.spacing=unit(0.1, "lines"),
  #                         legend.position = "none") +
  #   #                      axis.title.y = element_blank(),
  #   #                      axis.text.y = element_blank(),
  #   #                      strip.text.y = element_blank()) +
  #   labs( x = " ", y = sprintf("Variability Index (%s)", elements[i]))
  #
  # outfname <- sprintf( "%s/PanelCVEuclid-zoom%d.pdf", dirname, i)
  # ggsave( outfname, device = pdf(), width = 1.1, height = 6, units = "in", dpi = 300)
  # dev.off() # only 129kb in size
  # totPrinted <- totPrinted + 1

  df  <- filter(cvDF, Measure == "D2" & type == types[i])

  sp1 <- ggplot(df, aes(x = Genome, y = cv)) +
    geom_point(size = 0.8, aes(color = k)) +
    facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) +
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                          panel.spacing=unit(0.1, "lines"),
                          legend.position = "none") +
    #                      axis.text.y = element_blank(),
    #                      strip.text.y = element_blank()) +
    labs( x = " ", y = sprintf("Variability Index (%s)", elements[i]))

  outfname <- sprintf( "%s/PanelCVD2-zoom%d.pdf", dirname,i)
  ggsave( outfname, device = pdf(), width = 1.1, height = 6, units = "in", dpi = 300)
  dev.off() # only 129kb in size
  totPrinted <- totPrinted + 1
}

df  <- filter(cvDF, Measure != "D2" & Measure != "Euclidean" & type != "all" & type != "4:11")
df$Measure  <- factor( df$Measure, levels = sortedMeasures)

sp1 <- ggplot(df, aes(x = Genome, y = cv, fill = type)) +
  geom_point(size = 0.8, aes(color = type)) +
  facet_grid( rows = vars(k), cols = vars(Measure), labeller = labeller( k = label_both)) +
  theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                        axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                        panel.spacing=unit(0.1, "lines")) +
  #                      legend.position = "none") +
  #                       axis.text.y = element_blank(),
  #                       axis.title.y = element_blank()) +
  # scale_y_continuous(name = "CV")
  #                    breaks=c(0, 0.5, 1),
  #                    labels=c("0", "0.5", "1")) +
  labs( x = " ", y = "Variability Index")

outfname <- sprintf( "%s/PanelCV-all-zoomall.pdf", dirname)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off() # only 129kb in size
totPrinted <- totPrinted + 1


# grafico delle densità A/N x tutti i genomi
df  <- filter(tgtDF, Genome %in% sortedGenomes, Measure == "Jaccard")

df$Genome <- factor(df$Genome, levels = sortedGenomes)

sp1 <- ggplot(data=df, aes(x=Theta, y=density, label=density)) +
  geom_line(aes(color = k)) +
  geom_point() +
  geom_text(aes(label = round(density, 1)), size = 1.8, nudge_y = 0.2, show.legend = FALSE) +
  facet_grid( rows = vars(k), cols = vars( Genome), labeller = labeller( k = label_both)) +
  scale_y_continuous(name = "Genome A/N Values",
                     breaks=c(0, 0.5, 1),
                     labels=c("0", "0.5", "1")) +
  theme_light() + theme( panel.spacing=unit(0.2, "lines"),
                         legend.position = "none",
                         strip.text.x = element_text( size = 8, angle = 0),
                         axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                         axis.title.x = element_blank())
#                         axis.text.x=element_blank())

outfname <- sprintf( "%s/PanelDensities-all.pdf", dirname)
ggsave( outfname, device = pdf(), width = 6, height = 8, units = "in", dpi = 300)
dev.off()
totPrinted <- totPrinted + 1

# grafico delle densità (A+D)/N x tutti i genomi
df  <- filter(tgtDF, Genome %in% sortedGenomes, Measure == "Jaccard")

df$Genome <- factor(df$Genome, levels = sortedGenomes)

sp1 <- ggplot(data=df, aes(x=Theta, y=AD, label=density)) +
  geom_line(aes(color = k)) +
  geom_point() +
  geom_text(aes(label = round(AD, 1)), size = 1.8, nudge_y = 0.2, show.legend = FALSE) +
  facet_grid( rows = vars(k), cols = vars( Genome), labeller = labeller( k = label_both)) +
  scale_y_continuous(name = "Genome (A+D)/N Values",
                     breaks=c(0, 0.5, 1),
                     labels=c("0", "0.5", "1")) +
  theme_light() + theme( panel.spacing=unit(0.2, "lines"),
                         legend.position = "none",
                         strip.text.x = element_text( size = 8, angle = 0),
                         axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                         axis.title.x = element_blank())
#                         axis.text.x=element_blank())

outfname <- sprintf( "%s/PanelRari-all.pdf", dirname)
ggsave( outfname, device = pdf(), width = 6, height = 8, units = "in", dpi = 300)
dev.off()
totPrinted <- totPrinted + 1

cat(sprintf("CV plot Done. %d plot printed\n", totPrinted))
