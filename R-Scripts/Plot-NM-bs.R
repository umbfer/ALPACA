library(ggplot2)
library(dplyr)




setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")

gamma <- 0.10

dirPath <- sprintf("ModelAnalysis-G=%d", gamma * 100)
if (!dir.exists(dirPath)) {
  dir.create(dirPath)
}


# geneSize <- c(1, 3, 5, 7, 9, 11, 16, 32)
geneSize <- c(1, 8, 16, 32)

pltMeasures <- c("Jaccard", "D2", "Euclidean", "Euclid_norm")

alphaValues <- c( 0.01, 0.05, 0.10)

sortedLbl <- c()

df <- data.frame(  bs=numeric(), model = character(), seqLen = numeric(), k = numeric(), pairId = numeric(), 
                         distance = double(), stringsAsFactors=FALSE)

for (bs in geneSize) {
  cat(sprintf("-------------------------------------------------\nbs = %2d\n", bs))

  if (bs == 0) {
    # Sets the name of the file containing the input dataframe
    dfFilename <- "PresentAbsent-RawData.RDS"
    nullModel <- 'Uniform'
    nmLabel <- nullModel
    T1Model <- paste( sep='', nullModel, '-T1')
    dati <- readRDS(file = dfFilename)
    altModels = levels(dati$model)[1:2]
    amLabels <- altModels
    models <- c(nullModel, altModels)
    labels <- c(nmLabel, amLabels)
    sortedLbl <- c(sortedLbl, labels)
  } 
  else {
    dfFilename <- sprintf("%d,32/PresentAbsentEC-RawData-%d,32.RDS", bs, bs)
    nullModel <- 'ShuffledEColi'
    nmLabel <- sprintf("%s-%d", nullModel, bs)
    T1Model <- paste( sep='', nullModel, '-T1')
    dati <- readRDS(file = dfFilename)
    altModels = levels(dati$model)[1:2]
    amLabels <- c(sprintf("%s-%d", altModels[1], bs), sprintf("%s-%d", altModels[2], bs))
    models <- c(nullModel, altModels)
    labels <- c(nmLabel, amLabels)
    sortedLbl <- c(sortedLbl, labels)
  }

  for(i in 1:length(models))
    for (mes in pltMeasures) { # in models  
      # solo per alpha = 0.10
      NM <- filter(dati, dati$model == models[i] & (dati$gamma == 0 | dati$gamma == gamma) & dati$k != 4) # nullModel tutte le misure 
    
      tmpDF <- data.frame(seqLen = NM$seqLen)
      tmpDF$k <- NM$k
      tmpDF$pairId <- NM$pairId
      tmpDF$distance <- NM[[mes]]
      tmpDF$measure <- mes
      tmpDF$bs <- bs
      tmpDF$model <- labels[i]
      
      df <- rbind(df, tmpDF)
    }
}

# stop("break")
df$kf <- factor(df$k)
df$bsf <- factor(df$bs)
df$model <- factor(df$model, levels = sortedLbl)

print(levels(factor(df$model)))

for( mes in pltMeasures) {
  
  dfp <- filter(df, df$measure == mes) # nullModel tutte le misure 

  mesTitle <- sprintf("%s Distance", mes)
  
  sp <- ggplot( dfp, aes(x = model, y = distance)) +
    geom_boxplot( aes( color = kf, fill = kf), alpha=0.7, outlier.size = 0.25, lwd = 0) +
    facet_grid(rows = vars(seqLen), scales = "free_y") +
    # geom_boxplot( aes(fill = k), alpha=0.7, outlier.size = 0.25) +
    scale_y_continuous(name = mesTitle) + # , limits = c(0.70, 1)) +
    theme( axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # increase al font sizes
           axis.text.y = element_text(size = 8),
           legend.title = element_text(size = 14),
           legend.text = element_text(size = 12)) +
    labs( x = "Generative Model")
    # theme_light(base_size = 10) + labs(x = "") + # theme(legend.position = "none") +
    # scale_colour_brewer(palette = "Dark2")
    # scale_fill_grey(start = 0, end = .9)
    ggtitle("Null Model Analysis")
  
  # dev.new(width = 10, height = 10)
  outfname <- sprintf( "%s/ModelAnalysisByLen-%s-G=%d.pdf", dirPath, mes, gamma*100)
  ggsave( outfname, width = 9, height = 9, device = pdf(), dpi = 300)
  # print( sp)
  # dev.off() #only 129kb in size
  
  sp2 <- ggplot( dfp, aes(x = model, y = distance)) +
    geom_boxplot( aes( color = kf, fill = kf), alpha=0.7, outlier.size = 0.25, lwd = 0) +
    # geom_boxplot( aes(fill = k), alpha=0.7, outlier.size = 0.25) +
    scale_y_continuous(name = mesTitle) + # , limits = c(0.70, 1)) +
    theme( axis.text.x = element_text(size = 8, angle = 45, hjust = 1),  # increase al font sizes
           axis.text.y = element_text(size = 8),
           legend.title = element_text(size = 14),
           legend.text = element_text(size = 12)) +
    labs( x = "Generative Model")
  # theme_light(base_size = 10) + labs(x = "") + # theme(legend.position = "none") +
  # scale_colour_brewer(palette = "Dark2")
  # scale_fill_grey(start = 0, end = .9)
  ggtitle("Null Model Analysis By Lenghts")
  
  # dev.new(width = 10, height = 5)
  outfname <- sprintf( "%s/ModelAnalysis-%s-G=%d.pdf", dirPath, mes, gamma*100)
  ggsave( outfname, width = 9, height = 4, device = pdf(), dpi = 300)
  # print( sp2)
  # dev.off() #only 129kb in size
}