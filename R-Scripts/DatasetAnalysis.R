library(ggplot2)
library(dplyr)




setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")
bs <- c(1, 11, 16, 32)


# Sets the name of the file containing the input dataframe
dfFilename1 <- "PresentAbsent-RawData.RDS"
nullModel1 <- 'Uniform'
nmLabel1 <- nullModel1
T1Model1 <- paste( sep='', nullModel1, '-T1')

dfFilename2 <- sprintf("%d,32/PresentAbsentEC-RawData-%d,32.RDS", bs[2], bs[2])
nullModel2 <- 'ShuffledEColi'
nmLabel2 <- sprintf("%s-%d", nullModel2, bs[2])
T1Model2 <- paste( sep='', nullModel2, '-T1')

dfFilename3 <- sprintf("%d,32/PresentAbsentEC-RawData-%d,32.RDS", bs[3], bs[3])
nullModel3 <- 'ShuffledEColi'
nmLabel3 <- sprintf("%s-%d", nullModel3, bs[3])
T1Model3 <- paste( sep='', nullModel3, '-T1')

dfFilename4 <- sprintf("%d,32/PresentAbsentEC-RawData-%d,32.RDS", bs[4], bs[4])
nullModel4 <- 'ShuffledEColi'
nmLabel4 <- sprintf("%s-%d", nullModel4, bs[4])
T1Model4 <- paste( sep='', nullModel4, '-T1')

dati1 <- readRDS(file = dfFilename1)
dati2 <- readRDS(file = dfFilename2)
dati3 <- readRDS(file = dfFilename3)
dati4 <- readRDS(file = dfFilename4)

altModels1 = levels(dati1$model)[1:2]
amLabels1 <- altModels1
altModels2 = levels(dati2$model)[1:2]
amLabels2 <- c(sprintf("%s-%d", altModels2[1], bs[2]), sprintf("%s-b=%d", altModels2[2], bs[2]))
altModels3 = levels(dati3$model)[1:2]
amLabels3 <- c(sprintf("%s-%d", altModels3[1], bs[3]), sprintf("%s-b=%d", altModels3[2], bs[3]))
altModels4 = levels(dati4$model)[1:2]
amLabels4 <- c(sprintf("%s-%d", altModels4[1], bs[4]), sprintf("%s-b=%d", altModels4[2], bs[4]))

alphaValues <- c( 0.01, 0.05, 0.10)

for (kv in 1:8 * 4) {
  cat(sprintf("-------------------------------------------------\nk = %2d\n", kv))
  
  # solo per alpha = 0.10
  NM1 <- filter(dati1, dati1$k == kv & dati1$model == nullModel1) # tutte le misure per lo specifico NM del vecchio generative model

  NM2 <- filter(dati2, dati2$k == kv & dati2$model == nullModel2) # tutte le misure per lo specifico NM del nuovo generative model
  
  NM3 <- filter(dati3, dati3$k == kv & dati3$model == nullModel3) # tutte le misure per lo specifico NM del nuovo generative model

  NM4 <- filter(dati4, dati4$k == kv & dati4$model == nullModel4) # tutte le misure per lo specifico NM del nuovo generative model
  
  AM11 <- filter(dati1, dati1$k == kv & dati1$model == altModels1[1]) # tutte le misure per uno specifico AM1
  AM12 <- filter(dati1, dati1$k == kv & dati1$model == altModels1[2]) # tutte le misure per uno specifico AM2
  
  AM21 <- filter(dati2, dati2$k == kv & dati2$model == altModels2[1]) # tutte le misure per uno specifico AM1
  AM22 <- filter(dati2, dati2$k == kv & dati2$model == altModels2[2]) # tutte le misure per uno specifico AM2

  AM31 <- filter(dati3, dati3$k == kv & dati3$model == altModels3[1]) # tutte le misure per uno specifico AM1
  AM32 <- filter(dati3, dati3$k == kv & dati3$model == altModels3[2]) # tutte le misure per uno specifico AM2
  
  AM41 <- filter(dati4, dati4$k == kv & dati4$model == altModels4[1]) # tutte le misure per uno specifico AM1
  AM42 <- filter(dati4, dati4$k == kv & dati4$model == altModels4[2]) # tutte le misure 
  

  cat(sprintf("NullModel = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f\n",
              nmLabel1, kv, nrow(NM1), min(NM1$Jaccard), max(NM1$Jaccard), mean(NM1$Jaccard)))
  cat(sprintf("NullModel = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f\n",
              nmLabel2, kv, nrow(NM2), min(NM2$Jaccard), max(NM2$Jaccard), mean(NM2$Jaccard)))
  cat(sprintf("NullModel = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f\n",
              nmLabel3, kv, nrow(NM3), min(NM3$Jaccard), max(NM3$Jaccard), mean(NM3$Jaccard)))
  cat(sprintf("NullModel = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f\n\n",
              nmLabel4, kv, nrow(NM4), min(NM4$Jaccard), max(NM4$Jaccard), mean(NM4$Jaccard)))
  # Alternative model 1
  g = 0.01
  dff <- filter( AM11, AM11$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels1[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))

  dff <- filter( AM21, AM21$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels2[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM31, AM31$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels3[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))

  dff <- filter( AM41, AM41$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n\n",
              amLabels4[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  g = 0.05
  dff <- filter( AM11, AM11$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels1[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM21, AM21$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels2[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM31, AM31$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels3[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM41, AM41$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n\n",
              amLabels4[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  g = 0.1
  dff <- filter( AM11, AM11$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels1[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM21, AM21$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels2[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM31, AM31$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels3[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM41, AM41$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n\n",
              amLabels4[1], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  # Alternative model 2
  g = 0.01
  dff <- filter( AM12, AM12$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels1[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM22, AM22$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels2[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM32, AM32$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels3[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM42, AM42$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n\n",
              amLabels4[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  g = 0.05
  dff <- filter( AM12, AM12$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels1[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM22, AM22$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels2[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM32, AM32$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels3[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM42, AM42$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n\n",
              amLabels4[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  g = 0.1
  dff <- filter( AM12, AM12$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels1[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM22, AM22$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels2[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM32, AM32$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels3[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
  
  dff <- filter( AM42, AM42$gamma == g) # tutte le misure per uno specifico AM e valore di alpha
  cat(sprintf("AltModel  = %17.17s, k = %2d, rows = %d, min = %.3f, max = %.3f, mean = %.5f, gamma = %.2f\n",
              amLabels4[2], kv, nrow(dff), min(dff$Jaccard), max(dff$Jaccard), mean(dff$Jaccard), g))
}
