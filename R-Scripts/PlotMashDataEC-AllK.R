library(DescTools)
library(ggplot2)
library(dplyr)
library(stringr)


###### DESCRIPTION

setwd("~/Universita/Progetti/BioInformatica/Present-Absent")

# Defines the name of the file containing a copy of the dataframe created by this script
csvFilename <- "EscherichiaColi-Mash.csv"

df <- read.csv( file = csvFilename)

# df$kf <- as.character(df$k)

# calcola una colonna con A / (N - D)
# df$rapporto <- df$A / (df$Nmax - df$D)

# for (i in 1:nrow(df)) {
#  if ((df$Nmax[i] - df$D[i]) == 0) {
#    cat(sprintf("N - D == 0 -> k = %d,  N = %.0f, D = %.0f\n", df$k[i], df$Nmax[i], df$D[i]))
#  }
  # if (df$rapporto[i] > 1) {
  #   # message(df[i,])
  #   cat(sprintf("j=%.5f, A = %.0f, D = %.0f, N = %.0f, k = %d\n", df$rapporto[i], df$A[i], df$D[i], df$Nmax[i], df$k[i]))
  #   df$rapporto[i] <- 1
  # }
#}

### merge with hamming distances results
# df2 <- readRDS("HammingDistanceEC.df")
#
# for( i in 1:nrow(df2)) {
#   #        model,        gamma,     k, JaccardIndex, A B C D Nmax, kf,    rapporto
#   l <- list(df2$Name[i], df2$Gamma[i], 0, df2$Dist[i], 0,0,0,0,0, 'hamming', 1-df2$Dist[i] )
#   df[nrow(df) + 1,] <- l
# }
#
# df$seqLen <- df2$len[1]
df$model <- as.factor(df$model)
# df$seqLen <- as.factor(df$seqLen)
df$kf <- as.factor(df$k)


dfNM <- filter(df, df$k <= 32)
sp <- ggplot( dfNM, aes(x = gamma, y = MashDistance, color = kf)) +
            geom_line() +
            geom_point(aes(shape = kf)) +
            facet_grid(rows = vars(sketchSize)) +
            scale_shape_manual(values = 1:17) +
            labs(y = "Mash Distance")


# dev.new(width = 9, height = 9)
outfname <- sprintf("MashEscherichiacoli-AllK.pdf")
ggsave( outfname, device = pdf(), width = 9, height = 9, dpi = 300)
# print(sp)
# dev.off() #only 129kb in size


