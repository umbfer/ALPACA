library(DescTools)
library(ggplot2)
library(dplyr)
library(stringr)


###### DESCRIPTION

setwd("~/Universita/Progetti/BioInformatica/Present-Absent")

# Defines the name of the file containing a copy of the dataframe created by this script
csvFilename <- "DatiJaccard-us.csv"

df <- read.csv( file = csvFilename)

###### CODE
df$seqLen <- as.factor(df$seqLen)
df$kf <- as.factor(df$k)
df$model <- as.factor(df$model)

# calcola una colonna con A / (N - D)
df$rapporto <- df$A/ (df$Nmax - df$D)

for (i in 1:nrow(df)) {
  v <- if ((df$Nmax[i] - df$D[i]) == 0) 1.0 else  (df$A[i] / (df$Nmax[i] - df$D[i]))
  df$rapporto[i] <- v
  if (v > 1) {
    # message(df[i,])
    cat(sprintf("v=%.3f, A = %.3f, D = %.3f, N = %.3f, k = %d\n", v, df$A[i], df$D[i], df$Nmax[i], df$k[i]))
  }
}


#dfNM <- filter(df, df$Model == 'MotifRepl')
sp <- ggplot( df, aes(x = seqLen, y = rapporto, fill = kf)) + geom_boxplot() +
    # geom_boxplot( aes(color = k), outlier.size = 0.3) +
    # scale_y_continuous(sec.axis = sec_axis(~ . * 10))
    facet_grid(rows = vars( model)) +
    # facet_grid_sc(cols = vars( len), rows = vars( k), scales = list( y = scales_y)) +
    theme( axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + # axis.text.y = element_blank()) +
    labs(x = "Sequence length") + labs(y = "A/(N-D)")
    # ggtitle(sprintf("Distances for k = %d", kv))
    

# dev.new(width = 16, height = 9)
outfname <- sprintf("Jaccard-AllK.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)
# print(sp)
# dev.off() #only 129kb in size


dfNM <- filter(df, df$model == 'MotifRepl')
sp <- ggplot( df, aes(x = seqLen, y = rapporto, fill = kf)) + geom_boxplot() +
  # geom_boxplot( aes(color = k), outlier.size = 0.3) +
  # scale_y_continuous(sec.axis = sec_axis(~ . * 10))
  facet_grid(rows = vars( gamma)) +
  theme( axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + # axis.text.y = element_blank()) +
  labs(x = "Sequence length (MotifReplace)") + labs(y = "A/(N-D)")
# ggtitle(sprintf("Distances for k = %d", kv))


# dev.new(width = 16, height = 9)
outfname <- sprintf("MotifReplace-Jaccard-AllK.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)

dfNM <- filter(df, df$model == 'PatTransf')
sp <- ggplot( df, aes(x = seqLen, y = rapporto, fill = kf)) + geom_boxplot() +
  # geom_boxplot( aes(color = k), outlier.size = 0.3) +
  # scale_y_continuous(sec.axis = sec_axis(~ . * 10))
  facet_grid(rows = vars( gamma)) +
  theme( axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + # axis.text.y = element_blank()) +
  labs(x = "Sequence length (PatternTransfer)") + labs(y = "A/(N-D)")
# ggtitle(sprintf("Distances for k = %d", kv))


# dev.new(width = 16, height = 9)
outfname <- sprintf("Pattern Transfer-Jaccard-AllK.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)