library(DescTools)
library(ggplot2)
library(dplyr)
library(stringr)


###### DESCRIPTION

setwd("~/Universita/Progetti/BioInformatica/Present-Absent")

# Defines the name of the file containing a copy of the dataframe created by this script
csvFilename <- "DatiErrorRate-us.csv"

df <- read.csv( file = csvFilename)

###### CODE
df$len <- as.factor(df$len)
df$k <- as.factor(df$k)
df$Model <- as.factor(df$Model)

#dfNM <- filter(df, df$Model == 'MotifRepl')
sp <- ggplot( df, aes(x = len, y = error, fill = k)) + geom_boxplot() +
    # geom_boxplot( aes(color = k), outlier.size = 0.3) +
    # scale_y_continuous(sec.axis = sec_axis(~ . * 10))
    facet_grid(rows = vars( Model)) +
    # facet_grid_sc(cols = vars( len), rows = vars( k), scales = list( y = scales_y)) +
    theme( axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + # axis.text.y = element_blank()) +
    labs(x = "Sequence length") + labs(y = "delta/Hk")
    # ggtitle(sprintf("Distances for k = %d", kv))
    

# dev.new(width = 16, height = 9)
outfname <- sprintf("Errore(k)-AllK.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)
# print(sp)
# dev.off() #only 129kb in size


dfNM <- filter(df, df$Model == 'MotifRepl')
sp <- ggplot( df, aes(x = len, y = error, fill = k)) + geom_boxplot() +
  # geom_boxplot( aes(color = k), outlier.size = 0.3) +
  # scale_y_continuous(sec.axis = sec_axis(~ . * 10))
  facet_grid(rows = vars( G)) +
  theme( axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + # axis.text.y = element_blank()) +
  labs(x = "Sequence length (MotifReplace)", y = "delta/Hk")
# ggtitle(sprintf("Distances for k = %d", kv))


# dev.new(width = 16, height = 9)
outfname <- sprintf("MotifReplace-Errore-AllK.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)

dfNM <- filter(df, df$Model == 'PatTransf')
sp <- ggplot( df, aes(x = len, y = error, fill = k)) + geom_boxplot() +
  # geom_boxplot( aes(color = k), outlier.size = 0.3) +
  # scale_y_continuous(sec.axis = sec_axis(~ . * 10))
  facet_grid(rows = vars( G)) +
  theme( axis.text.x = element_text(size = 10, angle = 45, hjust = 1)) + # axis.text.y = element_blank()) +
  labs(x = "Sequence length (PatternTransfer)") + labs(y = "delta/Hk")
# ggtitle(sprintf("Distances for k = %d", kv))


# dev.new(width = 16, height = 9)
outfname <- sprintf("Pattern Transfer-Errore-AllK.pdf")
ggsave( outfname, device = pdf(), width = 6, height = 9, dpi = 300)