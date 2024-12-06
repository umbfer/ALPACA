library(DescTools)
library(ggplot2)
library(facetscales)
library(dplyr)


###### DESCRIPTION

# Produces a panel made of 6 boxplots (3x2), organized in three rows (i.e., $\gamma=0.01,0.05,0.1$)
# and two columns (i.e., PT and MR alternate models). Each boxplot reports the
# distribution of the proportion  of true positives cumulatively by length, for each AF and k. 

# Note: this script must be executed after Power+T1-Json2RDS.R

###### OPTIONS


# Sets the path of the directory containing the output of FADE
setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")

# Sets the name of the file containing the input dataframe
csvFilename <- 'GraficiEntropia-en.csv'
dfFilename <- "DatiEntropiaPresentAbsent.RDS"

###### CODE

if (!file.exists(dfFilename)) {
  columnClasses = c( "integer", "integer", "numeric", "numeric", "numeric", "numeric", "numeric")
  df <-read.csv( file = csvFilename, colClasses = columnClasses)
  # df$len = factor(df$len)
  df$k = factor(df$k)
  saveRDS( df, file = dfFilename)
} else {
  cat( sprintf("Data file %s exists. The csv file %s will NOT be loaded\n", dfFilename, csvFilename))
  df <- readRDS(file = dfFilename)
}


###### CODE
sp <- ggplot(df, aes(x=len, y=media)) +
  geom_line(aes(color = k, linetype = k)) +
  scale_x_continuous(name = 'seqLen', limits = c(1000, 10000000), trans='log10')

outfname <- sprintf('LineplotEntropyError-%s.pdf', 'average')
ggsave( outfname, device = pdf(), width = 9, height = 4, units = "in", dpi = 300)

sp <- ggplot(df, aes(x=len, y=varianza)) +
  geom_line(aes(color = k, linetype = k)) +
  scale_x_continuous(name = 'seqLen', limits = c(1000, 10000000), trans='log10')

outfname <- sprintf('LineplotEntropyError-%s.pdf', 'varianza')
ggsave( outfname, device = pdf(), width = 9, height = 4, units = "in", dpi = 300)

sp <- ggplot(df, aes(x = k, y = media)) +
  geom_boxplot( aes( color = k), alpha = 0.7, outlier.size = 0.3, width = 0.35)

outfname <- sprintf('BoxplotEntropyError-%s.pdf', 'average')
ggsave( outfname, device = pdf(), width = 9, height = 4, units = "in", dpi = 300)

sp <- ggplot(df, aes(x = k, y = varianza)) +
  geom_boxplot( aes( color = k), alpha = 0.7, outlier.size = 0.3, width = 0.35)

outfname <- sprintf('BoxplotEntropyError-%s.pdf', 'varianza')
ggsave( outfname, device = pdf(), width = 9, height = 4, units = "in", dpi = 300)

# dev.new(width = 9, height = 6)
# print(sp)
dev.off() #only 129kb in size
