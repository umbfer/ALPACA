library("ggpubr")

###### DESCRIPTION

setwd("~/Universita/Progetti/BioInformatica/Present-Absent")

# Defines the name of the file containing a copy of the dataframe created by this script
csvFilename <- "Escherichia-us.csv"

df <- read.csv( file = csvFilename)

cat(sprintf("Pearson Correlation( A/S, A/(N-D)) = %f (sketchSize = %d)\n", cor( df$A.S1, df$A..N.D., method="pearson"), 1000))

cat(sprintf("Pearson Correlation( A/S, A/(N-D)) = %f (sketchSize = %d)\n", cor( df$A.S2, df$A..N.D., method="pearson"), 10000))

cat(sprintf("Pearson Correlation( A/S, A/(N-D)) = %f (sketchSize = %d)\n", cor( df$A.S3, df$A..N.D., method="pearson"), 100000))

