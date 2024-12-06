library(rjson)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(facetscales)



###### DESCRIPTION

# Produces three different charts, reporting the results of the Type I check for
# each considered values of alpha (i.e., 0.01, 0.05, 0.1). 
# The output is a set PNG images with name T1Box-alpha=<x>.png
# where <x> reflects the actual value of alpha

# Note: this script must be executed after Power+T1-Json2RDS.R

###### OPTIONS

# Sets the path of the directory containing the input dataframe

setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset5-1000")
	

# Sets the output path for the images to be generated

dirname <- "T1BoxPlot"

# Sets the name of the file containing the input dataframe

dfFilename <- "Power+T1-Results.RDS"

###### CODE

if (!dir.exists(dirname)) {
	dir.create(dirname)
}


if (!file.exists(dfFilename)) {
  cat( sprintf("Input Dataframe (%s) does not exist. Exiting\n", dfFilename))
  quit(save = "ask")
}

# carica il dataframe dal file
dfAll <- readRDS(file = dfFilename)

# modifica i fattori di scala per ciascuna riga del pannello
scales_y <- list(
  	`0.01` = scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02)),
  	`0.05` = scale_y_continuous(limits = c(0, 0.20), breaks = seq(0, 0.20, 0.04)),
  	`0.10` = scale_y_continuous(limits = c(0, 0.30), breaks = seq(0, 0.30, 0.06)))


cat(sprintf("Data Frame loaded. %d rows\n", nrow(dfAll)))

# rename in a human readable format the measure names
measure_names <- function( measure) {
  ris <- c()
  for( m in measure) {
    ris <- c(ris , str_to_title( switch( m,
                                         'chisquare' = 'chi square',
                                         'd2star' = 'd2*',
                                         'harmonicmean' = 'harmonic mean',
                                         'squaredchord' = 'squared chord',
                                         'jensenshannon' = 'jensen shannon',
                                         m)))
  }
  return( ris)
}


for( a in c( 0.01, 0.05, 0.10)) { 
  
  MaxT1 <- switch( sprintf("%.2f", a), "0.01" = 0.050, "0.05" = 0.150, "0.10" = 0.3) # fattore di amplificazione del valore di T1
  cat(sprintf("%.3f - %.3f\n", a, MaxT1))

  dff <- filter(dfAll, dfAll$alpha == a & dfAll$gamma == 0.10 & as.character(dfAll$model) == 'Motif Replace') # T1 Error Check does not depend on gamma and Alternate Model

  dff$measure2 <- dff$measure
  levels(dff$measure2) <- measure_names(levels(dff$measure))
  
	sp <- ggplot( dff, aes(x = measure2, y = T1)) + 
	  # geom_boxplot( aes(color = k, fill = k), alpha=0.7, outlier.size = 0.25) +
	 	geom_boxplot( aes(fill = k), alpha=0.7, outlier.size = 0.25) +
		scale_y_continuous(name = "T1 Value", limits = c(0, MaxT1)) +
	  geom_hline(yintercept = a, linetype="dashed", color = "black") +
	 	theme_bw() + theme(  axis.text.x = element_text(size = 11, angle = 45, hjust = 1),  # increase al font sizes
							 axis.text.y = element_text(size = 12),
							 legend.title = element_text(size = 14),
							 legend.text = element_text(size = 13)) +
	 	labs(x = "") + # theme(legend.position = "none") +
	  scale_fill_grey(start = 0, end = .9)
	 	# ggtitle("Pannello risultati T1-Check") 
	
	# dev.new(width = 10, height = 5)
    outfname <- sprintf( "%s/T1Box-alpha=%.2f.png", dirname, a)
    ggsave( outfname, width = 9, height = 4, device = png(), dpi = 300) 
    # print( sp)
    # dev.off() #only 129kb in size
}