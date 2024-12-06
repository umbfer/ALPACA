library(rjson)
library(ggplot2)
library(dplyr)
library(RColorBrewer)



###### DESCRIPTION

# Produces two panels reporting the power trend, respectively for the PT and MR alternative models.
# In each panel, for each AF and gamma is reported the power level obtained across different
# values of n. It is also colored according to the value of k.
# The output is a set PNG images with name PanelPowerAnalysis-<AM>-A=0.10.png
# where <AM> reflects the alternative model being considered

# Note: this script must be executed after Power+T1-Json2RDS.R


###### OPTIONS

# Sets the path of the directory containing the input dataframe

setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset5-1000")

# Sets the name of the file containing the input dataframe

dfFilename <- "Power+T1-Results.RDS"

# Sets the output path for the images to be generated

dirname <- "PowerBoxPlot"


###### CODE

# rename in a human readable format the measure names
measure_names <- function( measure) {
  ris <- c()
  for( m in measure) {
    ris <- c(ris , str_to_title( switch( m,
                                         'chisquare' = 'chi square',
                                         'd2star' = 'd2*',
                                         'harmonicmean' = 'harmonic\nmean',
                                         'squaredchord' = 'squared\nchord',
                                         'jensenshannon' = 'jensen\nshannon',
                                         m)))
  }
  return( ris)
}

plot_labeller <- function(variable,value){
  # cat(sprintf("nome: %s, valore: %s\n", variable, as.character(value)))
  if (variable=='gamma') {
    # N.B. len e' un factor
    return(sprintf("G = %.2f", value))
  } else if (variable == 'measure') {
    # N.B. Measure e' un factor
    tr <- measure_names(as.character(value))
    # cat(sprintf("pre: %s\npost: %s\n", as.character(value), tr))
    return( tr)
    # return( as.character(value))
  } else {
    return(as.character(value))
  }
}

if (!dir.exists(dirname)) {
  dir.create(dirname)
}


if (!file.exists(dfFilename)) {
    cat( sprintf("Input Dataframe (%s) does not exist. Exiting\n", dfFilename))
	quit(save = "ask")
}

# carica il dataframe dal file
dati <- readRDS(file = dfFilename)
alphaTarget = 0.10

for (am in levels(dati$model)) {

	# solo per alpha = 0.10
    dff <- filter(dati, dati$alpha == alphaTarget & dati$model == am) # tutte le misure per uno specifico AM e valore di alpha

    sp <- ggplot( dff, aes( x = len, y = power, fill = k, alpha=0.8)) +
          geom_point( aes( color = k), alpha = 0.8, size = 0.08) +
          scale_x_continuous(name = NULL, breaks=c(200000, 2500000, 5000000, 7500000, 10000000),
                             labels=c("", "2.5e+6", "", "7.5e+6", ""), limits = c(200000, 10000000)) +
          scale_y_continuous(name = "Power", limits = c(0, 1)) +
          facet_grid( rows = vars( gamma), cols = vars( measure),  labeller = plot_labeller) +
          theme_bw() + theme(strip.text.x = element_text( size = 8, angle = 70),
                       axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                       panel.spacing=unit(0.1, "lines")) +
          guides(colour = guide_legend(override.aes = list(size=3)))
          # ggtitle( am)
    
	# dev.new(width = 9, height = 6)
	# print(sp)
	# stop("break")
	outfname <- sprintf( "%s/PanelPowerAnalysis-%s-A=%.2f.png", dirname, str_replace(am, " ", ""), alphaTarget)
	ggsave( outfname, device = png(), width = 9, height = 6, units = "in", dpi = 300)
	# dev.off() #only 129kb in size
}