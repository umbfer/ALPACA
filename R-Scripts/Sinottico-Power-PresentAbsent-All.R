library(DescTools)
library(ggplot2)
# library(facetscales)
library(dplyr)
library(rjson)


###### DESCRIPTION

# Produces a panel made of 6 boxplots (3x2), organized in three rows (i.e., $\gamma=0.01,0.05,0.1$)
# and two columns (i.e., PT and MR alternate models). Each boxplot reports the
# distribution of the proportion  of true positives cumulatively by length, for each AF and k. 

# Note: this script must be executed after Power+T1-Json2RDS.R

###### OPTIONS


# Sets the path of the directory containing the output of FADE
setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")

# Sets the name of the file containing the input dataframe
dfFilename <- "PresentAbsent-Power+T1.RDS"


###### CODE

# rename in a human readable format the measure names
measure_names <- function( measure) {
  ris <- c()
  for( m in measure) {
    ris <- c(ris , str_to_title( switch( m,
                                         'Mash.Distance.1000.' = 'Mash (sz=10^3)',
                                         'Mash.Distance.10000.' = 'Mash (sz=10^4)',
                                         'Mash.Distance.100000.' = 'Mash (sz=10^5)',
                                         m)))
  }
  return( ris)
}

plot_labeller <- function(variable,value){
  # cat(sprintf("nome: %s, valore: %s\n", variable, as.character(value)))
  if (variable=='gamma') {
    # N.B. len e' un factor
    return(sprintf("G = %s", value))
  } else if (variable == 'Measure') {
    # N.B. Measure e' un factor
    tr <- measure_names(as.character(value))
    # cat(sprintf("pre: %s\npost: %s\n", as.character(value), tr))
    return( tr)
    # return( as.character(value))
  } else {
    return(as.character(value))
  }
}


# output dirname
dirname <- "PowerBoxPlot"
if (!dir.exists(dirname)) {
  dir.create(dirname)
}


if (!file.exists(dfFilename)) {
  cat( sprintf("Input Dataframe (%s) does not exist. Exiting\n", dfFilename))
  quit(save = "ask")
}

# carica il dataframe dal file
dati <- readRDS(file = dfFilename)
# solo per alpha = 0.10
alphaTarget = 0.10

dff <- filter(dati, dati$alpha == alphaTarget) # tutte le misure per uno specifico valore di alpha
dff$Measure2 <- dff$Measure
levels(dff$Measure2) <- measure_names(levels(dff$Measure))

sp <- ggplot(dff, aes(x = Measure2, y = power, fill = k, alpha = 0.7)) +
  geom_boxplot( aes( color = k), alpha = 0.7, outlier.size = 0.3) +
  facet_grid(rows = vars(gamma), cols = vars(Model)) +
  labs(y = 'Power') + theme_bw() +
  scale_x_continuous(trans='log10') +
  theme(axis.title.x = element_blank(),
        axis.text.x = element_text( size = rel(0.8), angle = 45, hjust = 1))

dev.new(width = 9, height = 6)
print(sp)
stop("break")
outfname <- sprintf( "%s/SinotticoPowerAllMeasures-AllK-AllLength-A=%.2f.png", dirname, alphaTarget)
ggsave( outfname, device = png(), width = 9, height = 6, units = "in", dpi = 300)
# dev.off() #only 129kb in size
