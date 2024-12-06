library(ggplot2)
library(facetscales)
library(dplyr)


setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov")

len_names <- list(
  '200000' = "n = 200,000",
  '2000000' = "n = 2,000,000")

k_names <- list(
  '4' = "k = 4",
  '6' = "k = 6",
  '8' = "k = 8",
  '10' = "k = 10")

plot_labeller <- function(variable,value){
  if (variable=='len') {
    return(len_names[value])
  } else if (variable=='K') {
    return(k_names[value])
  } else {
    return(as.character(value))
  }
}

# modifica i fattori di scala per ciascuna riga del pannello
# N.B. l'etichetta del pannello deve essere alfanumerica non numerica
scales_y <- list(
    '4' = scale_y_continuous(limits = c(0, 0.15)),
    '6' = scale_y_continuous(limits = c(0, 0.05)),
    '8' = scale_y_continuous(limits = c(0, 0.01)),  
    '10' = scale_y_continuous(limits = c(0, 0.003)))


dfFilename <- 'AllKolmogorov.df'
ris <- readRDS(file = dfFilename)

ris <- filter( ris, ris$Name != 'C-NM')

# redefinisce l'ordine delle colonne della griglia
# levels(ris$len) <- c("n = 200,000","n = 2,000,000")

ris$len <- factor(ris$len)
ris$K <- factor(ris$K)

# levels(ris$K) <- c('k = 4', 'k = 6', 'k = 8', 'k = 10')

ris$Name <- factor( ris$Name, levels = c('NM', 'MR.G=0.010', 'MR.G=0.050', 'MR.G=0.100', 'PT.G=0.010', 'PT.G=0.050', 'PT.G=0.100'))


sp <- ggplot( ris, aes(x = Name,y = D, fill = Name, alpha=0.7)) + 
 	geom_boxplot( aes(color = Name), outlier.size = 0.3) +
 	facet_grid_sc(rows = vars( K), cols = vars( len), scales = list( y = scales_y), labeller = plot_labeller) +
 	theme_bw()+ theme( axis.text.x = element_text(size = 8, angle = 45, hjust = 1)) +
 	theme(legend.position = "none") + labs(y = "KS", x ="")
 	# ggtitle("Pannello risultati test di Kolmogorv-Smirnov") + labs(y= "D Value")


print( sp)
outfname <- 'AllKolmogorov.png'
ggsave( outfname, device = png(), dpi = 300)