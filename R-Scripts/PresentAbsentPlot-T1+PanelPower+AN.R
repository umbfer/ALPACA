library(rjson)
library(dplyr)
library(stringr)
library(RColorBrewer)
library(ggplot2)
library(r2r)
library(optparse)



###### DESCRIPTION

# Produces three different charts, reporting the results of the Type I check for
# each considered values of alpha (i.e., 0.01, 0.05, 0.1). 
# The output is a set PNG images with name T1Box-alpha=<x>.pdf
# where <x> reflects the actual value of alpha

# Note: this script must be executed after Power+T1-Json2RDS.R

###### OPTIONS

# Defines the name of the file containing a copy of the dataframe created by this script
# dfFilename <- "PresentAbsent-Power+T1.RDS"
# csvFilename <- 'PresentAbsentData-all.csv'
# nullModel <- 'Uniform'

#setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")


bs <- "uniform"
# Sets the name of the file containing the input dataframe

nullModel <- 'Uniform'
# nullModel <- 'ShuffledEColi'
T1Model <- paste( sep='', nullModel, '-T1')


option_list <- list(
  make_option(c("-d","--df"), type = "character",
              help = "Path to the output RDS file", metavar = "character"),
   make_option(c("--dirname"), type = "character",
              help = "Directory path for the output plots", metavar = "character")

)

opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (!is.null(opt$df)) {
    dfFilename <- opt$df
} else {
    dfFilename <- sprintf( "%s,32/PresentAbsentEC-Power+T1-%s,32.RDS", bs, bs)
    dfFilename <- 'uniform,32/PresentAbsentEC-Power+T1-uniform,32.RDS'
}

if (!is.null(opt$dirname)) {
  dirname <- opt$dirname
} else {
# Sets the output path for the images to be generated
dirname <- sprintf("%s,32/T1+Power-Plots", bs)
}

if (!dir.exists(dirname)) {
  dir.create(dirname)
}


###### CODE

if (!dir.exists(dirname)) {
	dir.create(dirname)
}


if (!file.exists(dfFilename)) {
  cat( sprintf("Input Dataframe (%s) does not exist. Exiting\n", dfFilename))
  quit(save = "ask")
}

# modifica i fattori di scala per ciascuna riga del pannello
TranslationTable  <- hashmap(default = 0)
TranslationTable[["Mash.Distance.10000."]] <- "Mash"


TerminologyServer <- function( key) {
  v = TranslationTable[[key]]
  if (v == 0) {
    return( key)
  } else {
    return( v)
  }
}


plot_labeller <- function(variable, value){
  # cat(sprintf("variable: <%s>, value: <%s>\n", variable, as.character(value)))
  if (variable == 'k') {
    return(sprintf("k=%s", as.character(value)))
  }  else if (variable == 'lenFac') {
    # lenFac è un factor
    return(formatC(as.numeric(as.character(value)), format="f", digits=0, big.mark="."))
  } else if (variable == "gamma") {
    return(sprintf("G=%.2f", as.double(value)))
  }else {
    return(as.character(value))
  }
}

scales_y <- list(
  `0.01` = scale_y_continuous(limits = c(0, 0.10), breaks = seq(0, 0.10, 0.02)),
  `0.05` = scale_y_continuous(limits = c(0, 0.20), breaks = seq(0, 0.20, 0.04)),
  `0.10` = scale_y_continuous(limits = c(0, 0.30), breaks = seq(0, 0.30, 0.06)))

# rename in a human readable format the measure names
measure_names <- function( measure) {
  ris <- c()
  for( m in measure) {
    if (as.character(m) != "D2") {
      ris <- c(ris , str_to_title( switch( m,
                                           'Mash.Distance.1000.' = 'Mash (sz=10^3)',
                                           'Mash.Distance.10000.' = 'Mash (sz=10^4)',
                                           'Mash.Distance.100000.' = 'Mash (sz=10^5)',
                                           m)))
    } else {
      ris <- c("D2", ris)
    }
  }
  return( ris)
}

# misure di riferimento
l1 <- c("D2", "Euclidean")
# misure analizzate
l2 <- c("Antidice", "Dice", "Jaccard", "Kulczynski", "Ochiai", "Russel")
# misure dominate da A e D (B e C diventano irrilevanti)
l3 <- c("Hamman", "Hamming", "Matching", "Sneath", "Tanimoto")
# misure escluse dal calcolo
l4 <- c("Anderberg", "Gower", "Yule", "Mash.distance.1000.", "Mash.distance.10000.", "Mash.distance.100000.")

# riordina la lista delle misure
sortedMeasures = c(l1, l2, l3)
#measureNanesDF <- data.frame( ref = l1, g1 = l2, alt = l3, no = l4, stringsAsFactors = FALSE)




# finally load input dataframe

###### CODE

# 'data.frame':	15120 obs. of  13 variables:
# 2 (Model) x 21 (Measure) x 5 (len) x 8 (k) x 3 (gamma) x 3 (alpha) = 15120
# $ Measure  : chr  "Anderberg" "Anderberg" "Anderberg" "Anderberg" ...
# $ Model    : Factor w/ 2 levels "MotifRepl-U",..: 1 2 1 2 1 2 1 2 1 2 ...
# $ len      : num  1000 1000 1000 1000 1000 1000 1000 1000 1000 1000 ...
# $ gamma    : num  0.01 0.01 0.01 0.01 0.01 0.01 0.05 0.05 0.05 0.05 ...
# $ k        : num  4 4 4 4 4 4 4 4 4 4 ...
# $ alpha    : num  0.01 0.01 0.05 0.05 0.1 0.1 0.01 0.01 0.05 0.05 ...
# $ threshold: num  0.337 0.337 0.411 0.411 0.461 ...
# $ power    : num  0 0.018 0.028 0.113 0.168 0.199 0 0.216 0.069 0.473 ...
# $ T1       : num  0.006 0.006 0.039 0.039 0.097 0.097 0.006 0.006 0.039 0.039 ...
# $ nmDensity: num  0.96 0.96 0.96 0.96 0.96 ...
# $ nmSD     : num  0.0129 0.0129 0.0129 0.0129 0.0129 ...
# $ amDensity: num  0.912 0.96 0.912 0.96 0.912 ...
# $ amSD     : num  0.0289 0.013 0.0289 0.013 0.0289 ...

df <-readRDS( file = dfFilename)
cat(sprintf("Dataset %s loaded. (%d rows).\n", dfFilename, nrow(df)))

# escludiamo le misure: euclidean norm, anderberg, gowel , phi e yule. 15120 -> 10800 observations
df <- filter( df, Measure != "Anderberg" & Measure != "Gower" & Measure != "Phi" & Measure != "Yule" &
  Measure != "Euclid_norm" & Measure != "Mash.Distance.1000." & Measure != "Mash.Distance.100000." &
  Measure != "Mash" )

cat(sprintf("Filtered measures: Anderberg, Gower, Phi, Yule, Euclid_norm, Mash.Distance.1000, Mash.Distance.10000., Mash.Distance.100000., Euclidean (%d rows).\n",nrow(df)))

# df$Measure <- replace( df$Measure, df$Measure == "Mash.Distance.10000.", "Mash")
df$Measure <- factor(df$Measure)
df$Model = factor(df$Model)
# df$len = factor(df$len)
df$lenFac = factor(df$len)

kValues = levels(factor(df$k))
lengths = levels(factor(df$len))
measures <- levels(factor(df$Measure))
altModels = levels(df$Model)[1:2]

#
# stampa 3 grafici per ciascun valore di alpha
#
AM = levels(df$Model)[1]
totPrinted <- 0
for( a in c( 0.01, 0.05, 0.10)) { # alpha values

  MaxT1 <- switch( sprintf("%.2f", a), "0.01" = 0.050, "0.05" = 0.150, "0.10" = 0.3) # fattore di amplificazione del valore di T1
  cat(sprintf("%.3f - %.3f\n", a, MaxT1))

  dff <- filter(df, df$alpha == a & df$gamma == 0.10 & df$Model == AM) # T1 Error Check does not depend on gamma and Alternate Model

  # riordina le misure
  # dff$measure2 <- factor(dff$Measure, levels = measure_names(levels(dff$Measure)))
  dff$Measure <- factor(dff$Measure, levels = sortedMeasures)
  dff$k <- factor(dff$k)

  # sp <- ggplot( dff, aes(x = Measure, y = T1)) +
  #   geom_boxplot( aes( fill = k), alpha=0.7, outlier.size = 0.25, lwd = 0) +
  #   # geom_boxplot( aes(fill = k), alpha=0.7, outlier.size = 0.25) +
  #   scale_y_continuous(name = "T1 Value", limits = c(0, MaxT1)) +
  #   geom_hline(yintercept = a, linetype="dashed", color = "black") +
  #   theme( axis.text.x = element_text(size = 11, angle = 45, hjust = 1),  # increase al font sizes
  #          axis.text.y = element_text(size = 12),
  #          legend.title = element_text(size = 14),
  #          legend.text = element_text(size = 13)) +
  #   labs( x = "") +
  #   # theme_light(base_size = 10) + labs(x = "") + # theme(legend.position = "none") +
  #   scale_colour_brewer(palette = "Dark2")
  # # scale_fill_grey(start = 0, end = .9)
  # # ggtitle("Pannello risultati T1-Check")
  sp <- ggplot( dff, aes( x = len, y = T1, alpha=0.8)) +
    geom_point( aes( color = k), alpha = 0.8, size = 1.1) +
    scale_x_log10(name = NULL, breaks=c(1000, 10000, 100000, 1000000, 10000000),
                  labels=c("10E3", "", "10E5", "", "10+E7"), limits = c(1000, 10000000)) +
    scale_y_continuous(name = "T1 Error", limits = c(0, MaxT1)) +
    geom_hline(yintercept = a, linetype="dashed", color = "black") +
    facet_grid( rows = vars( k), cols = vars( Measure),  labeller = plot_labeller) +
    theme_bw() + theme(strip.text.x = element_text( size = 8, angle = 70),
                       axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                       axis.text.y = element_text( size = rel( 0.7)),
                       legend.position = "none",
                       panel.spacing=unit(0.1, "lines")) +
    guides(colour = guide_legend(override.aes = list(size=3)))
  # ggtitle( am)

  # dev.new(width = 9, height = 6)
  # print(sp)
  # stop("break")
  # dev.new(width = 10, height = 5)
  outfname <- sprintf( "%s/T1Box-A=%.2f.pdf", dirname, a)
  ggsave( outfname, width = 9, height = 6, device = pdf(), dpi = 300)
  # print( sp)
  dev.off() #only 129kb in size
  totPrinted <- totPrinted + 1
}


# Plot 1 panel with power for each alternative model
l <- levels(factor(df$alpha))
alphaValues = unlist(lapply(l, as.numeric))
alphaTgt = alphaValues[2]

l <- levels(factor(df$gamma))
gammaValues = unlist(lapply(l, as.numeric))
gammaTgt = 0.05

for (gammaTgt in gammaValues) {

  for (alphaTgt in alphaValues) {
    cat(sprintf("Power for alpha = %.2f - gamma = %.2f\n", alphaTgt, gammaTgt))

    for (am in levels(df$Model)) {

      dff <- filter(df, df$alpha == alphaTgt & df$Model == am & df$gamma == gammaTgt) # tutte le misure per uno specifico AM, un valore di alpha ed un valore di gamma
      dff$k <- factor(dff$k)
      # riordina le misure
      # dff$Measure <- factor(dff$Measure, levels = measure_names(levels(dff$Measure)))
      dff$Measure <- factor(dff$Measure, levels = sortedMeasures)

      sp <- ggplot( dff, aes( x = len, y = power, alpha=0.8)) +
        geom_point( aes( color = k), alpha = 0.8, size = 1.1) +
        scale_x_log10(name = NULL, breaks=c(1000, 10000, 100000, 1000000, 10000000),
                           labels=c("10E3", "", "10E5", "", "10E7"), limits = c(1000, 10000000)) +
        scale_y_continuous(name = "Power", limits = c(0, 1), labels=c("0", "0.5", "1"), breaks = c(0, 0.5, 1)) +
        facet_grid( rows = vars( k), cols = vars( Measure),  labeller = plot_labeller) +
        theme_bw() + theme(strip.text.x = element_text( size = 8, angle = 70),
                           axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                           axis.text.y = element_text( size = rel( 0.7)),
                           legend.position = "none",
                           panel.spacing=unit(0.1, "lines")) +
        guides(colour = guide_legend(override.aes = list(size=3)))

      # ggtitle( am)

      # dev.new(width = 9, height = 6)
      # print(sp)
      # stop("break")
      outfname <- sprintf( "%s/PanelPowerAnalysis-%s-A=%.2f-G=%.2f.pdf", dirname, str_replace(am, " ", ""), alphaTgt,gammaTgt)
      ggsave( outfname, device = png(), width = 9, height = 4, units = "in", dpi = 300)
      dev.off() #only 129kb in size
      totPrinted <- totPrinted + 1
    }
  }
}

for (alphaTgt in alphaValues) {
  cat(sprintf("Power for alpha = %.2f - alla values of gamma\n", alphaTgt))
  for (am in levels(df$Model)) {

    dff <- filter(df, df$alpha == alphaTgt & df$Model == am) # tutte le misure per uno specifico AM, un valore di alpha ed un valore di gamma
    dff$k <- factor(dff$k)
    # riordina le misure
    # dff$Measure <- factor(dff$Measure, levels = measure_names(levels(dff$Measure)))
    dff$Measure <- factor(dff$Measure, levels = sortedMeasures)

    sp <- ggplot( dff, aes( x = len, y = power, alpha=0.8)) +
      geom_point( aes( color = k), alpha = 0.8, size = 1.1) +
      scale_x_log10(name = NULL, breaks=c(1000, 10000, 100000, 1000000, 10000000),
                    labels=c("10E3", "", "10E5", "", "10E7"), limits = c(1000, 10000000)) +
      scale_y_continuous(name = "Power", limits = c(0, 1), labels=c("0", "0.5", "1"), breaks = c(0, 0.5, 1)) +
      facet_grid( rows = vars( k, gamma), cols = vars( Measure),  labeller = plot_labeller) +
      theme_bw() + theme(strip.text.x = element_text( size = 8, angle = 70),
                         strip.text.y = element_text( size = 6),
                         axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                         axis.text.y = element_text( size = rel( 0.7)),
                         legend.position = "none",
                         panel.spacing=unit(0.1, "lines")) +
      guides(colour = guide_legend(override.aes = list(size=3)))
    # ggtitle( am)

    if (alphaTgt == 0.05) {
      if (as.character(am) == "MotifRepl-U") {
        sp <- sp + theme( axis.title.y = element_blank())
      } else {
        sp <- sp + theme(legend.position = "none")
      }
    }
    # dev.new(width = 9, height = 6)
    # print(sp)
    outfname <- sprintf( "%s/PanelPowerAnalysis-%s-A=%.2f.pdf", dirname, str_replace(am, " ", ""), alphaTgt)
    ggsave( outfname, device = png(), width = 6, height = 9, units = "in", dpi = 300)
    dev.off() #only 129kb in size
    totPrinted <- totPrinted + 1
  }
}

#
# Grafico del rapport A/N per il null model
#

# indipendente da alpha e da gamma fissati a piacere
# indipendente dall'alternative model (scegliamo MotifRepl-U)
# uguale per tutte le misure present absent (ne scegliamo una a caso == Jaccard)
dff <- filter(df, Measure == "Jaccard" & Model == "MotifRepl-U" & alpha == 0.05 & gamma == 0.05)
#
# ne restano 5 (len) x 8 (k) = 40 valori
#
if (nrow(dff) != 40) {
  stop("Errore nei dati input")
}
dff$k <- factor(dff$k)

sp <- ggplot( dff, aes( x = len, y = nmDensity, alpha=0.8)) +
  geom_point( aes( color = k), alpha = 0.8, size = 1.8) +
  scale_x_log10(name = NULL, breaks=c(1000, 10000, 100000, 1000000, 10000000),
                labels=c("10E3", "10E4", "10E5", "10E6", "10E7"), limits = c(1000, 10000000)) +
  scale_y_continuous(name = "Null Model A/N mean value", limits = c(0, 1), labels=c("0", "0.5", "1"), breaks = c(0, 0.5, 1)) +
  facet_grid( rows = vars( k), labeller = plot_labeller) +
  theme_bw() + theme(strip.text.x = element_text( size = 8, angle = 70),
                     strip.text.y = element_text( size = 10),
                     axis.text.x = element_text( size = rel( 1), angle = 45, hjust=1),
                     axis.text.y = element_text( size = rel( 1)),
                     legend.position = "none",
                     panel.spacing=unit(0.1, "lines")) +
  guides(colour = guide_legend(override.aes = list(size=3)))
# ggtitle( am)

# dev.new(width = 9, height = 6)
# print(sp)
outfname <- sprintf( "%s/Panel-NullModel-ANAnalysis.pdf", dirname)
ggsave( outfname, device = png(), width = 6, height = 9, units = "in", dpi = 300)
dev.off() #only 129kb in size
totPrinted <- totPrinted + 1


for (am in levels(factor(df$Model))) {
  # indipendente da alpha ma prendiamo tutti i gamma
  # uguale per tutte le misure present absent (ne scegliamo una a caso == Jaccard)
  dff <- filter(df, Measure == "Jaccard" & Model == am & alpha == 0.05 )
  #
  # ne restano 5 (len) x 8 (k) x 3 (gamma) = 120 valori
  #
  if (nrow(dff) != 120) {
    stop("Errore nei dati input")
  }
  dff$k <- factor(dff$k)
  yTitle = sprintf("%s Model A/N mean value", am)

  sp <- ggplot( dff, aes( x = len, y = amDensity, alpha=0.8)) +
    geom_point( aes( color = k), alpha = 0.8, size = 1.8) +
    scale_x_log10(name = NULL, breaks=c(1000, 10000, 100000, 1000000, 10000000),
                  labels=c("10E3", "10E4", "10E5", "10E6", "10E7"), limits = c(1000, 10000000)) +
    scale_y_continuous(name = yTitle, limits = c(0, 1), labels=c("0", "0.5", "1"), breaks = c(0, 0.5, 1)) +
    facet_grid( rows = vars( k), cols = vars( gamma), labeller = plot_labeller) +
    theme_bw() + theme(strip.text.x = element_text( size = 10, angle = 00),
                       strip.text.y = element_text( size = 10),
                       axis.text.x = element_text( size = rel( 1), angle = 45, hjust=1),
                       axis.text.y = element_text( size = rel( 1)),
                       legend.position = "none",
                       panel.spacing=unit(0.1, "lines")) +
    guides(colour = guide_legend(override.aes = list(size=3)))
  # ggtitle( am)

  # dev.new(width = 9, height = 6)
  # print(sp)
  outfname <- sprintf( "%s/Panel-%s-ANAvAnalysis.pdf", dirname, am)
  ggsave( outfname, device = png(), width = 6, height = 9, units = "in", dpi = 300)
  dev.off() #only 129kb in size
  totPrinted <- totPrinted + 1

  sp <- ggplot( dff, aes( x = len, y = nmSD, alpha=0.8)) +
    geom_point( aes( color = k), alpha = 0.8, size = 1.8) +
    scale_x_log10(name = NULL, breaks=c(1000, 10000, 100000, 1000000, 10000000),
                  labels=c("10E3", "10E4", "10E5", "10E6", "10E7"), limits = c(1000, 10000000)) +
    scale_y_continuous(name = "Standard Deviation Null Model A/N") +
    facet_grid( rows = vars( k), labeller = plot_labeller, scales = "free_y") +
    theme_bw() + theme(strip.text.x = element_text( size = 10, angle = 00),
                       strip.text.y = element_text( size = 10),
                       axis.text.x = element_text( size = rel( 1), angle = 45, hjust=1),
                       axis.text.y = element_text( size = rel( 1)),
                       legend.position = "none",
                       panel.spacing=unit(0.1, "lines")) +
    guides(colour = guide_legend(override.aes = list(size=3)))
  # ggtitle( am)

  # dev.new(width = 9, height = 6)
  # print(sp)
  outfname <- sprintf( "%s/Panel-NullModel-ANSDAnalysis.pdf", dirname)
  ggsave( outfname, device = png(), width = 6, height = 9, units = "in", dpi = 300)
  dev.off() #only 129kb in size

  totPrinted <- totPrinted + 1
}

cat(sprintf("%d plot printed", totPrinted))


