library(DescTools)
library(dplyr)
library(ggplot2)
library(hrbrthemes)
library(r2r)


###### DESCRIPTION

# compute Power Statistics and T1 error from raw data produced by PresenAbsent3.py script
# in CSV format


###### OPTIONS
###### CODE

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
  if (variable == 'kv') {
    # N.B. kv e' un factor
    return(sprintf("k = %s", as.character(value)))
  } else if (variable == 'k') {
    return(sprintf("k = %d", value))
  }  else if (variable == 'lenFac') {
    # lenFac è un factor
    return(formatC(as.numeric(as.character(value)), format="f", digits=0, big.mark="."))
  }else {
    return(as.character(value))
  }
}

# Sets the path of the directory containing the output of FADE
setwd("~/Universita/Src/IdeaProjects/power_statistics/data/PresentAbsent")

bs <- "uniform"

similarities = c('D2')

# Defines the name of the file containing a copy of the dataframe created by this script
# dfFilename <- "PresentAbsent-Power+T1.RDS"
# csvFilename <- 'PresentAbsentData-all.csv'
# nullModel <- 'Uniform'

dfFilename <- sprintf( "%s,32/PresentAbsentEC-Power+T1-%s,32.RDS", bs, bs)
dirname <- sprintf("%s,32/PlotDensity", bs)
if (!dir.exists(dirname)) {
  dir.create(dirname)
}


nullModel <- 'Uniform'
# nullModel <- 'ShuffledEColi'
T1Model <- paste( sep='', nullModel, '-T1')

###### CODE

# 'data.frame':	15120 obs. of  13 variables:
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
  Measure != "Euclid_norm" & Measure != "Mash.Distance.1000." & Measure != "Mash.Distance.100000.")

cat(sprintf("Filtered measures: Anderberg, Gower, Phi, Yule, Euclid_norm, Mash.Distance.1000, Mash.Distance.100000. (%d rows).\n",nrow(df)))

df$Measure <- replace( df$Measure, df$Measure == "Mash.Distance.10000.", "Mash")

df$Model = factor(df$Model)
df$kv = factor(df$k)
# df$len = factor(df$len)
df$lenFac = factor(df$len)

kValues = levels(factor(df$k))
lengths = levels(factor(df$len))
measures <- levels(factor(df$Measure))
altModels = levels(df$Model)[1:2]
classes <- c("rare", "normal", "saturated")

thresholds = c(0.01, 0.99)



for(i in 1:nrow(df)) {
  t <- df[i, 'nmDensity']
  if (t < thresholds[1])        { lblnm <- classes[1]}
  else if (t > thresholds[2])   { lblnm <- classes[3]}
  else                         { lblnm <- classes[2]}
  df[i, 'classNM'] <- lblnm
  t <- df[i, 'amDensity']
  if (t < thresholds[1])        { lblam <- classes[1]}
  else if (t > thresholds[2])   { lblam <- classes[3]}
  else                         { lblam <- classes[2]}
  df[i, 'classAM'] <- lblam
}

df$classNM <- factor(df$classNM, levels = classes)
df$classAM <- factor(df$classAM, levels = classes)

# write.csv(nmdf, "density-k-len.csv", row.names=FALSE)
a = 0.05
#
# primo grafico nmDensity solo per Measure = Jaccard e AM = MotifReplace e alpha = 0.05
df2 = filter( df, alpha == a & Model == "MotifRepl-U"  & Measure == "Jaccard")

sp <- ggplot(data=df2, aes(x=len, y=nmDensity, label=nmDensity)) +
  geom_line(aes(color = kv)) +
  geom_point() +
  geom_text(aes(label = round(nmDensity, 3)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
  facet_grid( cols = vars(classNM), rows = vars(k), labeller = labeller( k = label_both)) +
  scale_x_continuous(name = NULL,
                     breaks=c(1000, 10000, 100000, 1000000, 10000000),
                     labels=c("1e+3", "1e+4", "1e+5", "1e+6", "1e+7"),
                     limits = c(1000, 10000000),
                     trans='log10') +
  scale_y_continuous(name = "A/N (Null Model)",
                     breaks=c(0, 0.5, 1),
                     labels=c("0", "0.5", "1")) +
  theme_light() + theme( panel.spacing=unit(0.2, "lines"),
                         legend.position = "none")

outfname <- sprintf( "%s/PanelnmDensity-A=%.2f.pdf", dirname, a)
ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
dev.off()
totPrinted <- 1


for(m in altModels) {
  for(g in c(0.01, 0.05, 0.10)) {
#
# secondo grafico amDensity solo per Measure = Jaccard e AM = MotifReplace e alpha = 0.05
    df2 = filter( df, alpha == 0.05 & gamma == g & Model == m  & Measure == "Jaccard")
    sp <- ggplot(data=df2, aes(x=len, y=amDensity, label=amDensity)) +
      geom_line(aes(color=kv)) +
      geom_point() +
      geom_text(aes(label = round(nmDensity, 3)), size = 3, nudge_y = 0.2, show.legend = FALSE) +
      facet_grid( cols = vars(classAM), rows = vars(k), labeller = labeller( k = label_both)) +
      scale_x_continuous(name = NULL,
                         breaks=c(1000, 10000, 100000, 1000000, 10000000),
                         labels=c("1e+3", "1e+4", "1e+5", "1e+6", "1e+7"),
                         limits = c(1000, 10000000),
                         trans='log10') +
      scale_y_continuous(name = sprintf("A/N (AM = %s, G=%.2f)", m, g),
                         breaks=c(0, 0.5, 1),
                         labels=c("0", "0.5", "1")) +
      theme_light() + theme( panel.spacing=unit(0.2, "lines"),
                             legend.position = "none")

    outfname <- sprintf( "%s/PanelamDensity-%s-G=%.2f.pdf", dirname, m, g)
    ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
    totPrinted <- totPrinted + 1
    dev.off()
  }
}
# T1 results
cat(sprintf("Dataset all classes -> %d rows.\n",  nrow(df)))

# plot scarsi divisi in 3 alpha altrimenti illegibile
for( a in c(0.01, 0.05, 0.10)) {
  df1 = filter( df, alpha == a)
  sp1 <- ggplot( df1, aes(x = Measure, y = T1, fill = kv)) +
    geom_bar( width = 0.7, position = "dodge", stat = "identity") +
    geom_hline(yintercept = a, linetype="dashed", color = "black") +
    facet_grid( cols =vars(classNM), rows = vars(lenFac), labeller = plot_labeller) +
    theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                          axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                          panel.spacing=unit(0.1, "lines")) +
    labs(y = sprintf("T1 results (A = %.2f)", a)) +
    # scale_x_continuous(trans='log10') +
    guides(colour = guide_legend(override.aes = list(size=1)))
  # ggtitle( am)

  # dev.new(width = 6, height = 6)
  # print(sp1)
  outfname <- sprintf( "%s/PanelT1-A=%.2f.pdf", dirname, a)
  ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
  dev.off() # only 129kb in size
  totPrinted <- totPrinted + 1
}


cat(sprintf("T1 Done.\n"))

#
# plot della power per i due AM e per i 3 valori di gamma
#
for(m in altModels) {
  for(a in c(0.01, 0.05, 0.10)) {
    for(g in c(0.01, 0.05, 0.10)) {

      df2 = filter( df, gamma == g & alpha == a & Model == m)

      sp2 <- ggplot( df2, aes(x = Measure, y = power, fill = kv)) +
        geom_bar( width = 0.7, position = "dodge", stat = "identity") +
        facet_grid( cols = vars(classAM), rows = vars(lenFac), labeller = plot_labeller) +
        theme_light() + theme(strip.text.x = element_text( size = 8, angle = 0),
                              axis.text.x = element_text( size = rel( 0.7), angle = 45, hjust=1),
                              panel.spacing=unit(0.1, "lines")) +
        labs(y = sprintf("Power results (A = %.2f, G = %.2f)", a, g)) +
        guides(colour = guide_legend(override.aes = list(size=1)))
      # ggtitle( am)

      # dev.new(width = 6, height = 6)
      # print(sp1)
      outfname <- sprintf( "%s/PanelPower-%s-A=%.2f-G=%.2f.pdf", dirname, m, a, g)
      ggsave( outfname, device = pdf(), width = 9, height = 6, units = "in", dpi = 300)
      dev.off() # only 129kb in size
      totPrinted <- totPrinted + 1
    }
  }
}
cat(sprintf("%d plot printed", totPrinted))