library(DescTools)
library(dplyr)
library(rjson)


###### DESCRIPTION

# Converts the collection of JSON files containing the values about power statistics and
# T1 Error checks in a single dataframe containing a number of records equal to:
#   (# of sequence lenghts)x(# of values of k)x(# of AF measures)x(# of alphas)x(# of gammas)
#
# In this example, the number of records is:
#   50 x 4 x 15 x 3 x 3 = 54.000


###### OPTIONS

# Sets the path of the directory containing the JSON files to convert
setwd("/Users/pipp8/Universita/Src/IdeaProjects/PowerStatistics/data/results/dataset5-1000")


# Defines the name of the file containing a copy of the dataframe created by this script
dfFilename <- "Power+T1-Results.RDS"

###### CODE

sortedMeasures <- c('chebyshev', 'euclidean', 'manhattan',
                    'chisquare',
                    'canberra',
                    'd2', 'd2s', 'd2star', 'd2z',
                    'intersection', 'kulczynski2',
                    'harmonicmean', 'squaredchord',
                    'jeffrey', 'jensenshannon')

seqLengths <- 50 # 200.000 - 10.000.0000 step 200.000 # numero di risultati per ciascun plot

AMs <- c("Motif Replace", "Pattern Transfer")
alphaValues = c("010", "050", "100")

gammaValues = c(0.01, 0.05, 0.10)
nGamma <- length(gammaValues)

kValues = c(4, 6, 8, 10)
nK <- length(kValues)

nCol <- 6

if (file.exists(dfFilename)) {
  cat( sprintf("Data file %s exists. Do you want to overwrite (Y/N) ?\n", dfFilename))
  res <- readline()
  if (res != "Y") {
    quit(save="ask")
  }
}
  
dati <- data.frame() # Output unico dataframe 3 alpha x 2 AM x 15 misure x 600 = 54000 righe

for(alpha in alphaValues) {
  
  for (am in AMs) {
    
    jsonDir <- sprintf( "json-%s", alpha)
    prefix <- strsplit(am, ' ')
    pat <- sprintf("*-%s*", prefix[[1]][1])
    files <- list.files(jsonDir, pat)
    print( sprintf("Processing %d files from: %s/%s/%s", length(files), getwd(), jsonDir, pat))
    
    md = if (am == 'Motif Replace') 'MR' else 'PT'
    
    for (file in files) {
      
      if (grepl( "jaccard|mash", file, ignore.case = TRUE)) {
        cat(sprintf("skipping file: %s\n", file))
        next
      }
      
      cat( sprintf("Processing: %s, Alpha: 0.%s ... ", file, alpha))
    
      # ogni file json contiene i valori di Power e di T1 per 50 misure x 3 gamma x 4 k = 600 righe  
      exp <- fromJSON(file = paste(jsonDir, file, sep="/"))
      values <-exp[['values']]
      
      df <- lapply( values, function( p) { data.frame(matrix(unlist(p), ncol=nCol, byrow=T))})
      df <- do.call(rbind, df)
      colnames(df)[1:nCol] <- c("len", "alpha", "k", "power", "T1", "gamma")
      
      measureName <- gsub( "[dD]istance", "", exp$header$distanceName)
      
      # aggiunge una colonna con il nome della misura in ogni riga
      df$measure = measureName
      # ed una colonna con l'Alternate Model esteso
      df$model = am
      # ed una colonna con l'Alternate Model x gamma (non alpha !!!)
      for(r in 1:nrow(df)) {
        df$mds[r] = sprintf("%s.G=%.3f", md, df$gamma[r])
      }
      dati <- rbind(dati, df)
      
      cat("done.\n")
    } # for each file
  } # foreach AM
} # foreach alpha
# dati$measure <- factor(dati$measure, levels = sortedMeasures)
dati$measure <- factor(dati$measure, levels = sortedMeasures)
dati$model <- factor(dati$model)
dati$mds <- factor(dati$mds)
dati$k <- factor( dati$k, levels = c( 4, 6, 8, 10))
dati$alpha <- factor( dati$alpha)
saveRDS( dati, file = dfFilename)

cat(sprintf("Dataset %s %d rows saved.", dfFilename, nrow(dati)))
