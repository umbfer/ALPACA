library(DescTools)
library(dplyr)

###### DESCRIPTION

# Converts the collection of CSV files returned by the script runJaccard.sh/cafe in a single dataframe

###### OPTIONS

# Sets the path of the directory containing the output of FADE
setwd("~/Universita/Src/IdeaProjects/power_statistics/data/results/dataset5-1000/jaccard")

# Defines the name of the file containing a copy of the dataframe created by this script
dfFilename <- "RawDistances-Jaccard.RDS"			


###### CODE

ReadResults <- function( model, kVal, lenVal, g)
{
    ris <- data.frame( Distance = double(), Measure = character(), Model = character(), k = numeric(), len = numeric(), stringsAsFactors=TRUE)
  
    gVal <- if (g != 0) sprintf(".G=%.3f", g) else ""
    mdl <- switch(model,
               "MotifRepl-U" = {paste0("MR", gVal)},
               "PatTransf-U" = { paste0("PT", gVal)},
               "Uniform" = { 'NM'},
               "Uniform-T1" = { 'T1'},
               stop("Wrong model value")
        ) 
    # name of the CSV file with distances for all measures
    #		model = 'Uniform' -> dist-k=4_Uniform-1000.8600000.csv
    #		model = 'PatTransf-U' -> dist-k=4_PatTransf-U-1000.200000.G=0.010.csv
    f1 <- sprintf("k=%d/dist-k=%d_%s-%d.%d%s.csv", kVal, kVal, model, nPairs, lenVal, gVal)
 
    if (!file.exists(f1)) {
      cat(sprintf("Input csv: %s not found ... skipping.", f1))
      return(NULL)
    }
    tmp <- read.csv( file = f1)
    for(mes in sortedMeasures) {
    
      Distance <- get( mes, tmp)
      Measure <- mes
      Model <- mdl
      k <- kVal
      len <- lenVal	
      # append new data to the result data.frame
      ris <- rbind( ris, data.frame( Distance, Measure, Model, k, len, stringsAsFactors=TRUE))
      cat( '.')
    }
    cat( sprintf("%s ok\n", f1))
    
    return( ris)
}



if (file.exists(dfFilename)) {
  cat( sprintf("Data file %s exists. Do you want to overwrite (Y/N) ?\n", dfFilename))
  res <- readline()
  if (res != "Y") {
    quit(save="ask")
  }
}

models <- c('Uniform', 'Uniform-T1', 'MotifRepl-U', 'PatTransf-U')

nPairs = 1000
lengths = seq(200000, 10000000, 200000) # c(1000000, 5000000, 6000000, 10000000) # 
kValues = c(4, 6, 8, 10)

# measures = c( 'canberra', 'chebyshev', 'chisquare', 'd2', 'd2s', 'd2star', 'd2z', 'euclidean', 'harmonicmean', 'intersection', 'jeffrey', 'jensenshannon', 'kulczynski2', 'manhattan', 'squaredchord') # no jaccard e mash
# misure ordinate per famiglia
sortedMeasures <- c('jaccard')

# dataframe risultato
dati <- data.frame( Distance = double(), Measure = character(), Model = character(), k = numeric(), len = numeric(), stringsAsFactors=TRUE)

for( len in lengths)	{ # 50 lunghezze

	for( k in kValues) {  # 4 valori di k
	  
		for( model in models) { # 2 AM
			
		  if (substr(model, 1, 4) == 'Unif') {
		    ret <- ReadResults( model, k, len, 0)
		    if (!is.null(ret)) {
  		    dati <- rbind(dati, ret)
		    }
		  }
		  else {
  			for( g in c(0.01, 0.05, 0.10)) {               # 3 valori di gamma
  			  ret <- ReadResults( model, k, len, g)
  			  if (!is.null(ret)) {
  			    dati <- rbind( dati, ret)
  			  }
  			} # for all gamma
		  }
		} # for all models
		cat( sprintf("k = %d.  done.\n", k))
	} # for all k
	cat( sprintf("length = %d.  done.\n\n", len))
} # for all length

# Rename by name: change "jensenshannon" to "jensen shannon"
# levels(x)[levels(x)=="jensenshannon"] <- "jensen shannon"
# dati$Measure <- factor(dati$Measure) # sono già factor
# levels(dati$Measure) = sortedMeasures # ha prodotto effetti disastrosi
# dati$Model <- factor(dati$Model) # sono già factor
# levels(dati$Model) = c('NM', "T1", "MR.G=0.010", "MR.G=0.050", "MR.G=0.100", "PT.G=0.010", "PT.G=0.050", "PT.G=0.100")
# dati$k <- factor( dati$k, kValues) # non occorre che siano factor
# dati$len <- factor( dati$len, lengths) # non occorre che siano factor
saveRDS( dati, file = dfFilename)

cat(sprintf("Dataset %s %d rows saved.", dfFilename, nrow(dati)))
