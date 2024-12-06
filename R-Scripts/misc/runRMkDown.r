library(rmarkdown)

Sys.setenv(PATH = paste(Sys.getenv("PATH"), "/usr/local/texlive/2018/bin/x86_64-darwin", sep=.Platform$path.sep)) 

render("RisultatiTestKolmogorov.Rmd")
