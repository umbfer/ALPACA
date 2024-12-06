
setwd("/Users/pipp8/Universita/Src/IdeaProjects/power_statistics/data/results/Kolmogorov/k=4")

#setwd("/Volumes/Catalina/PowerStatistics/study/k=6")

files <- list.files(".", "*.dist")

i <- 1
for(file in files) {
	# df <- read.table(header=FALSE, file)
	# dist <- as.vector(df[,1])
	cat( sprintf("%d -> %s\n", i ,file))
	i <- i + 1
}

writeLines("")

x <- 7
y <- 8
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))

x <- 7
y <- 1
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n", files[x], files[y], ris[1], ris[2]))

x <- 7
y <- 2
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n", files[x], files[y], ris[1], ris[2]))

x <- 7
y <- 3
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))

#
x <- 7
y <- 4
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n", files[x], files[y], ris[1], ris[2]))
x <- 7
y <- 5
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n", files[x], files[y], ris[1], ris[2]))
x <- 7
y <- 6
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))
#
x <- 1
y <- 4
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n", files[x], files[y], ris[1], ris[2]))
x <- 2
y <- 5
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n", files[x], files[y], ris[1], ris[2]))
x <- 3
y <- 6
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))

x <- 1
y <- 2
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n", files[x], files[y], ris[1], ris[2]))
x <- 1
y <- 3
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))

x <- 4
y <- 5
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n", files[x], files[y], ris[1], ris[2]))
x <- 4
y <- 6
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))

x <- 4
y <- 1
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))
x <- 5
y <- 2
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))
x <- 6
y <- 3
df <- read.table(header=FALSE, files[x])
xv <- df[,1]
df <- read.table(header=FALSE, files[y])
yv <- df[,1]
ris <- ks.test(xv, yv)
cat( sprintf("%s vs %s -> D:%f, P:%f\n\n", files[x], files[y], ris[1], ris[2]))

