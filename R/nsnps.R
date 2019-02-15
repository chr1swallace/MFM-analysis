#!/usr/bin/env Rscript
## more signals from stoch vs stepwise

library(magrittr)
library(reshape)
library(GUESSFM)
library(ggplot2)
library(cowplot)

source("~/DIRS.txt")
comm <- paste("find",CVSROOT,"-mindepth 2 -maxdepth 2 |grep snpmod-99.RData")
comm <- comm[!grepl("004q-122973062-123565302",comm)]
snpmod.files <- system(comm,intern=TRUE)
source("~/Projects/cvs/runguess/R/functions.R")
source("~/Projects/cvs/runguess/R/plot2-functions.R")
#afiles <- list.files(CVSROOT,recursive=TRUE,pattern=paste0(".RData"),full=TRUE)
cols <- c(cols,Prior="#333333")

for(i in seq_along(snpmod.files)) {
    f <- snpmod.files[i]
    message(f)
    (load(f)) # SM2
    r <- basename(dirname(f)) 
    ## NSNP
    nsnp <- mypp.nsnp(SM2,expected=2)
    p <- plot(nsnp) + scale_colour_manual(values=cols) + background_grid()
    ggsave(p, file=file.path("~/scratch/MFM-output",r,"nsnp.pdf"))
}
