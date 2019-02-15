#!/usr/bin/env Rscript
## follows on from GUESS-qc.R, may need some code from that inserted

library(randomFunctions)
library(annotSnpStats)
library(magrittr)
source("~/DIRS.txt")
setwd("~/Projects/cd4chic/GUESSFM")

                                        #library(GUESSFM)
library(devtools)
load_all("~/RP/GUESSFM")

## files
## 1p-2406887-2785671 y
## 15q-67414055-67469568 y
## 10q-59823894-60140709 y
## 13q-42834735-43101763 x
## 20p-1497197-1689461 x
## 21q-16698439-16843298 --
## 20p-1497197-1689461 y
## 21q-43810084-43887145 y
## 22q-21712414-22022005 y
## 2    204733868    204737096

## 2q-191873553-192007734
## calculating marginal SNP inclusion probabilities
## Error in .local(x, i, j, ..., drop) : 
##   No match for one or more column selections
## Calls: ld -> [ -> [ -> .local
## Execution halted

args <- getArgs(default=list(d=file.path(CVSROOT,
                                         "2p-60914729-61892409"
                               ## "1q-192462312-192548041" #"10p-6030000-6220000" #         "10p-6030243-6169685"
                                        #14q-101290463-101328739"
                                         ## "1q-172650685-172940450"
                                         ),keep=FALSE))
## args <- getArgs(default=list(d=file.path(CVSROOT,"3p-45929800-46650993"),keep=FALSE))
args$d <- sub("/GUESS$","",args$d)
args

source("R/GUESS-common.R")

library(ggplot2)
theme_set(theme_bw())
library(gridExtra)

## dirs <- list.files(file.path(args$d,"GUESS"),full=TRUE)
## files <- unlist(lapply(dirs,list.files,full=TRUE,pattern="_features"))
##                                         # files <- list.files(file.path(args$d,args$ph),full=TRUE
## files <- grep("zzz|US|UK",files,invert=TRUE,value=TRUE)
## if(!length(files)) {
##     message("No output files found")
##     if(!interactive())
##         q("no")
## }
## files <- sub("_sweeps_features.txt","",files)
##                                         #files <- grep("_50000$",files,value=TRUE)
## names(files) <- basename(dirname(files))
## message("number of files found: ",length(files))
## cat(files,sep="\n")

## manhattans
file.conditional <- file.path(args$d,"conditional.RData")

condtest <- function(...) {
    best <- NULL
    while(TRUE) {
         newbest <- condtest.inner(...,best=best)
        ## (  newbest <- condtest.inner(X=snp.data, df=df, best=best))
        if(is.null(newbest))
            break
        best <- c(best,newbest)
    }
    best
}
condtest.inner <- function(X,df,best=NULL,stepwise.p.thr=1e-3,stepwise.max.predictors=NA, ...) {
 ## cs <- col.summary(X)
    ## X <- X[,cs[,"MAF"]>0.05]
    f <- "cc ~ PC1 + PC2 + PC3 + PC4"
    if(length(unique(df$country))>1)
        f <- paste(f,"+ country")
    if(is.null(best)) {
        Xtest <- X
        cond <- snp.rhs.tests(as.formula(f), snp.data=Xtest, data=df)
        p <- p.value(cond)
    } else {
        best <- names(best)
        df <- cbind(df,as(X[,best],"numeric"))
         LD <- ld(X,X[,best],stats="R.squared")
         maxLD <- apply(LD,1,max,na.rm=TRUE)
         drop <- unique(c(names(maxLD)[which(is.na(maxLD) | maxLD>0.9)],best))
        Xtest <- X[,setdiff(colnames(X),drop)]
        if(!ncol(Xtest))
            return(NULL)
        cond <- snp.rhs.tests(as.formula(paste(f, paste(best,collapse="+"),sep="+")),
                              snp.data=Xtest,
                              data=df) # binomial by default    
        p <- p.value(cond)
    }
    pmin <- min(p, na.rm=TRUE)
    if(pmin<stepwise.p.thr) {
        newbest <- colnames(Xtest)[ which.min(p) ]
        cs <- col.summary(Xtest[,newbest])
        message(newbest,"\tMAF=",signif(cs[1,"MAF"],2),"\tp=",format.pval(pmin))
        ret=structure(pmin,names=newbest)
        return(ret)
    }
    return(NULL)
}


(load(file.path(args$d,"GUESS/all-data.RData")))

phenotypes <- c(unique(samples(DATA)$phenotype),"iCEL","iRA")  %>%  setdiff(.,"CONTROL")
CONDITIONAL <- vector("list",length(phenotypes))
names(CONDITIONAL) <- phenotypes
ctl <- sm(DATA)[samples(DATA)$country=="UK" & samples(DATA)$phenotype=="CONTROL",]
cs <- col.summary(ctl)
## use <- cs$MAF>0.005 # avoid super rare
for(nm in phenotypes) {
    message("conditional testing for ",nm)
    if(grepl("^i",nm)) {
        ph <- sub("^i","",nm)
        countries <- samples(DATA)$country[ samples(DATA)$phenotype==ph ]  %>%  unique()
        wh <- which(samples(DATA)$phenotype %in% c(ph,"CONTROL") &
                                 samples(DATA)$country %in% countries)
    } else {
        ph <- nm
        wh <- which(samples(DATA)$phenotype %in% c(ph,"CONTROL") &
                                 samples(DATA)$country=="UK")
    }
    df <- samples(DATA)[wh,]
    df$cc <- df$affected - 1
    snp.data <- sm(DATA)[wh,]
    best <- condtest(X=snp.data, df=df)
    if(length(best))
        CONDITIONAL[[nm]] <- data.frame(snp=names(best),p=best,trait=nm,step=seq_along(best))
}
conditional <- do.call("rbind",CONDITIONAL)
save(conditional, file=file.conditional)

message("--END--")
