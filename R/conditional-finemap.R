#!/usr/bin/env Rscript
## follows on from GUESS-qc.R, may need some code from that inserted

library(randomFunctions)
library(annotSnpStats)
library(magrittr)
source("~/DIRS.txt")
setwd("~/Projects/cd4chic/GUESSFM")

## library(devtools)
## load_all("~/RP/GUESSFM")
args <- getArgs(default=list(d=file.path(CVSROOT,
                                         "6q-159322326-159541830"
                                        #14q-101290463-101328739"
                                         ## "1q-172650685-172940450"
                                         ),keep=FALSE))
args$d <- sub("/GUESS$","",args$d)
args

## source("R/GUESS-common.R")
## manhattans
file.conditional <- file.path(args$d,"conditional.RData")
file.pp <- file.path(args$d,"conditional-finemap.RData")
(load(file.conditional))
conditional <- subset(conditional,p<1e-6) # or no point fine mapping
if(nrow(conditional)==0){
    message("nothing to finemap")
    q("no")
}
conditional$snp  %<>%  as.character()
sc <- split(conditional,conditional$trait,drop=TRUE)

library(coloc)
calc.pp <- function(X,df,best=NULL, ...) {
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
        df <- cbind(df,as(X[,best],"numeric"))
        Xtest <- X[,setdiff(colnames(X),best)]
        if(!ncol(Xtest))
            return(NULL)
        cond <- snp.rhs.tests(as.formula(paste(f, paste(best,collapse="+"),sep="+")),
                              snp.data=Xtest,
                              data=df) # binomial by default    
        p <- p.value(cond)
    }
    pp <- finemap.abf(dataset=list(pvalues=p,
                                   N=nrow(df),
                                   MAF=col.summary(Xtest)[,"MAF"],
                                   type="cc",
                                   s=mean(df$cc),
                                   snp=colnames(Xtest)),
                      p1=1/ncol(Xtest))
    pp <- with(subset(pp,snp!="null"), structure(SNP.PP, names=as.character(snp)))
    pp <- pp/sum(pp) # prior that exactly one snp is causal
    return(pp)
}


(load(file.path(args$d,"GUESS/all-data.RData")))

library(parallel)
options(mc.cores=4)

phenotypes <- names(sc) 
PP <- vector("list",length(phenotypes))
names(PP) <- phenotypes
ctl <- sm(DATA)[samples(DATA)$country=="UK" & samples(DATA)$phenotype=="CONTROL",]
cs <- col.summary(ctl)
## use <- cs$MAF>0.005 # avoid super rare
for(nm in phenotypes) {
    message("finemapping ",nm)
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
    best <- NULL
    PP[[nm]] <- mclapply(1:nrow(sc[[nm]]), function(j) {
        if(sc[[nm]][j,"step"]!=j)
            return(NULL) # only fine map when every snp to this point has p < 1e-6
        calc.pp(X=snp.data,df=df,best=if(j==1) { NULL } else { sc[[nm]][1:(j-1),"snp"] })
    })
}
save(PP, file=file.pp)

message("--END--")
