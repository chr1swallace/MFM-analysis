#!/usr/bin/env Rscript
## follows on from GUESS-qc.R, may need some code from that inserted

library(randomFunctions)
library(annotSnpStats)
library(magrittr)
source("~/DIRS.txt")
setwd("~/Projects/cvs/runguess")


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
                                          "2q-204446380-204816382+rs117701653" # out of memory - 2883 snps!
                                        ##"6p-34507311-35327577/GUESS"
     #     "10p-6030243-6169685"
                                         ## "6p-34507311-35327577"
                                         ),keep=FALSE,cppthr=0.99)) #10p-6030243-6169685
#args <- getArgs(default=list(d=file.path(ROOT,"FM-GUESS/3p-45929800-46650993"),keep=FALSE,cppthr=0.99))
## failed
##22q-21712414-22022005 - LD, fixed
##4q-122973062-123565302
##3p-45929800-46650993
#args <- getArgs(default=list(d=file.path(ROOT,"FM-GUESS/6q-137882875-138275085"),keep=FALSE))
args$d <- sub("/GUESS","",args$d)
args
CPPTHR <- args$cppthr - 0.01 ## try
if(grepl("3p-45929800-46650993",args$d))
    CPPTHR <- CPPTHR - 0.05

source("R/GUESS-common.R")

library(ggplot2)
theme_set(theme_bw())
library(gridExtra)

dirs <- list.files(file.path(args$d,"GUESS"),full=TRUE)
files <- unlist(lapply(dirs,list.files,full=TRUE,pattern="_features"))
if(!length(files)) {
    message("No output files found")
    if(!interactive())
        q("no")
}
files <- sub("_sweeps_features.txt","",files)
names(files) <- basename(dirname(files))

message("number of files found: ",length(files))
cat(files,sep="\n")

## load guess 
dd <- lapply(files,read.snpmod)
names(dd) <- basename(dirname(files))
dd

## expand tags
(load(file.path(args$d,"GUESS","tags.RData")))

library(parallel)
NCORES <- if(interactive()) {1} else {4}
options(mc.cores=NCORES) 

tmpfile <- file.path(args$d,"GUESS","tmp-dx.RData")
if(file.exists(tmpfile)) {
    (load(tmpfile))
} else {
    dx <- mclapply(dd,expand.tags,tags=tags)
    save(dx,file=tmpfile)
}
rm(dd)
## bm <- best.models(dx)

##' FIRST, STOP HERE IF NULL MODEL HAS > 50% POSTERIOR ACROSS ALL TRAITS
dropmodels <- function(x,minpp=0.5) {
    if(!any(x@models$str==""))
        return(FALSE)
    w <- which(x@models$str=="")
    x@models[w,"PP"]
}

nullpp <- lapply(dx, dropmodels) %>% unlist()
message("null model posterior")
nullpp
message("null model >50% of posterior?")
print(drop <- nullpp>0.5)

if(all(drop)) {
    system(paste0("touch ",file.path(args$d,"GUESS/skip")))
    stop("all diseases null")
}

## GIVEN SOME TRAITS HAVE MODELS OF INTEREST, DROP ONLY THOSE THAT HAVE REALLY NO SUPPORT
drop <- nullpp>0.8

if(any(drop)) {
    drop <- names(drop)[drop]    
    message("dropping ",paste(drop,collapse=" "))
    ## dd <- dd[setdiff(names(dd),drop)]
    dx <- dx[setdiff(names(dd),drop)]
    files <- files[setdiff(names(files),drop)]
}
    
library(reshape)
library(plyr)

message("expanding for ",length(dx)," diseases:")
print(names(dx))

## refits
(load(f.geno))
message("region loaded from ",f.geno)
message(ncol(DATA)," SNPs in region")

library(speedglm)

## this is where we trim models, to capture CPPTHR (args$cpp.thr) proportion of the cumulative postprob
tmpfile <- file.path(args$d,"GUESS","tmp-union.RData")
if(file.exists(tmpfile)) {
    (load(tmpfile))
} else {
    tmp <- best.models(dx,cpp.thr=CPPTHR)
    lapply(tmp,nrow)
    if(length(tmp)>1) {
        union <- model.union(tmp[[1]]$str,tmp[[2]]$str,detail=FALSE)
        if(length(tmp)>2)
            for(j in 3:length(tmp))
                union <- model.union(union,tmp[[j]]$str)
        length(union)
    } else {
        union <- tmp[[1]]$str
    }
    save(union,file=tmpfile)
}
rm(tmp,dx); gc()

## if(grepl("22q-21712414-22022005",args$d))
##     union <- setdiff(union,c("rs5749495.21940310.A.G%rs5754234.21942978.A.T","rs34463645.21941164.CTT.C%rs5754234.21942978.A.T")) # in LD, causes errors in ICOELIAC with tag SNP, and low PP in COELIAC

phenotypes <- unique(intersect(sub("i","",names(files)),c("T1D","RA","MS","JIA","CEL","ATD")))
fits <- structure(vector("list",length(phenotypes)), names=phenotypes)
if("CEL" %in% names(fits))
    fits[c("iCEL")] <- NA
if("RA" %in% names(fits))
    fits[c("iRA")] <- NA

## SDATA <- as(DATA,"SnpMatrix")
snp.data0 <- sm( DATA[ samples(DATA)$phenotype=="CONTROL" & samples(DATA)$country=="UK",
                      tags@.Data ])
snp.alldata <-  sm(DATA[,tags@.Data ])
for(ph in phenotypes) {
    message(ph)
    of <- file.path(args$d,"GUESS",paste0("fit-",ph,".RData"))
    if(file.exists(of)) {
        fits[[ph]] <- readRDS(of)
    } else {
        use <- with(samples(DATA), which(country=="UK" & phenotype %in% c("CONTROL",ph)))
        snp.data <- snp.alldata[use,]
        covars <- samples(DATA)[use,c("PC1","PC2","PC3","PC4")]
        cc.ph <- samples(DATA)$affected[use] - 1
        rs <- row.summary(snp.data)
        use <- which(rs[,"Call.rate"]==1 & complete.cases(covars) & complete.cases(cc.ph))
        fit <- abf.calc(y=cc.ph[use],x=snp.data[use,],models=union,family="binomial",
                        snp.data=snp.data0,
                        q=covars[use,])
        saveRDS(fit,file=of)
        fits[[ph]] <- fit
    }
    if(ph %in% c("RA","CEL")) {
        ## international+UK
        iph <- paste0("i",ph)
        of <- file.path(args$d,"GUESS",paste0("fit-",iph,".RData"))
        if(file.exists(of)) {
            fits[[iph]] <- readRDS(of)
        } else {
            countries <- samples(DATA)$country[ samples(DATA)$phenotype==ph ]  %>%  unique()
            use <- with(samples(DATA), which(country %in% countries & phenotype %in% c("CONTROL",ph)))
            snp.data <- snp.alldata[use,]
            isnp.data0 <- sm( DATA[ samples(DATA)$phenotype=="CONTROL" &
                                                 samples(DATA)$country %in% countries,
                                   tags@.Data ])
            covars <- samples(DATA)[use,c("PC1","PC2","PC3","PC4","country")]
            cc.ph <- samples(DATA)$affected[use] - 1
            rs <- row.summary(snp.data)
            use <- which(rs[,"Call.rate"]==1 & complete.cases(covars) & complete.cases(cc.ph))
            fit <- abf.calc(y=cc.ph[use],x=snp.data[use,],models=union,family="binomial",
                            snp.data=isnp.data0,
                            q=covars[use,])
            saveRDS(fit,file=of)
            fits[[iph]] <- fit
        }
        ## ## only international
        ## iph <- paste0("o",ph)
        ## of <- file.path(args$d,"GUESS",paste0("fit-",iph,".RData"))
        ## if(file.exists(of)) {
        ##     fits[[iph]] <- readRDS(of)
        ## } else {
        ##     countries <- samples(DATA)$country[ samples(DATA)$phenotype==ph ]  %>%  unique()  %>%  setdiff(.,"UK")
        ##     use <- with(samples(DATA), which(country %in% countries & phenotype %in% c("CONTROL",ph)))
        ##     snp.data <- snp.alldata[use,]
        ##     isnp.data0 <- sm( DATA[ samples(DATA)$phenotype=="CONTROL" &
        ##                                          samples(DATA)$country %in% countries,
        ##                            tags@.Data ])
        ##     covars <- samples(DATA)[use,c("PC1","PC2","PC3","PC4","country")]
        ##     cc.ph <- samples(DATA)$affected[use] - 1
        ##     rs <- row.summary(snp.data)
        ##     use <- which(rs[,"Call.rate"]==1 & complete.cases(covars) & complete.cases(cc.ph))
        ##     fit <- abf.calc(y=cc.ph[use],x=snp.data[use,],models=union,family="binomial",
        ##                     snp.data=isnp.data0,
        ##                     q=covars[use,])
        ##     saveRDS(fit,file=of)
        ##     fits[[iph]] <- fit
        ## }
    }
}
allph <- names(files)

## make snpmods
SM2 <- lapply(fits,function(x) abf2snpmod(x,expected=2,overdispersion=1,nsnps=length(tags@.Data)))
of <-file.path(args$d,paste0("snpmod-",
                                      sub("0.","",as.character(args$cppthr)),
                                  ".RData"))
message("saving SM2 to ",of)
save(SM2,file=of)
