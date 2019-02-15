#!/usr/bin/env Rscript
library(methods)
library(randomFunctions)
library(annotSnpStats)
library(magrittr)
source("~/DIRS.txt")
setwd("~/Projects/cd4chic/GUESSFM")

                                 #library(GUESSFM)
library(devtools)
load_all("~/RP/GUESSFM")
load_all("~/RP/MTFM")


args <- getArgs(default=list(d=file.path(CVSROOT,
                                         ## "14q-69168821-69318062" # JIA not found
 ## "15q-38814377-38972177" # int only
                                         ## "2q-100544954-101038647" # v slow
                                         "10p-6030000-6220000" # even slower

                                         ## "18p-12738413-12924117"
                                         ## "10p-6030243-6169685"
                                         ## "2p-43310482-43362380" 
                                         ),cppthr=0.99,keep=FALSE))
args$d <- sub("/GUESS$","",args$d)
args

## do diseases with minp < 1e-6
library(data.table)
(load(file.path(args$d,"conditional.RData")))
conditional <- as.data.table(conditional)
conditional <- conditional[p<1e-6,]
source("~/Projects/cvs/runguess/R/functions.R")

##snpmod
thr=0.99 # or 0.9999
(load(file.path(args$d,paste0("snpmod-",sub("0.","",args$cppthr),".RData"))))

## check - this should be uk only

nd <- structure(c(12747L, 2772L, 7728L, 1214L, 4461L, 3870L, 6681L,11475L, 12041L),names= c("CONTROL", "ATD", "CEL", "JIA", "MS", "RA", "T1D","iRA","iCEL"))
dis <- intersect(names(SM2),conditional$trait) # min p < 1e-6
dis <- intersect(dis,names(nd)) # uk only
if(length(setdiff(dis,c("RA","CEL")))<2) { # otherwise count CEL+iCEL=2, for example
    system(paste("touch ",file.path(args$d,"skip-mtfm")))
    if(!interactive())
        q("no")
}

## why do this?
## if("iCEL" %in% dis && "CEL" %in% names(SM2) && !("CEL" %in% dis))
##     dis <- c(dis,"CEL")
## if("iRA" %in% dis && "RA" %in% names(SM2) && !("RA" %in% dis))
##     dis <- c(dis,"RA")
## dis <- c("GRAVES","ICOELIAC","MS","T1D")

message("found results for disease:")
print(dis)

## bestmod<- best.models(SM2[dis],cpp.thr=1.2)
## lapply(bestmod,dim)

bestmod.thr <- best.models(SM2,cpp.thr=0.99)
M <- lapply(bestmod.thr, "[[", "str") 
pr <- lapply(bestmod.thr, "[[", "prior")
abf <- lapply(bestmod.thr, "[[", "logABF")
PP <- lapply(bestmod.thr, "[[", "PP")
names(M)
p0=snpprior(n=1000,expected=3)["0"]
STR=M[dis]
ABF=abf[dis]
PP <- PP[dis]
pr=pr[dis]
message("\n\nCPP threshold = ",thr, "\n\tn.each (",paste(dis,collapse="/"),") = ",paste(sapply(M[dis],length),collapse="/"))

## todo
## f.geno <- file.path(args$d ,"GUESS","all-data.RData")
## with(DATA@samples,table(phenotype))
N0 <- nd["CONTROL"]
I0 <- structure(c(0L,0L,0L,0L,0L,0L,7443L,4814L),names= c("ATD", "CEL", "JIA", "MS", "RA", "T1D","iRA","iCEL"))
ND <- nd[dis]

ns <- nrow(SM2[[1]]@snps)
diseases <- dis
library(parallel); options(mc.cores=3)

f <- function(dis) {
kappa <- calckappa(nsnps=ns,p=2/ns,ndis=length(dis),target.odds=1)
if(length(dis)>2) {
    ndis <- length(dis)
    kappa <- c(dep=kappa,
               mid=calckappa(nsnps=ns,p=2/ns,ndis=length(SM2),target.odds=0.5^sqrt(ndis-1)/(1-0.5^sqrt(ndis-1))),
               ind=calckappa(nsnps=ns,p=2/ns,ndis=length(SM2),target.odds=0.5^(ndis-1))/(1-0.5^(ndis-1)))
}
if(length(kappa)>1) {
    mpp <- mclapply(kappa,function(k) { marginalpp(STR[dis],ABF[dis],pr[dis],k,p0,N0=N0,ND=ND[dis],nsnps=ns,I0=I0[dis]) })
    mpp <- do.call("c",mpp)
    ss <- strsplit(names(mpp),"\\.")  %>% lapply(.,rev)  %>% lapply(.,paste,collapse=".")  %>%  unlist()
    names(mpp) <- ss
} else {
    mpp <- marginalpp(STR[dis],ABF[dis],pr[dis],kappa,p0,N0=N0,ND=ND[dis],nsnps=ns,I0=I0[dis])
}
return(mpp)
}


## add groups
f.groups <- file.path(args$d,"snpmod-99-groups.RData")
(load(f.groups))

## no internationals
## if(!file.exists(file.path(args$d,"MTFM.RData"))) {
    dis <- setdiff(diseases,c("iRA","iCEL"))
    if(length(dis)>=2) {
        mpp <- f(dis)
        mpp <- lapply(mpp,function(x) {x$gr <- mod2group(rownames(x),groups=groups$groups); x})
        gpp <- lapply(mpp, function(x) { 
            tapply(1:nrow(x),x$gr,function(i)
                sum(x$shared.pp[i]))
        })
        message("saving to MTFM.RData")
        save(mpp,gpp,file=file.path(args$d,"MTFM.RData"))
    }
## }

## international 
dis <- diseases
if("iRA" %in% dis)
    dis <- setdiff(dis,"RA")
if("iCEL" %in% dis)
    dis <- setdiff(dis,"CEL")

## if(!file.exists(file.path(args$d,"iMTFM.RData")) && any(c("iRA","iCEL") %in% diseases)) {
if(any(c("iRA","iCEL") %in% diseases)) {
    message("running iMTFM for ",length(diseases)," diseases")
    impp <- f(dis)
    mpp <- lapply(impp,function(x) {x$gr <- mod2group(rownames(x),groups=groups$groups); x})
    gpp <- lapply(mpp, function(x) { 
        tapply(1:nrow(x),x$gr,function(i)
            sum(x$shared.pp[i]))
    })
    message("saving to iMTFM.RData")
    save(mpp,gpp,file=file.path(args$d,"iMTFM.RData"))
}
    
