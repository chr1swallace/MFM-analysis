#!/usr/bin/env Rscript
library(randomFunctions, quietly=TRUE)
library(annotSnpStats, quietly=TRUE)
library(snpStats)
library(magrittr)
source("~/DIRS.txt") ## gets CVSROOT, ROOT

## differs from analysis in cd4chic paper by dropping SNPs with MAF in controls < 0.001 rather than < 0.01.

args <- getArgs(defaults=list(
                  ## file=file.path(IMPROOT,"imputed-8q-126530818-126550818.RData"), #3p-45929800-46650993.RData"),
                  ## file=file.path(IMPROOT,"imputed-10p-6030000-6220000.RData"), #3p-45929800-46650993.RData"),
                  file=file.path(IMPROOT,"imputed-2q-204446380-204816382.RData"), #3p-45929800-46650993.RData"),
                  overwrite=FALSE,
                  MAF=0.001,
                  #region=75,
                  tag.r2=0.98, nsweep=50000, nsave=5000, nchains=5),
                numeric=c("tag.r2","nsweep","nsave","nchains","region"))

setwd(CVSRUN)

## files to read
d <- IMPROOT

f.geno <- file.path(IMPROOT,basename(args$file))

## create vector, snps.drop
source("R/drop-imputed-snps.R")
samples.drop <- c("5442596077_R01C02","5767876010_R03C02")

## files to create
##  library(devtools, quietly=TRUE)
## install_github("chr1swallace/GUESSFM")
#library(GUESSFM)
library(devtools)
load_all("~/RP/GUESSFM",quiet=TRUE)
gd.root <- file.path(CVSROOT,"GUESS")
file.root <- gsub("impute-|imputed-|.imputed.RData|.RData","",basename(f.geno))

if(file.root=="2q-204446380-204816382") { ## check soumya snp
    file.root <- paste0(file.root,"+rs117701653")
} 
gd <- file.path(CVSROOT,file.root,"GUESS")
if(!file.exists(gd))
  dir.create(gd,recursive=TRUE)
fd <- function(d,f) {
  if(!file.exists(d))
    dir.create(d,recursive=TRUE)
  file.path(d,f)
}

## load data and phenotypes
message("loading genotype data: ",f.geno)
print(load(f.geno))
colnames(samples(DATA)) <- sub("samples.","",colnames(samples(DATA)))

message(ncol(DATA)," SNPs in region")
phenotype(DATA) <- "affected"

## drop bad snps
snps.drop <- setdiff(snps.drop,"rs117701653.204628795.A.C")
wh <- which(annotSnpStats::snps(DATA)$rs_id %in% snps.drop |
            annotSnpStats::snps(DATA)$snp_id %in% snps.drop |
            rownames(annotSnpStats::snps(DATA)) %in% snps.drop)
message("dropping ",length(wh)," SNPs listed in snps.drop")
if(length(wh)) {
  DATA <- DATA[,-wh]
}

## drop bad samples
wh <- which(samples(DATA)$orig.rownames %in% samples.drop)
message("dropping ",length(wh)," SAMPLEs listed in samples.drop")
if(length(wh)) {
  DATA <- DATA[-wh,]
}

## filter snps
snps.drop3 <- rownames(DATA@snps)[(DATA@snps$info<0.3 & DATA@snps$type==0) |
                                  (DATA@snps$certainty<0.98 & DATA@snps$type==0) ]
wh <- which(annotSnpStats::snps(DATA)$rs_id %in% snps.drop3 |
            annotSnpStats::snps(DATA)$snp_id %in% snps.drop3 |
            rownames(annotSnpStats::snps(DATA)) %in% snps.drop3)
message("dropping ",length(wh)," SNPs with info<0.3 or certainty<0.98")
if(length(wh)) {
  DATA <- DATA[,-wh]
}

## rare snps - MAF < 0.005 in UK controls
use <- with(samples(DATA),country=="UK" & phenotype=="CONTROL")
cs <- col.summary(DATA[ which(use), ])
drop <- which(cs$MAF < 0.005 | is.na(cs$z.HWE))
if(length(drop)) {
    message("dropping snps with MAF < 0.005 in UK controls: ",length(drop))
    DATA <- DATA[,-drop]
}

## add pcs
f.samples <- file.path(ROOT,"FM-pca/PCS+.RData")
(load(f.samples))
rownames(pcs) <- make.names(rownames(pcs))
sinfo <- samples(DATA)
dim(sinfo)
dim(pcs)
missing <- pcs[ setdiff(rownames(pcs),rownames(sinfo)),] # ok - they are PBC (one region only) and 2 GRAVES that were deliberately dropped
DATA <- DATA[ rownames(DATA) %in% rownames(pcs),]
sinfo <- samples(DATA)
length(unique(rownames(sinfo)))
length(unique(rownames(pcs)))
sinfo2 <- pcs[rownames(sinfo),]
table(sinfo$phenotype,sinfo2$phenotype)
table(sinfo$country,sinfo2$country,sinfo$file) ## ICOELIAC has country messed up in samples(DATA)

## drop duplicate UK controls
wh <- which(sinfo2$phenotype=="CONTROL" & sinfo2$country=="UK" & sinfo2$file!="T1D")
with(sinfo2[wh,],table(country,file))
with(sinfo2[wh,],table(phenotype,file))
DATA <- DATA[ -wh, ]
sinfo2 <- sinfo2[ -wh, ]

samples(DATA) <- sinfo2
DATA <- DATA[ samples(DATA)$phenotype!="COELIAC", ]
transph <- function(ph) { # go from Nick's names to mine
  ## ph <- tolower(ph)
  ph <- sub("COELIAC","CEL",ph)
  ph <- sub("ICEL","CEL",ph)
  ph <- sub("GRAVES","ATD",ph)
  return(ph)
}
## (load(file.path(args$d,"GUESS/all-data.RData")))
tmp <- samples(DATA)
tmp$phenotype <- transph(tmp$phenotype)
samples(DATA) <- tmp


## loop over phenotypes
phenotypes <- setdiff(unique(samples(DATA)$phenotype),"CONTROL")


## look at differential certain call rates
PH <- c("CONTROL",phenotypes)
CS <- lapply(PH,function(ph) {
  col.summary(DATA[ which(samples(DATA)$phenotype==ph), ])
})
todrop <- function(x,y) {
  abs(log(x)-log(y))>0.05 | (y>0.99 & x<0.98) | (x>0.99 & y<0.98)
}
w1 <- which(PH=="CONTROL")
keep <- rep(TRUE,ncol(DATA))
for(ph in setdiff(PH,"CONTROL")) {
  w2 <- which(PH==ph)
  x <- CS[[w1]][,"Certain.calls"]
  y <- CS[[w2]][,"Certain.calls"]
  keep <- keep & !todrop(x,y)  
}
if(any(!keep)) {
  message("dropping ",sum(!keep)," SNPs with differential certain calls")
  DATA <- DATA[,which(keep)]
}

## regions/phenotypes to drop
rp <- read.table("drop-regions-phenotypes.txt",as.is=TRUE,header=TRUE)
wh <- which(rp$region==basename(gd) & rp$ph %in% phenotypes)
if(length(wh)) {
  message("dropping region/phenotypes from drop-regions-phenotypes.txt: ",length(wh))
  print(rp[wh,])
  phenotypes <- setdiff(phenotypes,rp[wh,"ph"])
}
if(!length(phenotypes))
    stop()

## non-missing phenotype, UK only
cc <- samples(DATA)$affected - 1
country <- samples(DATA)$country
use <- complete.cases(cc,country) #& country=="UK"
## iDATA <- DATA[complete.cases(cc,country) & samples(DATA)$phenotype %in% c("CONTROL","ICOELIAC","RA"),]

## extract UK subset of ICOELIAC
## library(data.table)
## clst <- fread("/scratch/wallace/IChip/Trynka2011_NatureGenetics_Immunochip/final.clst")
## head(clst)
## setnames(clst,c("pedigree","member","clst"))
## m <- merge(samples(DATA),clst,by=c("pedigree","member"))
## dim(clst)
## dim(m)
## table(clst$clst)
## table(m$clst)
## head(m)
## with(m,table(file))
## with(m,table(clst,affected))
## m <- m[m$affected==2 & m$clst==2,] # uk cases only
## use <- use & (samples(DATA)$phenotype!="ICOELIAC" | samples(DATA)$orig.rownames %in% m$orig.rownames)
## table(use)
## table(samples(DATA)$phenotype,use)

cc <- cc[which(use)]
country <- country[which(use)]
DATA <- DATA[which(use),]
message("subsetting to ",length(which(use))," SAMPLEs with phenotype information")
print(table(cc,country))

with(samples(DATA),table(file,phenotype))
with(samples(DATA),table(phenotype,cc))
with(samples(DATA),table(phenotype,country))

## save this version of data for all analyses
save(DATA,file=file.path(gd,"all-data.RData"))

## what is min p?  Only go ahead if mip < 10-5
ss <- lapply(phenotypes, function(ph) {
    countries <- samples(DATA)$country[ samples(DATA)$phenotype==ph ]  %>%  unique()
    use <- which(samples(DATA)$phenotype %in% c("CONTROL",ph) & samples(DATA)$country %in% countries)
    df <- samples(DATA)
    df$cc <- cc
    snp.rhs.tests(cc ~ PC1 + PC2 + PC3 + PC4 + country,snp.data=sm(DATA)[use,],data=df) %>% p.value(.)
}) %>% do.call("cbind",.)
head(ss)    
head(DATA@snps)

pvars <- paste("p",phenotypes,sep=".")
colnames(ss) <- pvars
print(minp <- apply(ss[,pvars],2,min,na.rm=TRUE))
ss <- cbind(DATA@snps,ss)
## write.table(ss,file=file.path(CVSROOT,file.root,"single-snp-pvalues.csv"))

phenotypes <- phenotypes[ minp < 1e-6 ]
sapply(phenotypes, function(ph) message("INTENDED ",basename(args$file)," ",ph))

phenotypes <- transph(phenotypes)
skips <- setdiff(c("T1D", "RA", "MS", "JIA", "CEL", "ATD"), phenotypes)
for(ph in skips) 
    paste("touch",file.path(gd,paste0(ph,"-skip"))) %>% system()

if(all(minp>1e-6) || !length(phenotypes)) {
  message()
  message("--------------------------------------------")
  message("no significant assoc to the region, quitting")
  message("--------------------------------------------")
  paste("touch",file.path(gd,"skip")) %>% system()
  message()
  q("no")
}

## if(length(phenotypes<6)) {
## }

## check if output already exists for selected phenotypes
f.data <- file.path(gd,phenotypes,"data.RData")
phenotypes <- phenotypes[ !file.exists(f.data) ]
if(!length(phenotypes))
    stop("All required output files exist already")

message("Running GUESS for ",length(phenotypes)," diseases: ",paste(phenotypes,collapse=", "))
  options(scipen=1000000000)
## tags
f.tags <- fd(file.path(gd),"tags.RData")
sinfo <- samples(DATA)
if(args$tag.r2<1 && !file.exists(f.tags)) {
  message("Preparing to tag.")
    wh.DATA <- which(samples(DATA)$phenotype == "CONTROL" & samples(DATA)$country=="UK")
    snp.data1 <- as(DATA,"SnpMatrix")[-wh.DATA,]
    snp.data0 <- as(DATA,"SnpMatrix")[wh.DATA,]
    message("Input matrix has ",ncol(snp.data0)," SNPs.")
    cs0 <- col.summary(snp.data0)
    cs1 <- col.summary(snp.data1)
    wh <- which(is.na(cs0[,"z.HWE"]) |
                cs0[,"MAF"]<0.005 |
                cs0[,"Call.rate"]<0.99 | cs1[,"Call.rate"]<0.99 |
                cs0[,"Certain.calls"]<0.75 | cs1[,"Certain.calls"]<0.75 |
                abs(cs0[,"z.HWE"])>4)
    if(length(wh)) {
      message("Dropping ",length(wh)," SNPs with |z.HWE|>5, MAF < 0.005 in controls or call rate <0.99 in cases or controls")
      snp.data0 <- snp.data0[,-wh]
    }
    tags <- tag(snp.data0, method="single", tag.threshold=args$tag.r2)
    message("after tagging at ",args$tag.r2,", matrix now has ",length(unique(tags(tags)))," SNPs.")
    save(tags, file=f.tags)
} else {
    (load(f.tags))
}


coms <- structure(rep("",length(phenotypes)),
                  names=c(phenotypes))
if("CEL" %in% names(coms))
    coms[c("iCEL")] <- ""
if("RA" %in% names(coms))
    coms[c("iRA")] <- ""

sinfo <- samples(DATA)
cc <- sinfo$affected - 1
for(ph in phenotypes) {
  message(ph)

  f.data <- function(ph) fd(file.path(gd,ph),"data.RData")
  f.par <- function(ph) fd(file.path(gd,ph), "par.xml")
  ## if(file.exists(f.par(ph)) & !args$overwrite) {
  ##   message("output file already exists, skipping: ",f.par(ph))
  ##   next
  ## }
  wh.DATA <- which(sinfo$phenotype %in% c(ph,"CONTROL") &
                   sinfo$country=="UK" &
                   complete.cases(sinfo[,c("PC1","PC2","PC3","PC4")]))
  cc.ph <- cc[wh.DATA]
  snp.data <- as(DATA,"SnpMatrix")[wh.DATA,]
  snp.data <- snp.data[, unique(tags(tags))]

  ## PCs
  covars <- sinfo[wh.DATA,c("PC1","PC2","PC3","PC4")]
  fix <- which(sapply(covars,class)=="matrix")
  if(length(fix)) {
    for(j in fix) 
      covars[,j] <- as.vector(covars[,j])
  }
  
  ## save data
  save(snp.data,cc.ph,covars,file=f.data(ph))

  message(ph)
  message("Cases/controls: ", sum(cc.ph==1), " ", sum(cc.ph==0))
  message("SNPs: ",ncol(snp.data))
  #message("Tags: ",length(unique(tags(tags))))

  ## command=paste("qCom.sh -L", file.path(gd,ph,"2log"), "/home/chrisw/local/bin/GUESS")
  ## command=paste("/home/chrisw/local/bin/GUESS")
  coms[[ph]] <- run.bvs(X=snp.data,Y=cc.ph,covars=covars,
                       gdir=file.path(gd,ph),nsweep=args$nsweep, family="binomial",tag.r2=NA,
                       nsave=args$nsave,nchains=args$nchains,nexp=2,run=FALSE)

 if(ph %in% c("RA","CEL")) {
        ## international+UK
     iph <- paste0("i",ph)
     countries <- samples(DATA)$country[ samples(DATA)$phenotype==ph ]  %>%  unique()
     wh.DATA <- which(sinfo$phenotype %in% c(ph,"CONTROL") &
                      sinfo$country %in% countries &
                      complete.cases(sinfo[,c("PC1","PC2","PC3","PC4")]))
     cc.ph <- cc[wh.DATA]
     snp.data <- as(DATA,"SnpMatrix")[wh.DATA,]
     snp.data <- snp.data[, unique(tags(tags))]
     
     ## PCs
     covars <- sinfo[wh.DATA,c("country","PC1","PC2","PC3","PC4")]
     fix <- which(sapply(covars,class)=="matrix")
     if(length(fix)) {
         for(j in fix) 
             covars[,j] <- as.vector(covars[,j])
     }
  
     ## save data
     save(snp.data,cc.ph,covars,file=f.data(iph))

     message(iph)
     message("Cases/controls: ", sum(cc.ph==1), " ", sum(cc.ph==0))
     message("SNPs: ",ncol(snp.data))
     
     coms[[iph]] <- run.bvs(X=snp.data,Y=cc.ph,covars=covars,
                            gdir=file.path(gd,iph),nsweep=args$nsweep, family="binomial",tag.r2=NA,
                            nsave=args$nsave,nchains=args$nchains,nexp=2,run=FALSE)
    } 
}

