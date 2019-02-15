#!/usr/bin/env Rscript
library(randomFunctions)
library(GUESSFM)

source("~/DIRS.txt")
## files
args <- getArgs(default=list(d=file.path(CVSROOT,"13q-42834735-43101763/GUESS")))
region <- basename(dirname(args$d))

print(args)

dirs <- list.files(args$d,full=TRUE)
files <- unlist(lapply(dirs,list.files,full=TRUE,pattern="_features"))
length(files)
if(!length(files)) {
  message("No output files found")
  if(!interactive())
    q("no")
}

ss.path <- file.path(ROOT,"FM-ss/impdata")
#files <- sub("_features.txt","",files)
#files <- files[grep("200000",files)]
#files <- grep("_50000$",files,value=TRUE)
## use R2GUESS to check runs -
## TODO - NEEDS CHECK THAT read.ess ALLOWS COVARS - NO, covars now regressed out
names(files) <- basename(dirname(files))

message("number of files found: ",length(files))
cat(files,sep="\n")

library(R2GUESS)
f.log <- function(f) { file.path(dirname(f),"log") }
f.ess <- function(f) { file.path(dirname(f),"ess.RData") }
f.plot <- function(f) { file.path(dirname(f),"ess.png")}
f.snpsum <- function(f) { file.path(dirname(f),"snpsum.csv")}
f.guessout <- function(f) { file.path(dirname(f),"guess-approx.RData")}
f.data <- file.path(args$d,"all-data.RData")

library(snpStats)
(load(f.data))

## for(f in files) {
##   if(file.exists(f.ess(f)))
##     unlink(f.ess(f))
## }


ESS <- lapply(files, function(f) {
#  if(!file.exists(f.ess(f)) || file.info(f.ess(f))$mtime < file.info(f.log(f))$mtime) {
    message("reading ess object from ",f)
    ess <- try(read.ess(f))
    if(class("ess")=="try-error") {
      print(attr(ess,"condition")$message)
      return(NULL)
    } else {
      save(ess,file=f.ess(f))
    }
  ## } else {
  ##   load(f.ess(f))
  ## }
  return(ess) })


#str(ESS)
use <- sapply(ESS,class)!="try-error"
if(any(!use)) {
  message(sum(use),"/",length(use)," runs readable:")
  print(use)
  files <- files[use]
  ESS <- ESS[use]
}
message("ESS datasets read: ",length(ESS))
message(paste(names(ESS),collapse=", "))

## generate summary plots
for(i in seq_along(files)) {
  f <- files[[i]]
  fp <- f.plot(f)
  if(!file.exists(fp) || file.info(fp)$mtime < file.info(f.log(f))$mtime ) {
    message("creating summary plot at ",fp)
    png(f.plot(f), height=1200,width=1200)
    plot(ESS[[i]])
    dev.off()
  }
}
message("now go look at the plots!")
sapply(files,f.plot)
# sapply(files, function(f) system(paste("display",f.plot(f))))

## load guess 
dd <- lapply(files,function(f) read.snpmod(dirname(f)))
pp <- pp.nsnp(dd,expected=2)

ff <- file.path(args$d,"qc3-plots.pdf")
pdf(ff)
plot(pp)
plot_diffuse(dd)
dev.off()
message("now look at ",ff)

tmp <- qc(pp)
tmp$region <- region
write.table(tmp,file=file.path(args$d,"qc3-ppnsnp.csv"),row.names=FALSE)

##             snp_id                      rs_id position exp_freq_a1  info
## 22 imm_12_54714511            12:56428244:G:A 56428244       0.062 1.000
## 23 imm_12_54721679            12:56435412:G:A 56435412       0.345 1.000
## 24             --- 12:56464088:C:CTTTTTTTTTTT 56464088       0.071 0.885
## 25             ---            12:56471334:C:T 56471334       0.076 0.979
## 26             ---            12:56454298:A:G 56454298       0.081 0.980
## 27 imm_12_54766915            12:56480648:G:A 56480648       0.349 1.000
##    certainty type info_type0 concord_type0 r2_type0 A1           A2
## 22     1.000    2      0.980         0.997    0.973  G            A
## 23     1.000    2      0.993         0.998    0.996  G            A
## 24     0.981    0     -1.000        -1.000   -1.000  C CTTTTTTTTTTT
## 25     0.996    0     -1.000        -1.000   -1.000  C            T
## 26     0.996    0     -1.000        -1.000   -1.000  A            G
## 27     1.000    2      0.987         0.996    0.990  G            A
##                   region  ph n.corr.snps n.tot.snps PP.sum.corr
## 22 12q-56369506-56482180 T1D          78        680  0.02066202
## 23 12q-56369506-56482180 T1D          78        680  0.02066202
## 24 12q-56369506-56482180 T1D          78        680  0.02066202
## 25 12q-56369506-56482180 T1D          78        680  0.02066202
## 26 12q-56369506-56482180 T1D          78        680  0.02066202
## 27 12q-56369506-56482180 T1D          78        680  0.02066202

maxld <- vector("list",length(files))
names(maxld) <- names(files)
snps.to.check <- maxld
for(ph in names(files)) {
  message(ph)
  dfile <- file.path(dirname(files[[ph]]),"data.RData")
  (load(dfile))  
  maxld[[ph]] <- qc(dd[[ph]],snp.data)
  y <- maxld[[ph]]$SNPs
  y <- y[y>=0.2]
  if(!is.null(y))
    snps.to.check[[ph]] <- cbind(DATA@snps[names(y),],
                                 region=region,
                                 ph=ph,
                                 as.data.frame(t(maxld[[ph]]$summary)),
                                 snp.freq.in.bad.models=y)
}

snps.to.check <- do.call("rbind",snps.to.check)
snps.to.check <- subset(snps.to.check,sumPP.models.corr.snps>0.01)
if(nrow(snps.to.check))
  write.table(snps.to.check,file=file.path(args$d,"qc-snps.csv"),row.names=FALSE)
