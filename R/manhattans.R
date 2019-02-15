#!/usr/bin/env Rscript

library(randomFunctions, quietly=TRUE)
library(annotSnpStats, quietly=TRUE)
library(snpStats)
library(magrittr)
source("~/DIRS.txt") ## gets CVSROOT, ROOT

## differs from analysis in cd4chic paper by dropping SNPs with MAF in controls < 0.001 rather than < 0.01.

args <- getArgs(default=list(d=file.path(CVSROOT,
                                        "10p-6030243-6169685"
                                         ## "4q-122973062-123565302"
#"5q-158515199-158836107"
 # "1p-117030982-117142239"
                                        #"6p-34507311-35327577/GUESS"
                                         ),keep=FALSE,cppthr=0.99)) #10p-6030243-6169685
args$d <- sub("/GUESS","",args$d)
args

## save this version of data for all analyses
(load(file=file.path(args$d,"GUESS/all-data.RData")))

## manhattans
phenotypes <- unique(samples(DATA)$phenotype)  %>%  setdiff(.,"CONTROL")

## uk
ss <- lapply(phenotypes, function(ph) {
    use <- which(samples(DATA)$phenotype %in% c("CONTROL",ph) & samples(DATA)$country == "UK")
    df <- samples(DATA)
    df$cc <- df$affected - 1
    snp.rhs.tests(cc ~ PC1 + PC2 + PC3 + PC4 ,snp.data=sm(DATA)[use,],data=df) %>% p.value(.)
}) %>% do.call("cbind",.)

## international
iss <- lapply(c("RA","CEL"), function(ph) {
    countries <- samples(DATA)$country[ samples(DATA)$phenotype==ph ]  %>%  unique()
    use <- which(samples(DATA)$phenotype %in% c("CONTROL",ph) & samples(DATA)$country %in% countries)
    df <- samples(DATA)
    df$cc <- df$affected-1
    snp.rhs.tests(cc ~ PC1 + PC2 + PC3 + PC4 + country,snp.data=sm(DATA)[use,],data=df) %>% p.value(.)
}) %>% do.call("cbind",.)
head(iss)    

colnames(ss) <- phenotypes
colnames(iss) <- paste0("i",c("RA","CEL"))
ss <- cbind(ss,iss)
summary(ss)

apply(ss,2,min)
colnames(ss)  %<>%  paste0("p.",.)
head(DATA@snps)
ss <- cbind(DATA@snps,ss)
head(ss)

write.table(ss,file=file.path(args$d,"single-snp-pvalues.csv"))
