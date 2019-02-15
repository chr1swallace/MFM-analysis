#!/usr/bin/env Rscript
library(randomFunctions)
library(magrittr)
source("~/DIRS.txt")
args <- getArgs(default=list(d=file.path(CVSROOT,
                                         "10p-6030000-6220000"
                                         ## "2q-181794257-182263951"
 #"2q-204446380-204816382+rs117701653" # out of memory - 2883 snps!
                                          ## "2q-204446380-204816382" # out of memory - 2883 snps!
                                         ## "3q-188069360-188139629"
## "5q-55414956-55447909" #"12q-56623347-56753822"
                                        #"3p-45929800-46650993"                                         #"20q-44594228-44784336"
                                        # "4q-122973062-123565302"
                                        #"6q-167344898-167547954"
                                        ## "19q-49094856-49278082"
                                        ## "10p-6030243-6169685"
                                         ),keep=FALSE,cppthr=0.99))
                                        #args <- list(d="/home/cew54/scratch/FM-GUESS/5p-40286967-40818088",cppthr=0.99) #2q-204446380-204816382")
                                        #args <- getArgs(defaults=args)
prefix <- paste0("snpmod-",sub("0.","",as.character(args$cppthr)))

library(annotSnpStats)
library(magrittr)
setwd("~/Projects/cvs/runguess")
library(devtools)
load_all("~/RP/GUESSFM")

                                        #args <- getArgs(default=list(d=file.path(ROOT,"FM-GUESS/10p-6030243-6169685"),keep=FALSE))
args
list.files(args$d)
source("R/GUESS-common.R")

################################################################################
## pp, mppi, r

library(ggplot2)
theme_set(theme_bw())

(load(file=file.path(args$d,paste0(prefix,".RData")))) ## SM2

### data to Calculate r2
dfile <- file.path(args$d,"GUESS","all-data.RData")
(load(dfile))
use <- with(samples(DATA), which(country=="UK" & phenotype=="CONTROL"))
snp.data <- sm(DATA)[use,]
dim(snp.data)

################################################################################

## cluster all diseases together

gfile <-file.path(args$d,paste0(prefix,"-groups.RData"))
## if(file.exists(gfile)) {
##     load(gfile)
## } else {
    groups <- group.multi(SM2,snp.data)
    ## order by position of tags
    o <- order(DATA@snps[tags(groups$groups),"position"])
    groups$groups@tags <- groups$groups@tags[o]
    groups$groups@.Data <- groups$groups@.Data[o]
## }

if(basename(args$d)=="2q-181794257-182263951") {
    s1 <- grep("rs12620692",unlist(groups$groups),value=TRUE) # B
    s2 <- grep("rs12693272",unlist(groups$groups),value=TRUE) # C
    snpin(c(s1,s2),groups$groups) # groups 2, 3
    ## ld = 0.9 between 2 and 8 (group D)
    groups$groups[[2]]
    groups$groups[[3]]
    lapply(SM2, check.merge, groups$groups[c(2,3)])
    groups$groups <- groups.merge(groups$groups,groups$groups@tags[c(2,3)])
}
if(basename(args$d)=="6q-126442383-127412057") {
    s1 <- grep("rs34452120",unlist(groups$groups),value=TRUE) # B
    s2 <- grep("rs34137445",unlist(groups$groups),value=TRUE) # C
    snpin(c(s1,s2),groups$groups) # groups 2, 3
    ## ld = 0.9 between 2 and 8 (group D)
    groups$groups[[2]]
    groups$groups[[1]]
    lapply(SM2, check.merge, groups$groups[c(2,1)])
    groups$groups <- groups.merge(groups$groups,groups$groups@tags[c(2,1)])
}

## refine
if(basename(args$d)=="10p-6030000-6220000") {
    s1 <- grep("rs7090530",unlist(groups$groups),value=TRUE)
    s2 <- grep("rs706779",unlist(groups$groups),value=TRUE)
    snpin(c(s1,s2),groups$groups) # group 6
    groups$groups@tags <- c(groups$groups@tags,s2)
    groups$groups@.Data <- c(groups$groups@.Data,list(c(s1,s2)))
    groups$groups@.Data[[6]] <- setdiff(groups$groups@.Data[[6]], c(s1,s2))

    ## ld = 0.9 between 2 and 8 (group D)
    groups$groups[[2]]
    groups$groups[[8]]
    lapply(SM2, check.merge, groups$groups[c(2,8)])
    groups$groups <- groups.merge(groups$groups,groups$groups@tags[c(2,8)])
}
if(basename(args$d)=="20p-1497197-1689461") {
     s1 <- grep("rs202535",unlist(groups$groups),value=TRUE)
     s2 <- grep("rs202536",unlist(groups$groups),value=TRUE)
     snpin(c(s1,s2),groups$groups) # groups 2, 3
     groups$groups@tags[c(2,3)]
     groups$groups@.Data[c(2,3)]
     # nothing to do
}
if(basename(args$d)=="2q-100544954-101038647") {
     s1 <- grep("rs10209110",unlist(groups$groups),value=TRUE)
     s2 <- grep("rs4851252",unlist(groups$groups),value=TRUE) 
     snpin(c(s1,s2),groups$groups) # groups 1, 3
     groups$r2[ groups$groups@.Data[[1]], groups$groups@.Data[[3]]] # all about 0.5
     ## by RA sub pop
     pops <- with(samples(DATA), country)  %>%  unique()
     for(p in pops) {
         w <- which(with(samples(DATA),country==p & phenotype=="CONTROL"))
         message(p)
         print(summary(c(ld(DATA[w,groups$groups@.Data[[1]]], DATA[w, groups$groups@.Data[[3]]], stat="R.squared"))))
     }
     # nothing to do
}
if(basename(args$d)=="2q-204446380-204816382+rs117701653") {
     s1 <- grep("rs34050244",unlist(groups$groups),value=TRUE) # G
     s2 <- grep("rs3087243",unlist(groups$groups),value=TRUE) # K
     snpin(c(s1,s2),groups$groups) # groups 7
     groups$r2[ groups$groups@.Data[[8]], groups$groups@.Data[[11]]] # 0.75 - 0.89
     ## by RA sub pop
     pops <- with(samples(DATA), country)  %>%  unique()
     for(p in pops) {
         w <- which(with(samples(DATA),country==p & phenotype=="CONTROL"))
         message(p)
         print(summary(c(ld(DATA[w,groups$groups@.Data[[8]]], DATA[w, groups$groups@.Data[[11]]], stat="R.squared"))))
     } # 0.6 min - 0.94 max
     groups$groups[[8]]
    groups$groups[[11]]
    lapply(SM2, check.merge, groups$groups[c(8,11)]) ## ok - max all is 0.001 for ATD
    groups$groups <- groups.merge(groups$groups,groups$groups@tags[c(8,11)])
}
if(basename(args$d)=="18p-12738413-12924117") {
     s1 <- grep("rs1893217",unlist(groups$groups),value=TRUE) # C
     s2 <- grep("rs12967678",unlist(groups$groups),value=TRUE) # F
     snpin(c(s1,s2),groups$groups) # groups 3,6
     ## suggests rs12967678 should move from 6 to 3
     ## first check whether 3 and 6 should just merge - no - joint appearance in ~ 0.03 iCEL
     lapply(SM2, check.merge, groups$groups[c(3,6)])
     ## check move by making it a group on its own, and testing a merge
     groups$groups@tags[[6]] # snp to move is not the tag :)
     groups$groups@tags <- c(groups$groups@tags,s2)
     groups$groups@.Data <- c(groups$groups@.Data,list(c(s2)))
     groups$groups@.Data[[6]] <- setdiff(groups$groups@.Data[[6]], c(s2))
     lapply(SM2, check.merge, groups$groups[c(3,10)]) # looks ok, max all is ~ e-05
     lapply(SM2, check.merge, groups$groups[c(6,10)]) # looks ok
     ## what is LD?
     summary(groups$r2[groups$groups@.Data[[3]],s2]) # with new proposed group, mean=0.66, min=0.65
     summary(groups$r2[groups$groups@.Data[[6]],s2]) # with old group, mean=0.63, min=0.55
     ## conclude: support move
     gnew <- groups.merge(groups$groups,groups$groups@tags[c(3,10)])
     pattern.plot(SM2,gnew,groups$r2)
     groups$groups <- gnew
}

## merge
mfun <- function(s1,s2) {
    s1 <- grep(s1,unlist(groups$groups),value=TRUE)
    s2 <- grep(s2,unlist(groups$groups),value=TRUE)
    if(length(s1)!=1 || length(s2)!=1)
        stop("problem identifying snps")
    m <- which(apply(snpin(c(s1,s2),groups$groups),2,any))
    if(length(m)!=2)
        stop("problem with snpin for ",s1, " ",s2)
    message("r2")
    print(groups$r2[s1,s2])
    print(lapply(SM2, check.merge, groups$groups[m]))
    gnew <- groups.merge(groups$groups,groups$groups@tags[m])
    print(pattern.plot(SM2,gnew,groups$r2))
    gnew
}
    
if(basename(args$d) %in% c("2q-204446380-204816382","2q-204446380-204816382+rs117701653")) { # merge F, J
    gnew <- mfun(s1="rs112783914",s2="rs11571297")
    groups$groups <- gnew
}
if(basename(args$d)=="11q-118341921-118765600") { # merge A, B
    gnew <- mfun(s1="rs71960818",s2="rs10790268")
    groups$groups <- gnew
}
if(basename(args$d)=="22q-21712414-220220050") { # merge A, B
    gnew <- mfun(s1="rs7444",s2="rs5749495")
    groups$groups <- gnew
}
if(basename(args$d)=="3q-188069360-188139629") { # merge A, B
    gnew <- mfun(s1="rs715199",s2="rs2030519")
    groups$groups <- gnew
}
if(basename(args$d)=="7p-50900900-51133350") {
    gnew <- mfun(s1="rs7780389",s2="rs34414436")
    groups$groups <- gnew
}
if(FALSE) {
checks <- c("rs34050244","rs3087243","rs2162610","rs3116499","rs112783914","rs11571297")
c2 <- sapply(checks, function(x) grep(x, unlist(groups$groups), value=TRUE))
snpin(c2,groups$groups) # 7, 7, 9, 5, 7, 7
groups$r2[c2,c2] ## first two should be merged

min(groups$r2[ groups$groups@.Data[[7]], groups$groups@.Data[[7]] ])

mfun("rs34050244","rs3087243")
mfun("rs2162610","rs3116499") ## not high enough ld to be in same group
}


save(groups,file=gfile)

pattern.plot(SM2,groups$groups,groups$r2)
ggsave(file=file.path(args$d,paste0(prefix,"-groups.pdf")))
       
