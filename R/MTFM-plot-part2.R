#!/usr/bin/env Rscript
library(methods)
library(randomFunctions)
library(annotSnpStats)
library(magrittr)
source("~/DIRS.txt")
setwd("~/Projects/cd4chic/GUESSFM")
## library(RColorBrewer)
## colsd <- brewer.pal(8,"Dark2")
## colsp <- brewer.pal(8,"Paired")
## cols <- c(colsp[1:4], colsd[c(2,3,4,8)])
## names(cols) <- c("CEL","iCEL","RA","iRA","ATD","JIA","MS","T1D")
## myord <- c("ATD","CEL","iCEL","JIA","MS","RA","iRA","T1D")
## cols <- cols[myord]
## dput(cols)
                                        #library(GUESSFM)
library(devtools)
load_all("~/RP/GUESSFM")
load_all("~/RP/MTFM")
library(ggplot2)
source("~/Projects/cvs/runguess/R/plot2-functions.R")
devtools::load_all("~/RP/snpHaps")

args <- getArgs(default=list(d=file.path(CVSROOT,
                                         ## "1p-113863087-114527968"
                                         # "18p-12738413-12924117" # PTPN2
                                           ## "10p-6030000-6220000"
                                         ## "1p-113863087-114527968" # empty groupld
                                         ## "9p-34649442-34974974"
                                         ## "13q-42834735-43101763"
                                         ## "1p-113863087-114527968"
                                         ## "2p-60914729-61892409"
  ## "10p-6030000-6220000"
                                         ## "3q-119104981-119300655"
                                         #"2q-100544954-101038647"
  ## "1p-117030982-117142239"
                                         ## "1p-113863087-114527968" # One SNP (same) across multiple traits
                                        "2q-204446380-204816382"
                                         ##  "4q-122973062-123565302"
                                         ## "10p-6030243-6169685"
                                         ),cppthr=0.99,keep=FALSE))
args$d <- sub("/GUESS$","",args$d)
args
plot.dir <- file.path(args$d,"plots")

## use newer ctla4
## if(basename(args$d)=="2q-204446380-204816382+rs117701653")
##     plot.dir  %<>%  sub("+rs117701653","",.)
## if(basename(args$d)=="2q-204446380-204816382")
##     q("not processing old 2q")

ord <- function(x) {
    if(is.factor(x))
        return( factor(as.character(x),levels=intersect(myord, as.character(x))) )
    if(is.list(x)) {
        names(x) <- sub(".mid","",names(x))
        neword <- intersect(myord,names(x))
        return(x[neword])
    }
    warning("x not a class to be ordered")
    x
}

## load genotype data
(load(file.path(args$d,"GUESS/all-data.RData"))) # DATA, pcs

##load manhattans
library(data.table)
manhattan <- fread(file.path(args$d,"single-snp-pvalues.csv"),skip=1)
nm <- scan(file.path(args$d,"single-snp-pvalues.csv"),what="",nlines=1)
manhattan <- manhattan[,-1]
setnames(manhattan,nm)

## do diseases with minp < 1e-6
(load(file.path(args$d,"conditional.RData")))
conditional <- as.data.table(conditional)
cond.orig=conditional
## snp            p trait step
##  1:    rs61839660.6094697.C.T 3.599946e-34   T1D    1
##  2:    rs11594656.6122009.T.A 5.853643e-12   T1D    2
##  3:    rs12220852.6133371.T.C 8.788855e-10   T1D    3
##  4:    rs41295159.6148535.C.G 4.240106e-06   T1D    4
##  5:     rs7893467.6079035.G.T 3.759859e-05   T1D    5
##  6:    rs11256899.6168832.T.C 2.176727e-04   T1D    6
##  7:    rs12722491.6101430.G.T 1.954979e-04   T1D    7
##  8:    rs41295049.6109676.G.A 4.089382e-04   T1D    8
##  9:      rs706778.6098949.C.T 4.685687e-05    RA    1
## 10:     rs2104286.6099045.T.C 1.126146e-13    MS    1
## 11:    rs12722502.6093139.C.T 6.094697e-04    MS    2
## 12:    rs12722496.6096667.A.G 6.126800e-06   JIA    1
## 13:    rs62626308.6040028.C.A 9.037799e-04   JIA    2
## 14:     rs3134882.6102725.G.A 5.561912e-04   JIA    3
## 15:      rs706779.6098824.T.C 4.629063e-08   ATD    1
## 16:   rs41295091.6120051.TG.T 4.994410e-04   ATD    2
## 17:     rs4747846.6074451.G.C 2.964853e-06   CEL    1
## 18: rs78296807.6131411.T.TTCC 3.172801e-04   CEL    2
## 19:    rs35285258.6118770.C.T 3.964439e-05  iCEL    1
## 20:    rs12264061.6073487.A.G 5.447823e-04  iCEL    2
## 21:      rs706778.6098949.C.T 7.554578e-08   iRA    1
conditional <- conditional[p<1e-6,]
conditional[,i:=1:.N,by="trait"] ## avoid occasionally conditionally on sub-sig to get a sig next step
conditional <- conditional[i==step, names(conditional) , by="trait", with=FALSE]
conditional$trait %<>% ord()
conditional
source("~/Projects/cvs/runguess/R/functions.R")

##snpmod
(load(file.path(args$d,paste0("snpmod-",sub("0.","",args$cppthr),".RData"))))
SM2 <- ord(SM2[ intersect(conditional$trait,names(SM2)) ])
conditional <- conditional[ trait %in% names(SM2), ] ## SM2 output has been cleaned to remove a couple of regions that are only significant because of a rogue SNP

## summx
(load(file.path(args$d,paste("summx2-99.RData"))))
summx$trait  %<>% ord()

## groups
(load(file.path(args$d,paste("snpmod-99-groups.RData"))))

## rename groups?
library(yaml)
grename <- read_yaml("~/Projects/cvs/runguess/groups.yaml")
r <- basename(args$d)
remap <- grename[[r]]
    
## MTFM
int <- file.exists(file.path(args$d,"iMTFM.RData"))
if(int) {
    (load(file.path(args$d,paste("iMTFM.RData"))))
    impp <- ord(mpp); igpp <- ord(gpp)
}
uk <- file.exists(file.path(args$d,"MTFM.RData"))
if(uk){
    (load(file.path(args$d,paste("MTFM.RData"))))
    mpp %<>% ord()
    gpp %<>% ord()
}

bm <- best.models(SM2,cpp.thr=0.99)
bm <- lapply(bm,function(x) {x$gr <- mod2group(x$str,groups=groups$groups); x})
bm <- lapply(bm, function(x) { 
    tapply(1:nrow(x),x$gr,function(i) sum(x$PP[i]))
})
conditional$gr <- mod2group(as.character(conditional$snp),groups=groups$groups)
cond.orig$gr <- mod2group(as.character(cond.orig$snp),groups=groups$groups)
cond.orig$gr <- regroup(toupper(letters[ as.numeric(cond.orig$gr) ]))
cond.orig

message("individual")
print(lapply(bm, function(x) x[x>0.1]))
message("MTFM")

## ## for now, only present mid results, the others can be used as bounds in tables
if(uk) {
##     if(any(grepl("mid",names(gpp)))) {
##     ## gpp <- gpp[ grep(".mid",names(mpp)) ]
##     ## names(gpp) %<>% sub(".mid","",.)
##     ## mpp <- mpp[ grep(".mid",names(mpp)) ]
##     ## names(mpp) %<>% sub(".mid","",.)
##     if(int) {
##         igpp <- igpp[ grep(".mid",names(impp)) ]
##         names(igpp) %<>% sub(".mid","",.)
##         impp <- impp[ grep(".mid",names(impp)) ]
##         names(impp) %<>% sub(".mid","",.)
##     }
## }
print(lapply(gpp, function(x) x[x>0.1]))
for(i in seq_along(mpp)) {
    mpp[[i]]$var <- rownames(mpp[[i]])
}
}
if(int) {
    print(lapply(igpp, function(x) x[x>0.1]))
    for(i in seq_along(impp))
        impp[[i]]$var <- rownames(impp[[i]])
}

################################################################################
## tables

library(Hmisc)

conditional$gr <- regroup(toupper(letters[ as.numeric(conditional$gr) ]))
## conditional$trait %<>% as.character()
conditional <- conditional[order(conditional$trait,conditional$step),]

f <- function(n) {
    n2 <- n
    nulls <- which(n=="-")
    if(length(nulls))
        n2[nulls] <- "--"
    nots <- which(n!="--")
    if(length(nots)) {
        n2[nots] <- strsplit(n,"-")  %>%  lapply(., function(x) {
            ifelse(x=="0", "0", regroup(toupper(letters[as.numeric(x)])))
        } %>%  sort()  %>%  paste(.,collapse="+" ))  %>% unlist()
    }
    if(any(n2==""))
        n2[n2==""] <- "NULL"
    n2
}

ff <- function(x) {
    y <- x[x>0.1]
    names(y) <- f(names(y))
    y[order(y,decreasing=TRUE)]
}

fff <- function(x,...) {
    v <- lapply(x,ff)  %>% unlist()
    x <- as.data.frame(v)
    colnames(x) <- "PP"
    x$trait <- sub("\\..*","",names(v)) 
    x$model <- sub(".*\\.","",names(v))
    split(x,x$trait)
    ## latex(x,rowname=NULL,...)
}
## kk <- latex(conditional,caption="stepwise",file=file.path(plot.dir,"tables.tex"))
## kk <- fff(bm,caption="GUESSFM",file=file.path(plot.dir,"tables.tex"),append=TRUE)
## if(uk)
##     kk <- fff(gpp,caption="UK MTFM",file=file.path(plot.dir,"tables.tex"),append=TRUE)
## if(int)
##     kk <- fff(igpp,caption="International MTFM",file=file.path(plot.dir,"tables.tex"),append=TRUE)

conditional
lc <- split(conditional,conditional$trait)
lg <- fff(bm)
if(uk)
    luk <- fff(gpp)
if(int)
    lint <- fff(igpp)

library(purrr)
ph <- "ATD"
tabph <- function(ph) {
    blank2 <- data.frame(model=NA,PP=NA) #("",nrow=1,ncol=2,dimnames=list(ph,c("model","PP")))
    ##rownames(blank2) <- ""
    l <- list(as.data.frame(lc[[ph]][,.(snp=paste(gr, sub("\\..*","",snp), sep="/"),p)]),
              lg[[ph]][,c("model","PP")])
    if(uk) {
        tmp <- if(ph %in% names(luk)) { luk[[ph]][,c("model","PP")] } else { blank2 }
        l <- c(l, list(tmp))
    }
    if(int) {
        tmp <- if(ph %in% names(lint)) { lint[[ph]][,c("model","PP")] } else {blank2}
        l <- c(l, list(tmp))
    }
    ret <- purrr::reduce(l, cbindpad)
    names(ret) <- sub("\\.1","",names(ret))
    ret
}

allph <- intersect(names(lc),names(lg))
tab <- lapply(allph, tabph) #%>% reduce(., rbind)
rtab <- purrr::reduce(tab,rbind)
nrep <- 1 + uk + int

## html output
library(rmarkdown)
library(DT)
hh <- rtab
hh <- cbind(disease=rep(allph,times=sapply(tab,nrow)),hh)
rownames(hh) <- NULL
for(j in which(colnames(hh) %in% c("p","PP")))
    hh[[j]] %<>% signif(.,3)

                                        # a custom table container
sketch = htmltools::withTags(table(
  class = 'display',
  thead(
    tr(
      th(rowspan = 2, 'Disease'),
      th(colspan = 2, 'Stepwise'),
      th(colspan = 2, 'Stochastic'),
      if(uk) {th(colspan = 2, 'Joint (uk)')} else {NULL},
      if(int) {th(colspan = 2, 'Joint (int)')} else {NULL}
    ),
    tr(
      lapply(c("SNP","P",rep(c('Mod', 'PP'), nrep)), th)
    )
  )
))
print(sketch)

table.html <- datatable(hh, container = sketch, rownames = FALSE,
                        options = list(pageLength = nrow(hh), autoWidth = TRUE),
                        caption="Stepwise and stochastic search results. Stepwise are shown to p>10-4, to allow more complete comparison with stochastic models, although only results with p<10-6 are considered 'significant' in our formal comparisons.")

## latex output
rtab$p <- format(rtab$p,digits=2)  %>%  sub("NA","",.)  %>%
  sub("([0-9.]+)e-([0-9]+)","$\\1 \\\\times 10^{-\\2}$",.)
rtab$p
cgroups <- c("Stepwise","Stochastic",
             if(uk) { "Joint (uk)" } else { NULL },
             if(int) { "Joint (int)" } else { NULL })
    kk <- latex(rtab,
          rowlabel="",
          cgroup=cgroups,n.cgroup=c(2,rep(2,nrep)),
          colheads=c("SNP","p",rep(c("Mod","PP"),nrep)),
          rgroup=allph,n.rgroup=sapply(tab,nrow),
          ## rgroupTexCmd="rgroup",
          rowname=rep("",nrow(rtab)),
          dcolumn=TRUE,dec=3,
          booktabs=TRUE,
          table.env=FALSE,
          file=file.path(plot.dir,"tables.tex"))
    

################################################################################

## group LD plot
gr <- c(conditional$gr, rtab[,grep("model",colnames(rtab))]) %>%  unlist()  %>%  strsplit(.,"+")  %>% unlist()  %>%  unique()  %>%  setdiff(.,c("-","+",NA))  %>% sort() # might be remapped
gr0 <- setdiff(gr,"0")
if(length(gr0)>1) {
    ngr <- sapply(unregroup(gr0), function(x) which(toupper(letters)==x))
    tgr <- snps(groups$groups)[ngr]
    r2 <- matrix(0,length(tgr),length(tgr),dimnames=list(gr0,gr0))
    for(i in seq_along(tgr)) {
        for(j in seq_along(tgr)) {
            r2[i,j] <- mean(groups$r2[ tgr[[i]], tgr[[j]] ])
        }
    }
    p <- (ldplot(r2))
    if(!is.null(p)) {
        pdf(file.path(plot.dir,"groupld.pdf"),height=3,width=6)
        print(p)
        dev.off()
    }
} else {
    p <- NULL
}

## final html
tmp <- tempfile(fileext=".Rmd")

if(!is.null(p)) {
paste(
"---
title: ",r," 
---
## Overall summary of results (top) and average LD ($r^2$) between SNP groups (bottom)
```{r,echo=FALSE}
table.html
print(p)
```
")  %>% cat(., file=tmp)
} else {
paste(
"---
title: ",r," 
---
## Overall summary of results
```{r,echo=FALSE}
table.html
```
")  %>% cat(., file=tmp)
}

render(tmp,
       output_format="html_document",
       output_file=file.path("/home/cew54/scratch/MFM-output",r,"results-table.html"))

################################################################################
## other plots
library(cowplot)
library(ggplot2)

mpp.marg <- function(x,PP="shared.pp"){
    marg.snps.vecs(rownames(x),x[[PP]])
}
if(uk)
    mmpp <- lapply(mpp,mpp.marg) 
if(int)
    immpp <- lapply(impp,mpp.marg) 

chr <- sub("[pq].*","",basename(args$d))
## gene locations
if(!file.exists(file.path(args$d,"genes.RData"))) {
    any.genes <- FALSE
} else {
    (load(file.path(args$d,"genes.RData")))
    any.genes <- length(genes)>0
    df <- as.data.table(as.data.frame(genes))
    df2 <- df[,.(start=min(start),end=max(end)),by="genename"]
    p.genes <- ggplot(df) + ggplot2::geom_rect(aes(xmin=start,xmax=end,ymin=0,ymax=1)) +
      geom_segment(aes(x=start,xend=end,y=0.5,yend=0.5),data=df2) +
      facet_grid(genename ~ .) + theme(axis.title.y = element_blank(), 
                                       axis.text.y = element_blank(), 
                                       axis.ticks.y= element_blank(),
                                       axis.line=element_blank(),
                                       strip.text.y = element_text(size=8,angle=0),
                                       strip.text=element_text(face="italic"),
                                       panel.spacing = unit(0, "lines"),
                                       panel.grid=element_blank())
}

pdf(file.path(args$d,"plots","summx-ind.pdf"),height=7,width=11)
s1 <- bigplotter(SM2)
dev.off()
if(uk) {
    pdf(file.path(args$d,"plots","summx-mpp.pdf"),height=7,width=11)
    s2 <- bigplotter(mmpp)
    dev.off()
}
if(int) {
    pdf(file.path(args$d,"plots","summx-impp.pdf"),height=7,width=11)
    s3 <- bigplotter(immpp)
    dev.off()
}

## snps for haplotypes
## first pick groups
if(int) {
    g <- lapply(c(bm,gpp,igpp), function(x) names(x)[x>0.05])  %>% unlist()  %>%  unique()  %>%
      strsplit(.,"-")  %>%  unlist()  %>% unique()
} else if(uk) {
    g <- lapply(c(bm,gpp), function(x) names(x)[x>0.05])  %>% unlist()  %>%  unique()  %>%
      strsplit(.,"-")  %>%  unlist()  %>% unique()
} else {
    g <- lapply(bm, function(x) names(x)[x>0.05])  %>% unlist()  %>%  unique()  %>%
      strsplit(.,"-")  %>%  unlist()  %>% unique()
} 
g <- unique(c(g,which(toupper(letters) %in% unregroup(conditional$gr))))

## then a good tag for each group
hsnps <- lapply(SM2, function(x) {
    tmp <- as.data.table(x@snps)
    tmp$gr <- (mod2group(tmp$var,groups=groups$groups))
    tmp <- tmp[gr %in% g & Marg_Prob_Incl>0.01,]
    tmp <- tmp[order(gr,Marg_Prob_Incl,decreasing=TRUE),]
    tmp[,.SD[1,.(var)],by="gr"]
})  %>%  rbindlist()

toadd <- NULL
if(uk)
    toadd <- mmpp
if(int)
    toadd <- c(toadd,immpp)
if(uk|int) {
    multi <- lapply(toadd, function(x) {
        tmp <- as.data.table(x)
        tmp$gr <- mod2group(tmp$var,groups=groups$groups)
        tmp <- tmp[gr %in% g & Marg_Prob_Incl>0.01,]
        tmp <- tmp[order(gr,Marg_Prob_Incl,decreasing=TRUE),]
        tmp[,.SD[1,.(var)],by="gr"]
    })  %>%  rbindlist()
    hsnps <- unique(rbind(hsnps,multi))
}
hsnps <- hsnps[var!="1",]
hsnps

## customize IL2RA
if(args$d=="/rds/project/cew54/rds-cew54-wallace-share/Projects/cvs/10p-6030000-6220000") {
    ## In the extended data please replace the D SNP with 
    ## rs41295055 (as in 6b).
    hsnps[var=="rs7089861.6110326.C.G", var:="rs41295055.6111622.C.T"]
}
## customize ctla4
if(args$d=="/rds/project/cew54/rds-cew54-wallace-share/Projects/cvs/2q-204446380-204816382") {
## replace rs231779 with rs231775, which will make it consistent with figure 5d)
     hsnps[var=="rs231779.204734487.C.T", var:="rs231775.204732714.A.G"]
##for rs78847176/J and rs146377774/A?  These two SNP groups are not noted in the table on page 59 and there is no minor allele shown.  Are they very rare variants with an effect on one of the diseases?
    ## hsnps <- hsnps[!(var %in% c("rs78847176.204725944.C.T","rs146377774.204450286.A.G")),]
    hsnps <- hsnps[!(var %in% c("rs146377774.204450286.A.G")),]
}
    
## new - 15/10/18 - response to Linda comment re PTPN2
## ensure conditionals are included
csnps <- data.table(gr=mod2group(as.character(conditional$snp),groups=groups$groups),
                    var=as.character(conditional$snp))
hsnps <- unique(rbind(hsnps,csnps))


## order
snps <- DATA@snps
  tsplit <- split(hsnps$var, hsnps$gr) #[ cl ]
  tsplit <- lapply(tsplit, function(x) x[ order(snps[x, "position"]) ])
  ## order tag groups by leftmost snp
  first <- sapply(tsplit,"[",1)
  tsplit <- tsplit[ order(snps[first, "position"]) ]
hsnps <- data.table(gr=rep(names(tsplit),sapply(tsplit,length)),
                    var=unlist(tsplit))

                                        #o <- order( DATA@snps[hsnps$var,]$position )
## o <- order(hsnps$gr,DATA@snps[hsnps$var,]$position)
## hsnps <- hsnps[o,]


## allele frequencies
df <- DATA@samples
phenotypes <- names(SM2)
ph <- unique(df[,c("phenotype","country")])
ph <- subset(ph,!is.na(phenotype))

ph <- ph[order(ph$phenotype,ph$country),]
maf <- matrix(0,nrow(ph),nrow(hsnps), dimnames=list(apply(ph,1,paste,collapse="/"),
                                          hsnps$var))
for(i in 1:nrow(ph)) {
  maf[i,] <- col.summary(DATA[df$country==ph[i,"country"] & df$phenotype==ph[i,"phenotype"], hsnps$var])[,"RAF"]
}
write.csv(maf,file=file.path(plot.dir,"allele-frequencies.csv"))

## add PCs to samples data.frame
df$y <- df$affected - 1

#' generate haplotypes
X <- as(DATA[,hsnps$var],"SnpMatrix")
snphap.dir <- file.path(args$d,"snphap")
if(!file.exists(snphap.dir))
    dir.create(snphap.dir,recursive = TRUE)
  d <- genhaps(X,slist=hsnps$var,
               snps=DATA@snps,
               samples=df,
               A1="A1",A2="A2",
               cols.copy=c("phenotype","PC1","PC2","PC3","PC4","country"),
               f.in=file.path(snphap.dir,"haps.in"),
               redo=TRUE,
               by=df$country)
  d$cc <- as.numeric(d$phenotype!="CONTROL")

#' haplotype analysis
base <- NULL
if(args$d=="/rds/project/cew54/rds-cew54-wallace-share/Projects/cvs/10p-6030000-6220000") 
    base <- "ATCAATTTATCTTCCC"

RESULTS <- vector("list",length(phenotypes))
names(RESULTS) <- phenotypes
for(ph in phenotypes) {
    uph <- sub("i","",ph)
    if(grepl("i",ph)) {
        countries <- samples(DATA)$country[ samples(DATA)$phenotype==uph ]  %>%  unique()
        use <- with(samples(DATA), which(country %in% countries & phenotype %in% c("CONTROL",uph)))
        covars <- c("PC1","PC2","PC3","PC4","country")
        RESULTS[[ph]] <- model.mi(dir=snphap.dir,df[use,], haps.pattern="haps.out2",covars=covars,by=df[use,]$country,
                                  base=base)
        RESULTS[[ph]]$disease <- ph
    } else {
        use <- with(df, phenotype %in% c("CONTROL",ph) & country=="UK")
        covars <- c("PC1","PC2","PC3","PC4")
        RESULTS[[ph]] <- model.mi(dir=snphap.dir,df[use,], haps.pattern="haps.UK.out2",covars=covars,
                                  base=base)
        RESULTS[[ph]]$disease <- ph
    }
}


results <- lapply(RESULTS, function(x) {
  x <- x[grep("hap",rownames(x)),]
  x$hap <- factor(rownames(x),levels=rownames(x))
  return(x)
}) %>% do.call("rbind",.)
results$disease <- factor(results$disease)
write.csv(results,file=file.path(plot.dir,"haplotype-results.csv"),row.names=FALSE)

if(args$d=="/rds/project/cew54/rds-cew54-wallace-share/Projects/cvs/10p-6030000-6220000") {
    ## reset iRA fq column which is used only for subsetting so it matches uk and the same minor haplotypes are shown for all diseases
    ra="hapACCATACCCTCCTCCC"
    notra="hapATCATACCCTCTTCTG"
    results[ results$hap %in% c(ra,notra),]
    results$Fq[ results$hap==notra & results$disease=="iRA" ] <- 0.92
    results$Fq[ results$hap==ra & results$disease=="iRA" ] <- 0.84
}

#
#' haplotype frequencies
hfreq <- hapfreq(dir=snphap.dir,df,haps.pattern="haps.out2",covars=c("country","phenotype"))
hfreq <- hfreq / matrix(rowSums(hfreq),nrow=nrow(hfreq),ncol=ncol(hfreq))
hfreq <- hfreq[!grepl("/NA",rownames(hfreq)), ]
write.csv(hfreq,file=file.path(plot.dir,"haplotype-frequencies.csv"))

## thr st every snp selected appears variable in plot
f=hfreq["UK/CONTROL",,drop=FALSE]
ss <- setdiff(colnames(f),"OTHER")  %>% strsplit(.,"") %>% do.call("rbind",.)
tss <- t(ss)
m <- tss==tss[,1]
m <- t(m)
allsame <- apply(m,2,all)
if(any(allsame)) {
    warning(basename(args$d),": hsnp with no variability: ")
    print(hsnps[allsame,])
}

#' plot results
plots <- snpHaps::plotter(results,
                          hsnps=addlabel(hsnps$var,groups$groups,grename[[r]]),
                          hfreq=hfreq["UK/CONTROL",,drop=FALSE],
                          thr=if(r=="10p-6030000-6220000") {
                                  0.009
                              } else if(r== "2q-204446380-204816382") {
                                  0.005
                              } else 0.01)

library(cowplot)
library(magrittr)
plist <- list(plots$odds + scale_colour_manual(values=cols),
              plots$freq + geom_hline(yintercept=1,width=0.01,linetype="dotted") + theme(legend.position = "none") + scale_colour_grey(),
              plots$haps)   %>% lapply(., function(p)
    p + xlim(0.75,nrow(plots$hlab)+1) +
    theme(axis.text.x=element_blank(),axis.title.x=element_blank()))
p <- plot_grid(plotlist=plist,
          align="v",ncol=1,rel_heights=c(3,1,3))
p

pdf(file.path(plot.dir,"haplotypes.pdf"),height=10,width=8)
print(p)
dev.off()

