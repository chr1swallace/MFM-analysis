#!/usr/bin/env Rscript
## follows on from GUESS-qc.R, may need some code from that inserted

library(randomFunctions) # devtools::install_github("chr1swallace/random-functions")
library(annotSnpStats)
library(magrittr)
library(RColorBrewer)
## library(reshape)
## library(plyr)
## library(igraph)
library(ggplot2)
theme_set(theme_bw())
## library(gridExtra)
library(cowplot)
setwd("~/Projects/cvs/runguess")
source("~/DIRS.txt") #CVSROOT = "/scratch/wallace/cvs"

## files
args <- getArgs(default=list(d=file.path(CVSROOT,
                                         "12q-57626582-58488667"
       ## "2q-204446380-204816382+rs117701653"#                   "13q-42834735-43101763"
   ## "5q-141414460-141641922" # times out
                                        #"19q-49094856-49278082"
    #"1q-155019710-156045183"
                                        #"17q-38718890-38878827"
                                       ## "4q-122973062-123565302"  #"14q-101290463-101328739"
## "6q-167344898-167547954"
                                       ## "16q-85988360-86024104"
    #"10p-6030243-6169685"
                                         ),keep=FALSE,cppthr=0.99))
## args <- getArgs(default=list(d=file.path(CVSROOT,"16p-11017058-11307024"),
##                              cppthr=0.999)) # 2q-204446380-204816382 #10p-6030243-6169685 1q-172650685-172940450
args$d <- sub("/GUESS$","",args$d)
args

################################################################################

## functions and set some parameters
devtools::load_all("~/RP/GUESSFM")
source("~/Projects/cvs/runguess/R/functions.R")
source("~/Projects/cvs/runguess/R/GUESS-common.R")

## start.plot <- function(filebase,height=8,width=8,type=c("pdf","png"),...) {
##   type <- match.arg(type)  
##   suffix <- switch(type, png=".png", pdf=".pdf")
##   f <- file.path(plot.dir,paste0(filebase,suffix))
##   message("printing to ",f)
##   if(type=="png")
##     png(f,height=height,width=width,units="in",res=300,...)
##   if(type=="pdf")
##     pdf(f,height=height,width=width)
##   return(f)
## }
## doplot <- function(filebase,thing,display=interactive(),...) {
##   ## if("ggplot" %in% class(thing))
##   ##   thing <- thing + ggtitle(paste(REGION,NAME))
##   if(display)
##     print(thing)
##   f <- start.plot(filebase,...)
##   print(thing)
##   dev.off()
##   return(f)
## }

## lookup <- c(IL2RA="10p-6030243-6169685",
##             CTLA4="2q-204446380-204816382",
##             IFIH1="2q-162960873-163361685",
##             IL2="4q-122973062-123565302",
##             IL7R="5p-35799296-36034866",
##             ETS1="11q-128249430-128480513",
##             CCR7="17q-38718890-38878827",
##             IRF8="16q-85988360-86024104",
##             DEXI="16p-11017058-11307024",
##             RNLS="10q-90008047-90220370",
##             SH2B3="12q-111699146-113030487" ,
##             IKZF4="12q-56369506-56482180" ,
##             PTPN2="18p-12738413-12924117" ,
##             IL10="1q-206802440-207032751" ,
##             TYK2="19p-10396336-10628468" ,
##             LPP="3q-188069360-188139629",
##             BACH2="6q-90806835-91039808",
##             IKZF1="7p-50366637-50694384"
##             )
## rev.lookup <- structure(names(lookup), names=lookup)
## REGION <- basename(args$d)
## NAME <- if(REGION %in% names(rev.lookup)) { rev.lookup[REGION] } else { "" }

## message("Running REGION: ",REGION,"\tNAME: ",NAME)

## plot.dir <- file.path(args$d,"plots")
## if(NAME!="") { plot.dir <- paste0(plot.dir,"-",NAME) }
if(!file.exists(plot.dir))
  dir.create(plot.dir,recursive=TRUE)

################################################################################

## load genotype data
(load(file.path(args$d,"GUESS/all-data.RData"))) # DATA, pcs

## load models
(load(file.path(args$d,paste0("snpmod-",sub("0.","",args$cppthr),".RData"))))
DISEASES <- names(SM2) # intersect(c("ICOELIAC","RA","T1D","GRAVES","JIA","MS", "RA"),
                      #c("ICOELIAC","RA","T1D","GRAVES","MS","JIA"),
                      # names(SM2))
SM2 <- SM2[DISEASES]
best.models(SM2)

## load groups
(load(file.path(args$d,paste0("snpmod-",sub("0.","",args$cppthr),"-groups.RData"))))

groups2 <- groups$groups

##load manhattans
library(data.table)
manhattan <- fread(file.path(args$d,"single-snp-pvalues.csv"),skip=1)
nm <- scan(file.path(args$d,"single-snp-pvalues.csv"),what="",nlines=1)
manhattan <- manhattan[,-1]
setnames(manhattan,nm)

## remove extraneous "sample" from sample ids if needed
if(all(grepl("samples",colnames(samples(DATA))))) {
  s <- samples(DATA)
  colnames(s) <- sub("samples.","",colnames(s))
  samples(DATA) <- s
}
## data.uk <- as(DATA[samples(DATA)$country=="UK" & samples(DATA)$phenotype=="CONTROL", ], "SnpMatrix")
man <- melt(manhattan,c("rs_id","position"),grep("^p\\.",colnames(manhattan),value=TRUE))
man[,trait:=sub("p.","",variable)]
table(man$trait)
DISEASES
head(man)
setnames(man,"value","p")
#man <- man[trait %in% DISEASES & !is.na(p),]

## summary of fine mapping results
summx <- guess.summ(SM2,groups=groups2,snps=DATA@snps,position="position")
summx <- scalepos(summx,position="position",prange=c(min(manhattan$position),max(manhattan$position)))
chr <- sub("[pq].*","",basename(args$d))
summx$chr <- paste0("chr",chr)
save(summx, file=file.path(args$d,paste0("summx2-",sub("0.","",args$cppthr),".RData")))

## add to man
man <- man[trait %in% unique(summx$trait) & !is.na(p),]
man <- merge(man,summx[,c("rs_id","snp","trait","tag")],all.x=TRUE,by=c("rs_id","trait"))
man <- man[c(which(is.na(man$tag)),which(!is.na(man$tag))),]
head(man)

library(ggplot2)
## gene locations
if(!file.exists(file.path(args$d,"genes.RData"))) {
    any.genes <- FALSE
} else {
    (load(file.path(args$d,"genes.RData")))
    any.genes <- length(genes)>0
}

if(any.genes) {
    ## library(ggbio)
    ## p.genes <- ggplot(genes) + geom_alignment(aes(group=genename)) #+ scale_x_sequnit("Mb")
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
## revaxis <- function(p) {
## # extract gtable
## g <- ggplot_gtable(ggplot_build(p))
##     xl <- which(g$layout$name == "axis-l")
##     xr <- which(g$layout$name == "axis-r")
##     ax <- g$grobs[[xl]]
##     ## ax <- g$grobs[[xl]]$children[[2]]
##     ## ax$widths <- rev(ax$widths)
##     ## ax$grobs <- rev(ax$grobs)
##     #ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
##     g$grobs[c(xl,xr)] <- g$grobs[c(xr,xl)]
## #    g$grobs[xr]$childrenOrder <- rev(g$grobs[xr]$childrenOrder)
##     g$grobs[xr] <- t(g$grobs[xr])
##     g$widths[c(2,3,5,6)] <- g$widths[c(6,5,3,2)]
##     ## gtable_show_layout(g)
##     ## plot(g)
## ## ia <- which(g$layout$name == "axis-l")
## ## ax <- g$grobs[[ia]]$children[[2]]
## ## ax$widths <- rev(ax$widths)
## ## ax$grobs <- rev(ax$grobs)
## ## ax$grobs[[1]]$x <- ax$grobs[[1]]$x - unit(1, "npc") + unit(0.15, "cm")
## ## pp <- c(subset(g$layout, name == "panel", select = t:r))
## ## g <- gtable_add_cols(g, g$widths[g$layout[ia, ]$l], length(g$widths) - 1)
## ## g <-  gtable_add_grob(g, ax, pp$t, length(g$widths) - 1, pp$b)
## ## g$grobs[[ia]]$children[[2]] <- NULL
## ## ##############################
## ## ia <- which(g$layout$name == "ylab-l")
## ## ylab <- g$grobs[[ia]]
## ## #g <- gtable_add_cols(g, g$widths[g$layout[ia, ]$l], length(g$widths) - 1)
## ## g <-  gtable_add_grob(g, ylab, pp$t, length(g$widths) - 2, pp$b)
## ## g$grobs[[ia]]$label = ''
##     g
## }


## X <- sm(DATA)
## colnames(X) <- addlabel(colnames(x),groups$groups)
plotter <- function(man,summ) {

    xl <- min(man$position,na.rm=TRUE)
    xm <- max(man$position,na.rm=TRUE)

    pvals <- ggplot(man,aes(x=position,y=-log10(p) ,col=tag)) +
      geom_point() +
      facet_grid(trait ~ .,scales="free_y") +
      addlines(summ[!is.na(summ$x.scale),]) +
      geom_hline(yintercept=8,linetype="dotted",col="grey") + 
      geom_hline(yintercept=0,col="grey") + 
      theme(strip.text.y = element_text(size = 9, angle = 0)) +
      scale_y_continuous(breaks=function(l) c(round(l[1]),mean(round(l)),max(round(l))),
                         labels=function(b) c("",b[-1]))
    ## pvals
    
    signals <- signal.plot(summ) + scale_y_continuous(breaks=c(0,0.5,1),labels=c("","","1"))
    chr <- ggchr(summ)
    snps <- ggsnp(summ) + addlines(summ,xvar="x.scale") #+ theme(text=element_text(size=4))
    lds <- ggld(sm(DATA), summ)
    
    nd <- length(DISEASES)
    G <- list(pvals,chr,signals,snps)
    h <- c(nd,nd/4,nd,nd/2)
    if(any.genes) {
        G <- c(list(p.genes +addlines(summ) #+ ggtitle(basename(args$d))
                    ), G)
        h <- c(nd/2,h)
    }
    if(!is.null(lds)) {
        h <- c(h,nd)
        G <- c(G,list(lds))
    }
    G %<>% lapply(., function(p) p + xlim(xl-1e+4,xm+1e+4) +
                                 theme(legend.position="none",
                                       axis.title.x = element_blank(), 
                                       axis.text.x = element_blank(), 
                                       axis.ticks.x= element_blank(), 
                                       strip.background=element_rect(fill="white"),
                                       panel.spacing = unit(0, "lines"),
                                       panel.border=element_rect(colour="grey20"),
                                       plot.margin= unit(c(0, 1, 0, 0.5), "lines")))
    plot_grid(plotlist=G,ncol=1,align="v",axis="lr",
              rel_heights=h)
    
}

summx$snp <- addlabel(summx$snp,groups$groups)
summx$tag <- addlabel(summx$tag,groups$groups)
man$tag <- addlabel(man$tag,groups$groups)
man <- man[order(man$position),]
plotter(man,summx)

doplot("summx2", plotter(man,summx), height=6,width=8)

snps.use <- subset(summx,
                   ## trait %in% DISEASES &
                   pp>0.01,
                   select="snpnum",drop=TRUE) %>%
  unique() %>% setdiff(., NA)

if(length(snps.use)) {
    summx2 <- subset(summx,snpnum %in% snps.use) # & trait %in% DISEASES)
    summx2 <- GUESSFM:::summ.setminmax(summx2)
    summx2 <- scalepos(summx2,position="position",prange=c(min(man$position),max(man$position)))
    man[ !(rs_id %in% summx2$rs_id), tag:=NA]
    man[ !(rs_id %in% summx2$rs_id), snp:=NA]
    man[,tag:=NA]
    man2 <- merge(man,as.data.table(summx2)[,.(rs_id,tag)],all.x=TRUE)
    man2 <- man2[order(man2$position),]
    doplot("summx2-001", plotter(man2,summx2), height=6)
}


## summx2 <- subset(summx2,summx2$trait %in% c("ATD","CELIAC","MS","T1D"))
## summx3 <- subset(summx3,summx3$trait %in% c("ATD","CELIAC","MS","T1D"))
## man <- subset(man,man$trait %in% c("ATD","CELIAC","MS","T1D"))
bthr <- 0.02
hsnps <- character()
while(!length(hsnps)) {
    bm <- best.snps(SM2,bthr)
    hsnps <- lapply(bm,rownames) %>% unlist() %>% unique()
    bthr <- bthr/2
}
if(length(hsnps) && hsnps!="1") {
    summx3 <- subset(summx,make.names(rs_id) %in% hsnps)
    summx3 <- GUESSFM:::summ.setminmax(summx3)
    summx3 <- scalepos(summx3,position="position",prange=c(min(man$position),max(man$position)))
    man[ !(rs_id %in% summx3$rs_id), tag:=NA]
    man[ !(rs_id %in% summx3$rs_id), snp:=NA]
    man[,tag:=NA]
    man3 <- merge(man,as.data.table(summx3)[,.(rs_id,tag)],all.x=TRUE)
    man3 <- man3[order(man3$position),]
    doplot("summx2-best", plotter(man3,summx3), height=8.5,width=12)
}

tmp <- summx2[,c("tag","snp","position","trait","pp","ppsum")]
w1 <- reshape::recast(tmp, tag + snp + position ~ trait,  measure.var="pp")
w2 <- reshape::recast(tmp, tag + snp + position ~ trait,  measure.var="ppsum")
colnames(w2)[-c(1:3)] <- paste0(colnames(w2)[-c(1:3)],".sum")
w <- merge(w1,w2)
max.sum <- apply(w[,grep(".sum$",colnames(w)),drop=FALSE],1,max)
w <- w[order(w$tag,w$position),]
write.table(w, file=file.path(plot.dir,"summw2.csv"))

