library(magrittr)
library(DT)
source("~/DIRS.txt")
regions <- scan("/home/cew54/scratch/cvs-table/regions.txt",what="")  %>% setdiff(.,c(#"10p-6388071-6545104",
                                                                                        "10p-6030243-6169685","2q-204446380-204816382+rs117701653"))
d <- "/home/cew54/scratch/MFM-output"
setwd(d)

## haplotypes
comm <- paste0("ls */haplotypes.pdf")
ofiles <- system(comm,intern=TRUE)
haplinks <- paste0("<a href=\"",regions,"/haplotypes.pdf\">link</a>")

## standard links
makelinks <- function(fname) {
    ofiles <- system(paste0("ls */",fname),intern=TRUE)
summxfiles <- paste0(regions,"/",fname)
summxlinks <- ifelse(summxfiles %in% ofiles,
                     paste0("<a href=\"",regions,"/",fname,"\">link</a>"),
                     "")
}

## genes
gfiles <- system(paste0("ls ",CVSROOT,"/*/genes.RData"),intern=TRUE)
genenames <- sapply(regions, function(r) {
    f <- paste0(CVSROOT,"/",r,"/genes.RData")
    if(!(f %in% gfiles))
        return("")
    (load(f))
    unique( genes@elementMetadata$genename)  %>% sort()  %>% paste(collapse=", ")
})
names(genenames) <- NULL

## diseases
dfiles <- system(paste0("ls ",CVSROOT,"/*/plots/haplotype-results.csv"),intern=TRUE)
library(data.table)
ddata <- lapply(dfiles,function(f) {
    tmp <- fread(f)
    tmp$region <- f  %>% dirname()  %>% dirname()  %>% basename()
    tmp })  %>% rbindlist()  %>% unique(.,by=c("disease","region"))
diseases <- split(ddata$disease,ddata$region)  %>% sapply(., paste,collapse=", ")
diseases <- ifelse(is.na(diseases[regions]),"",diseases[regions])

## classes
library(yaml)
cl <- read_yaml(file.path(d,"classes.yaml"))
nullstr <- function(x) 
    ifelse(is.null(x),"",paste(x,collapse="/"))
clmat <- sapply(regions, function(r) {
    sapply(names(cl), function(nm) nullstr(names(cl[[nm]][[r]])))
}) %>% t()
colnames(clmat) <- paste0("Flag.",1:3)
clmat <- as.data.frame(clmat)

fsub <- function(str) { str }
##     sub("2q-204446380-204816382","2q-204446380-204816382+rs117701653",str)
## }

m <- cbind(data.frame(region=regions,
                genes=genenames,
                diseases=diseases[regions],
                results=makelinks("results-table.html"),
                ncvs=makelinks("nsnp.pdf"),
                haplotypes=fsub(haplinks),
                GUESSFM=makelinks("summx-ind.pdf"),
                MFM.UK=makelinks("summx-mpp.pdf"),
                MFM.int=makelinks("summx-impp.pdf"),stringsAsFactors=FALSE),
           clmat)
rownames(m) <- NULL
save(m,file=file.path(d,"m.RData"))

library(rmarkdown)
render(file.path(d,"index.Rmd"))

