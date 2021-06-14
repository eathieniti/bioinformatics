##### RNA-seq analysis is easy as 1-2-3 with limma, Glimma and edgeR

#' This is a tutorial to perform DE for rna seq data 
library(limma)
library(Glimma)
library(edgeR)
library(Mus.musculus)


url <- "https://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE63310&format=file"
utils::download.file(url, destfile="GSE63310_RAW.tar", mode="wb") 
utils::untar("GSE63310_RAW.tar", exdir = ".")
files <- c("GSM1545535_10_6_5_11.txt", "GSM1545536_9_6_5_11.txt", "GSM1545538_purep53.txt",
           "GSM1545539_JMS8-2.txt", "GSM1545540_JMS8-3.txt", "GSM1545541_JMS8-4.txt",
           "GSM1545542_JMS8-5.txt", "GSM1545544_JMS9-P7c.txt", "GSM1545545_JMS9-P8c.txt")

#9 FILES READ IN WITH RAW GENE LEVEL COUNTS 

for(i in paste(files, ".gz", sep=""))
  R.utils::gunzip(i, overwrite=TRUE)

x <- readDGE(files, columns=c(1,3))
class(x)

x <- readDGE(files, columns=c(1,3))


#### oRGANIZE sample information 
# add information that might have an efffects on expression levels 
#eg. basal, LP, ML. genotupe, wt vs knockout phenotype disease status, sex, age, 
# sample treatment, and batch infom

#data contains samples df that stopes cell type and batch info 
samplenames <- substring(colnames(x), 12, nchar(colnames(x)))
colnames(x) <- samplenames
group <- as.factor(c("LP", "ML", "Basal", "Basal", "ML", "LP", 
                     "Basal", "ML", "LP"))
x$samples$group <- group
lane <- as.factor(rep(c("L004","L006","L008"), c(3,4,2)))
x$samples$lane <- lane
x$samples

##### Organize gene annotations 
#

eneid <- rownames(x)
genes <- select(Mus.musculus, keys=geneid, columns=c("SYMBOL", "TXCHROM"), 
                keytype="ENTREZID")
head(genes)

genes <- genes[!duplicated(genes$ENTREZID),]

x$genes <- genes

##### Data preprocessing 

#### Transform from the raw cae. Conver raw counts to CPM and log-CPM using the cpm function 
# 

cpm <- cpm(x)
lcpm <- cpm(x, log=TRUE)


L <- mean(x$samples$lib.size) * 1e-6
M <- median(x$samples$lib.size) * 1e-6
c(L, M)


##### Remove genes that are lowly expressed 

# genes that are lowly expressed in all samples  - in both normal and disease!!

keep.exprs <- filterByExpr(x, group=group)
x <- x[keep.exprs,, keep.lib.sizes=FALSE]
dim(x)



library('GiNA')
lcpm <- cpm(x, log=TRUE)
par(mfrow=c(1,2))
col.group <- group
levels(col.group) <-  brewer.pal(nlevels(col.group), "Set1")
col.group <- as.character(col.group)
col.lane <- lane
levels(col.lane) <-  brewer.pal(nlevels(col.lane), "Set2")
col.lane <- as.character(col.lane)
plotMDS(lcpm, labels=group, col=col.group)
title(main="A. Sample groups")
plotMDS(lcpm, labels=lane, col=col.lane, dim=c(3,4))
title(main="B. Sequencing lanes")



##### Normalize gene expression distributions 

