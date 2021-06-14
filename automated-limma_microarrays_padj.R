### Limma and curation of resulting Gene Lists ###
#####
#
# Input Files
# a. Gene Expression Data, tab-separated (colnames, rownames and data)
# b. Classlabels, tab separated (0 or 1, last line empty)
# c. Platform, tab separated (2 columns; probe ids and gene symbols)
#
# Pre checking - do manually:
# a. Set respective working directory (setwd())
# b. boxplot and decide on setting toNormalize = TRUE or FALSE
#
### ~ ###

#####
### FUNCTIONS ###
#####

# reading input files into matrices
parseFile <- function(filename, h){
  matrix <- as.matrix(read.delim(filename, header = h))
  return(matrix)
}

checkLog <- function(intensities){
  toLog <- FALSE
  if (max(intensities) > 25) toLog <- TRUE
  return(toLog)
}

# executed if negative values exist in intensities
cureNegative <- function(intensities){
  min <- abs(min(intensities)) + 0.001 # + 0.001 to avoid 0s
  intensities <- intensities + min
  return(intensities)
}

# handle data based on toLog, toNormalize or both, or un-log, normalize and re-log
handleData <- function(intensities, norm, log){
  if ( norm && log ){
    cat(sprintf("Normalizing and doing a log2 on the data.\n"))
    intensities <- normalizeQuantiles(intensities)
    intensities <- log2(intensities)
  } else if ( !norm && log ){
    cat(sprintf("Doing a log2 only on the data.\n"))
    intensities <- log2(intensities)
  } else if ( norm && !log ){ # need to un-log, normalize and then re-log
    cat(sprintf("Doing an Un-log, then normalizing and finally re-doing a log2 on the data.\n"))
    intensities <- 2 ^ intensities
    intensities <- normalizeQuantiles(intensities)
    intensities <- log2(intensities)
  }
  return(intensities)
}

runLimma <- function(intensities, classes){
  # transposing classes, as need for limma's model design
  classes <- t(classes)
  colnames(classes) <- "Group"
  design <- cbind(Intercept = 1, classes)
  fit <- lmFit(intensities, design)
  # Moderated t-statistic
  fit <- eBayes(fit)
  topTable <- topTable(fit, coef = 2, n = nrow(intensities))
  limma <- cbind(rownames(topTable), topTable$P.Value, topTable$logFC)
  colnames(limma) <- cbind("probeIDs","p-value","logFC")
  cat(sprintf("Limma finished.\n"))
  return(limma)
}

# Run Limma including p-adjusted  in results
runLimma_padj <- function(intensities, classes){
  # transposing classes, as need for limma's model design
  classes <- t(classes)
  colnames(classes) <- "Group"
  design <- cbind(Intercept = 1, classes)
  fit <- lmFit(intensities, design)
  # Moderated t-statistic
  fit <- eBayes(fit)
  topTable <- topTable(fit, coef = 2, n = nrow(intensities))
  limma <- cbind(rownames(topTable), topTable$P.Value, topTable$logFC,topTable$adj.P.Val)
  colnames(limma) <- cbind("probeIDs","p-value","logFC","adj.p-value")
  cat(sprintf("Limma (adj-p included) finished.\n"))
  return(limma)
}

pvalueClear <- function(l, cut){
  # data type change
  l <- data.frame(l)
  l[, "probeIDs"] <- as.matrix(as.character(l[, "probeIDs"]))
  l[, "p.value"] <- as.matrix(as.double(as.character(l[, "p.value"])))
  l[, "logFC"] <- as.matrix(as.double(as.character(l[, "logFC"])))
  # actual cutoff
  l <- l[l[, "p.value"] < cut, ] # delete where pvalue > 0.05
  cat(sprintf("pvalue entries > %f removed.\n", cut))
  return(l)
}

pvalueClear_padj <- function(l, cut){
  # data type change
  l <- data.frame(l)
  l[, "probeIDs"] <- as.matrix(as.character(l[, "probeIDs"]))
  l[, "p.value"] <- as.matrix(as.double(as.character(l[, "p.value"])))
  l[, "logFC"] <- as.matrix(as.double(as.character(l[, "logFC"])))
  l[, "adj.p.value"] <- as.matrix(as.double(as.character(l[, "adj.p.value"])))
  # actual cutoff
  l <- l[l[, "adj.p.value"] < cut, ] # delete where pvalue adj > 0.05
  cat(sprintf("pvalue  adj entries > %f removed.\n", cut))
  return(l)
}

# if remaining Gene Symbols found both up- and down- regulated remove from list
removeDualSignProbes <- function(l){
  signMatrix <- matrix("" , nrow = 0, ncol = 2) # gene, FC value
  positions <- vector() # initialize empty vector of row positions which will be finally removed
  for (i in 1:nrow(l)){
    position <- match(l[i, "GeneSymbol"], signMatrix[, 1])
    if (is.na(position)){ # if not checked yet add in signMatrix
      signMatrix <- rbind(signMatrix, c(as.character(l[i, "GeneSymbol"]), l[i, "logFC"]))
    } else if (sign(as.double(signMatrix[position, 2])) + sign(l[i, "logFC"]) == 0){
      # remove all entries of this gene from l
      pos <- which(l[, "GeneSymbol"] == l[i, "GeneSymbol"])
      positions <- c(positions, pos)
      cat(sprintf("Dual sign Gene found and will be removed: %s\n", l[i, "GeneSymbol"])) # debug
    }
  }
  # if is necessary else it deletes everything
  if (length(positions) > 0) l <- l[-as.integer(positions), ]
  return(l)
}

removeDuplicates <- function(l){
  # sort by pvalue
  l <- l[order(l[, "p.value"], decreasing = FALSE), ]
  # remove duplicates on gene symbol
  l <- l[!duplicated(l[, "GeneSymbol"]),]
  cat(sprintf("Duplicate Gene Symbol entries removed.\n"))
  return(l)
}

removeDuplicates_padj <- function(l){
  # sort by pvalue adj
  l <- l[order(l[, "adj.p.value"], decreasing = FALSE), ]
  # remove duplicates on gene symbol
  l <- l[!duplicated(l[, "GeneSymbol"]),]
  cat(sprintf("Duplicate Gene Symbol entries removed.\n"))
  return(l)
}

sortByAbsLogFC <- function(l){
  l <- l[order(abs(l[, "logFC"]), decreasing = TRUE), ]
  cat(sprintf("Sorted by descending absolute logFC values.\n"))
  return(l)
}

#####
### MAIN ###
#####

# library
library(limma)

# directory
setwd("C:/Users/andreak/Desktop/GSE45848") 

# input variables
folder <- "C:/Users/andreak/Desktop/GSE45848"
series_matrix_file <- paste(folder, "series_matrix.txt", sep = "/")
classlabels_file <- paste(folder, "class_labels.txt", sep = "/")
platform_file <- paste(folder, "GPL887-20438_platforms.txt", sep = "/")

# output filename
#dataset <- "GSE85347_patients_vs_controls"
outfile <- paste(folder, "GSE45848_limma_v1.tsv", sep = "/")

# parse files, filename + boolean header
series_matrix <- parseFile(series_matrix_file, TRUE) # character or double
classlabels <- parseFile(classlabels_file, FALSE) # integer
platform <- parseFile(platform_file, TRUE) # character (watch out for dates!)
colnames(platform)=c("probeIDs","GeneSymbol")

#remove rows with empty entries 
series_matrix <- series_matrix[complete.cases(series_matrix), ]


# breaking series_matrix components to probe IDS and gene intensities
probeIDs <- matrix(as.character(series_matrix[, 1]))
colnames(probeIDs) <- "probeIDs"
temp <- as.matrix(series_matrix[,-1])
gene_intensities <- matrix(as.double(temp),
                           nrow = nrow(temp), ncol = ncol(temp))
colnames(gene_intensities) <- colnames(temp)
rownames(gene_intensities) <- probeIDs
remove(temp) # clearing memory

#check the max value in your gene intensities  
max(gene_intensities) 
summary(gene_intensities)

#if (min(gene_intensities) < 0) gene_intensities <- cureNegative(gene_intensities)

# boxplot and decide toNormalize = TRUE or FALSE
boxplot(gene_intensities) #, ylim=c(-0.1,0.1)
toNormalize <- FALSE

# check if intensities are log-ed
toLog <- checkLog(gene_intensities) # TRUE or FALSE

# handle input Data based on their values
gene_intensities <- handleData(gene_intensities, toNormalize, toLog)

gene_intensities <- gene_intensities[complete.cases(gene_intensities), ]

# check for negative values, 2 approaches
# 1. add abs min

if (min(gene_intensities) < 0) gene_intensities <- cureNegative(gene_intensities)
# or 2. change to threshold
# if (min(gene_intensities) < 0) gene_intensities[gene_intensities < 0] <- 0.001


# # # VERSION 1 with p-value

# run Limma with gene intensities and classlabels
limma <- runLimma(gene_intensities, classlabels)

# Check boxplot before running limma 
boxplot(gene_intensities) #, ylim=c(-0.1,0.1)

# clear p-value > 0.05
cutoff <- 0.05
limma <- pvalueClear(limma, cutoff)

#calculate adjusted p-value
#tT <- topTable(limma, adjust.method="fdr", number=nfow(limma))

# assigning gene symbols
colnames(platform)[1] <- "probeIDs"
limma <- merge(limma, platform)
# removing values without Gene Symbols first
limma <- limma[limma[, "GeneSymbol"] != "", ] # --- or " " or ""

# experimental function
limma <- removeDualSignProbes(limma)

# removing lower pvalue duplicates
limma <- removeDuplicates(limma)

# sorting by descending absolute logFC
limma <- sortByAbsLogFC(limma)

# printing outfile
write.table(limma, outfile, quote = FALSE, row.names = FALSE, sep = "\t")



# # # VERSION 2 with p-value adjusted

limma_a <- runLimma_padj(gene_intensities, classlabels)

# clear p-value > 0.05
cutoff <- 0.05
limma_a <- pvalueClear_padj(limma_a, cutoff)

#calculate adjusted p-value
#tT <- topTable(limma, adjust.method="fdr", number=nfow(limma))

# assigning gene symbols
colnames(platform)[1] <- "probeIDs"
limma_a <- merge(limma_a, platform)
# removing values without Gene Symbols first
limma_a <- limma_a[limma_a[, "GeneSymbol"] != "", ] # --- or " " or ""

# experimental function
limma_a <- removeDualSignProbes(limma_a)


# removing lower pvalue adj duplicates
limma_a <- removeDuplicates_padj(limma_a)

# sorting by descending absolute logFC
limma_a <- sortByAbsLogFC(limma_a)

# printing outfile
outfile_a <- paste(folder, "GSE45848_limma_a.tsv", sep = "/")
write.table(limma_a, outfile_a, quote = FALSE, row.names = FALSE, sep = "\t")

length(intersect(limma$GeneSymbol,limma_a$GeneSymbol))
#####
### SCRIPT END ###

