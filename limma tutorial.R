
library(limma)
library(BiocManager)
library('ChimpHumanBrainData')
library(affy)
celfileDir = system.file('extdata',package='ChimpHumanBrainData')
celfileNames = list.celfiles(celfileDir)
brainBatch=ReadAffy(filenames=celfileNames,celfile.path=celfileDir,compress=TRUE)

targets <- readTargets()


# TUTORIAL  
# dataset with 7 brain regions . 
# 3 human brains and 3 chimpanzee
# 

sampleNames(brainBatch)
sampleNames(brainBatch)=paste(rep(c("CH","HU"),each=12),rep(c(1:3,1:3),each=4), 
                              rep(c("Prefrontal","Caudate","Cerebellum","Broca"),6),sep="")

brain.expr=exprs(brainBatch)
library(hexbin)
plot(hexplom(log2(brain.expr[,paste(rep(c("CH","HU"),each=3),
                                    c(1:3,1:3),rep("Prefrontal",6),sep="")])))

trts=factor(paste(rep(c("CH","HU"),each=12),
                  rep(c("Prefrontal","Caudate","Cerebellum","Broca"),6),sep=""))
blocks=factor(rep(1:6,each=4))

# Normalize the expression values and combine into probeset summaries 
brain.rma=rma(brainBatch)


# LIMMA analysis 
# 1. Compute the s2p pooled variance. Compute within region variance for each gene 
# and the correlation among regions from the same brain 

design.trt=model.matrix(~0+trts)

# if there are blocks or technical replicates the correlation of genes 
# within the blocks needs to be computed 

# in our case the blocking factor is called blocks
library(statmod)
corfit <- duplicateCorrelation(brain.rma, design.trt, block = blocks)
hist(tanh(corfit$atanh.correlations))


# Comments; correlations are mainly positive with mode around 0.6 

# COmpute pool sample variance for each gene 
fitTrtMean <- lmFit(brain.rma, design.trt, block = blocks, cor = corfit$consensus.correlation)
fitTrtMean$coefficients
fitTrtMean$sigma


# 2 components : coefficients: contains the mean expression for each gene in each treatmet
# sigma: the estimate of Sp the pooled SD


# 2. Create the coefficient matrix for the contrasts 
# 1. Average chimpanzee versus average human 
# 2. Chimpanzee versus human for each region 
# 3. The interaction between species and the comparison of cerebellum to Broca’s region

colnames(design.trt)
contrast.matrix=makeContrasts( (trtsCHBroca+trtsCHCaudate+trtsCHCerebellum+trtsCHPrefrontal)/4 -(trtsHUBroca+trtsHUCaudate+trtsHUCerebellum+trtsHUPrefrontal)/4, trtsCHBroca-trtsHUBroca, trtsCHCaudate-trtsHUCaudate, trtsCHCerebellum-trtsHUCerebellum, trtsCHPrefrontal-trtsHUPrefrontal, (trtsCHCerebellum-trtsHUCerebellum)-(trtsCHBroca-trtsHUBroca), levels=design.trt) 
colnames(contrast.matrix)= c("ChVsHu","Broca","Caudate","Cerebellum","Prefrontal","Interact")


# 3. Compute the estimated contrasts 

fit.contrast=contrasts.fit(fitTrtMean,contrast.matrix)



# 4. Compute the moderated contrast t-test for each gene 
# compute consensus pooled variance, and use it to compute empirical Bayes (moderated pooled variance for each gene )


efit.contrast=eBayes(fit.contrast)


# 5. Plot the histogram of p-values for each contrast for each gene
par(mfrow=c(2,3))
for (i in 1:ncol(efit.contrast$p.value)) { 
  hist(efit.contrast$p.value[,i],main=colnames(efit.contrast$p.value)[i])
  }


# 6. Create the list of significant genes based on th ep-values 

genes=geneNames(brainBatch)
topTable(efit.contrast,coef=1,adjust.method="BY",n=10,p.value=1e-5,genelist=genes)



write.table(file="fits.txt",  cbind(genes,fitTrtMean$coefficients,efit.contrast$coefficients,efit.contrast$p.value),  row.names=F,  col.names=c("GeneID",colnames(fitTrtMean$coefficients),colnames(efit.contrast$p.value), paste("p",colnames(efit.contrast$coefficients))),sep=",")



