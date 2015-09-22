library(DESeq2)
directory = 'C:/RNASeq/ReadCounts/'
setwd(directory)

# Load data directly from HTSeq Count
sampleFiles=grep('.geneIDCounts.txt',list.files(directory),value=TRUE)
sampleNames=sub('.geneIDCounts.txt','',sampleFiles)
sampleCondition=sub('[1-3]','',sampleNames)
sampleTable=data.frame(sampleName=sampleNames, fileName=sampleFiles, condition=sampleCondition)
ddsHTSeq = DESeqDataSetFromHTSeqCount(sampleTable=sampleTable, directory=directory, design= ~ condition)

# Differential expression analysis
dds <- DESeq(ddsHTSeq)
dds$condition = relevel(dds$condition, 'Control')
res <- results(dds)

# Extract the significantly differentially expressed genes
resOrdered<- res[order(res$padj),]
resSig0.05 <- subset(resOrdered, padj<0.05)

# Print results to file
setwd('C:/RNASeq/DESeq2Output/')
write.table(as.data.frame(res), file='ExperimentalvsControl_DEResults.txt',sep='\t',quote=FALSE)
setwd('C:/RNASeq/DESeq2DEGenes/')
write.table(as.data.frame(resSig0.05), file='ExperimentalvsControl_DE_pAdj0.05.txt',sep='\t',quote=FALSE)
