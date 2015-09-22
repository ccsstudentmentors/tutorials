library(DESeq2)
directory = 'C:/RNASeq/ReadCounts/'
setwd(directory)

# load each sample individually
Control1 = read.table('Control1.geneNameCounts.txt',sep='\t',header=FALSE)
Control2 = read.table('Control2.geneNameCounts.txt',sep='\t',header=FALSE)
Control3 = read.table('Control3.geneNameCounts.txt',sep='\t',header=FALSE)
Experimental1 = read.table('Experimental1.geneNameCounts.txt',sep='\t',header=FALSE)
Experimental2 = read.table('Experimental2.geneNameCounts.txt',sep='\t',header=FALSE)
Experimental3 = read.table('Experimental3.geneNameCounts.txt',sep='\t',header=FALSE)

# combine data sets into a matrix
geneCounts = data.frame(Control1[,2],Control2[,2],Control3[,2],Experimental1[,2],Experimental2[,2],Experimental3[,2])
row.names(geneCounts)=Control1[,1]
condition = c(rep('Control',3),rep('Experimental',3))
sampleNames=c('Control1','Control2','Control3','Experimental1','Experimental2','Experimental3')
colnames(geneCounts) = sampleNames
colData= data.frame(condition,row.names=sampleNames)
dds = DESeqDataSetFromMatrix(countData = geneCounts, colData=colData, design = ~ condition)

# Differential expression analysis
dds <- DESeq(dds)
dds$condition <- relevel(dds$condition, 'Control')
res <- results(dds)

# Extract the significantly differentially expressed genes
resOrdered<- res[order(res$padj),]
resSig <- subset(resOrdered, padj<0.05)

# Print results to file
setwd('C:/RNASeq/DESeq2Output/')
write.table(res, file='ExperimentalvsControl_DEResults.txt',sep='\t',quote=FALSE)
setwd('C:/RNASeq/DESeq2DEGenes/')
write.table(resSig, file='ExperimentalvsControl_DE_pAdj0.05.txt',sep='\t',quote=FALSE)
