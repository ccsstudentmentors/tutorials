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
resA <- results(dds, contrast=c('condition','TreatmentA','Control'))
resB <- results(dds, contrast=c('condition','TreatmentB','Control'))

# Extract the significantly differentially expressed genes
resAOrdered<- resA[order(resA$padj),]
resASig <- subset(resAOrdered, padj<0.05)
resBOrdered<- resB[order(resB$padj),]
resBSig <- subset(resBOrdered, padj<0.05)

# Convert Gene IDs to Gene Names
conversionTable=read.table('C:\RNASeq\GeneIDtoNameConversionTable.txt',sep='\t',header=FALSE,row.names=1)
colnames(conversionTable) = c('GeneName')
nameResA=merge(conversionTable,as.data.frame(resA),by.x=0,by.y=0)
nameResASig=merge(conversionTable,as.data.frame(resASig),by.x=0,by.y=0)
nameResB=merge(conversionTable,as.data.frame(resB),by.x=0,by.y=0)
nameResBSig=merge(conversionTable,as.data.frame(resBSig),by.x=0,by.y=0)

setwd('C:/RNASeq/DESeq2Output/')
write.table(nameResA, file='TreatmentAvsControl_DEResults.txt',sep='\t',quote=FALSE)
write.table(nameResB, file='TreatmentBvsControl_DEResults.txt',sep='\t',quote=FALSE)
setwd('C:/RNASeq/DESeq2DEGenes/')
write.table(nameResSigA, file='TreatmentAvsControl_DE_pAdj0.05.txt',sep='\t',quote=FALSE)
write.table(nameResSigB, file='TreatmentBvsControl_DE_pAdj0.05.txt',sep='\t',quote=FALSE)
