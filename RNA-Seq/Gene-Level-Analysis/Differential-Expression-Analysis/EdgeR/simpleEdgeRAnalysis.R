library(edgeR)
directory = 'C:/RNASeq/ReadCounts/'
setwd(directory)

# load each sample individually
Control1 = read.table('Control1.geneIDCounts.txt',sep='\t',header=FALSE)
Control2 = read.table('Control2.geneIDCounts.txt',sep='\t',header=FALSE)
Control3 = read.table('Control3.geneIDCounts.txt',sep='\t',header=FALSE)
Experimental1 = read.table('Experimental1.geneIDCounts.txt',sep='\t',header=FALSE)
Experimental2 = read.table('Experimental2.geneIDCounts.txt',sep='\t',header=FALSE)
Experimental3 = read.table('Experimental3.geneIDCounts.txt',sep='\t',header=FALSE)

# combine data sets into a matrix
geneCounts = data.frame(Control1[,2],Control2[,2],Control3[,2],Experimental1[,2],Experimental2[,2],Experimental3[,2])
row.names(geneCounts) = Control1[,1]
sizeGeneCounts = dim(geneCounts)
geneCounts = geneCounts[1:(sizeGeneCounts[1]-5),]
condition = c(rep('Control',3),rep('Experimental',3))
sampleNames = c('Control1','Control2','Control3','Experimental1','Experimental2','Experimental3')
colnames(geneCounts) = sampleNames

# build the generalized linear model that will be used for differential expression testing
dge <- DGEList(counts=allGeneCounts, group=condition)
design <- model.matrix(~group+0, data=dge$samples)
colnames(design) = gsub("group","",colnames(design))
disp <- estimateGLMCommonDisp(dge, design)
disp <- estimateGLMTrendedDisp(disp, design)
disp <- estimateGLMTagwiseDisp(disp, design)
fit <- glmFit(disp, design)

# tell edgeR what comparison you want to perform
ExperimentalVsControl = makeContrasts(Experimental-Control, levels=design)

# perform the differential expression testing for that comparison
lrt.ExperimentalVsControl = glmLRT(fit, contrast=ExperimentalVsControl)
res = as.data.frame(lrt.ExperimentalVsControl$table)

# extract the significantly differentially expressed genes
resOrdered = res[order(res$PValue),]
resSig = subset(resOrdered, PValue<0.05)

# Convert Gene IDs to Gene Names (if applicable)
conversionTable=read.table('C:\RNASeq\GeneIDtoNameConversionTable.txt',sep='\t',header=FALSE,row.names=1)
colnames(conversionTable) = c('Gene Name')
nameRes=merge(conversionTable,res,by.x=0,by.y=0)
nameResSig=merge(conversionTable,resSig,by.x=0,by.y=0)

# Print results to file
setwd('C:/RNASeq/EdgeROutput/')
write.table(nameRes, file='ExperimentalvsControl_DEResults.txt',sep='\t',quote=FALSE)
setwd('C:/RNASeq/EdgeRDEGenes/')
write.table(nameResSig, file='ExperimentalvsControl_DE_pVal0.05.txt',sep='\t',quote=FALSE)
