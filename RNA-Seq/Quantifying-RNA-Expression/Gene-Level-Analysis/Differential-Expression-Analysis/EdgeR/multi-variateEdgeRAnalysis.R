library(edgeR)
directory = 'C:/RNASeq/ReadCounts/'
setwd(directory)

# load each sample individually
Control1 = read.table('Control1.geneIDCounts.txt',sep='\t',header=FALSE)
Control2 = read.table('Control2.geneIDCounts.txt',sep='\t',header=FALSE)
Control3 = read.table('Control3.geneIDCounts.txt',sep='\t',header=FALSE)
TreatmentA1 = read.table('TreatmentA1.geneIDCounts.txt',sep='\t',header=FALSE)
TreatmentA2 = read.table('TreatmentA2.geneIDCounts.txt',sep='\t',header=FALSE)
TreatmentA3 = read.table('TreatmentA3.geneIDCounts.txt',sep='\t',header=FALSE)
TreatmentB1 = read.table('TreatmentB1.geneIDCounts.txt',sep='\t',header=FALSE)
TreatmentB2 = read.table('TreatmentB2.geneIDCounts.txt',sep='\t',header=FALSE)
TreatmentB3 = read.table('TreatmentB3.geneIDCounts.txt',sep='\t',header=FALSE)

# combine data sets into a matrix
geneCounts = data.frame(Control1[,2],Control2[,2],Control3[,2],TreatmentA1[,2],TreatmentA2[,2],TreatmentA3[,2],TreatmentB1[,2],TreatmentB2[,2],TreatmentB3[,2])
row.names(geneCounts) = Control1[,1]
sizeGeneCounts = dim(geneCounts)
geneCounts = geneCounts[1:(sizeGeneCounts[1]-5),]
condition = c(rep('Control',3),rep('TreatmentA',3),rep('TreatmentB',3))
sampleNames = c('Control1','Control2','Control3','TreatmentA1','TreatmentA2','TreatmentA3','TreatmentB1','TreatmentB2','TreatmentB3')
colnames(geneCounts) = sampleNames

# build the generalized linear model that will be used for differential expression testing
dge <- DGEList(counts=geneCounts, group=condition)
design <- model.matrix(~condition+0, data=dge$samples)
colnames(design) = gsub("condition","",colnames(design))
disp <- estimateGLMCommonDisp(dge, design)
disp <- estimateGLMTrendedDisp(disp, design)
disp <- estimateGLMTagwiseDisp(disp, design)
fit <- glmFit(disp, design)

# tell edgeR what comparisons you want to perform
TreatmentAVsControl = makeContrasts(TreatmentA-Control, levels=design)
TreatmentBVsControl = makeContrasts(TreatmentB-Control, levels=design)

# perform the differential expression testing for those comparisons
lrt.TreatmentAVsControl = glmLRT(fit, contrast=TreatmentAVsControl)
res1 = as.data.frame(lrt.TreatmentAVsControl$table)
lrt.TreatmentBVsControl = glmLRT(fit, contrast=TreatmentBVsControl)
res2 = as.data.frame(lrt.TreatmentBVsControl$table)

# extract the significantly differentially expressed genes
res1Ordered = res1[order(res1$PValue),]
res1Sig = subset(res1Ordered, PValue<0.05)
res2Ordered = res2[order(res2$PValue),]
res2Sig = subset(res2Ordered, PValue<0.05)

# Convert Gene IDs to Gene Names (if applicable)
conversionTable=read.table('C:\RNASeq\GeneIDtoNameConversionTable.txt',sep='\t',header=FALSE,row.names=1)
colnames(conversionTable) = c('GeneName')
nameRes1=merge(conversionTable,res1,by.x=0,by.y=0)
nameRes1Sig=merge(conversionTable,res1Sig,by.x=0,by.y=0)
nameRes2=merge(conversionTable,res2,by.x=0,by.y=0)
nameRes2Sig=merge(conversionTable,res2Sig,by.x=0,by.y=0)

# Print results to file
setwd('C:/RNASeq/EdgeROutput/')
write.table(nameRes1, file='TreatmentAvsControl_DEResults.txt',sep='\t',quote=FALSE)
write.table(nameRes2, file='TreatmentBvsControl_DEResults.txt',sep='\t',quote=FALSE)
setwd('C:/RNASeq/EdgeRDEGenes/')
write.table(nameResSig1, file='TreatmentAvsControl_DE_pVal0.05.txt',sep='\t',quote=FALSE)
write.table(nameResSig2, file='TreatmentBvsControl_DE_pVal0.05.txt',sep='\t',quote=FALSE)
