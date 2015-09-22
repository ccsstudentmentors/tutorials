So you have chosen to work with your data at the gene level rather than the isoform level? 
Great choice! The interpretation of your results will be far more straightforward and will require much less difficult validation than if you had ventured down to the isoform level.

There are many tools to perform Gene-Level analyses of RNA-Seq data, but I will just cover two of the most popular options here: edgeR and DESeq2. 

Here is a basic overview of what Gene-Level analysis entails:
We need to sort our bam file of aligned reads (we will use samtools sort)
We need to count how many reads were aligned to each gene in our genomic annotation (we will use HTSeq-count)
We can then calculate the differentially expressed genes for our dataset using either edgeR or DESeq2 (or many other programs that I won't discuss)

So, get started by heading to the folder 'Sorting Aligned Reads'. 
