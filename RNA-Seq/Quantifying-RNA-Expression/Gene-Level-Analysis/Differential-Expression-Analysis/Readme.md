Now that you have counted how many reads aligned to each gene for each sample, you can perform your gene-level differential expression testing!

You have a choice for this final step though -- edgeR or DESeq2 (or many others that I won't cover)?

Honestly, both programs are great. They both use relatively similar statistical models to analyze the data. Many people run both analyses and simply use the differentially expressed genes found by both programs. So proceed however you see fit.

If you are looking for the easiest one to start with, I'd go with DESeq2. If you followed my earlier advice and used geneIDs for the read counting step, then running DESeq2 requires very, very few lines of code. 

Both analyses start with the set of files containing the read counts per gene for each sample. 

So, just head to the folder of the analytical program you want to use. 
