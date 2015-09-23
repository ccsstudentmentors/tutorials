Now that you have sorted your aligned read bam files, it is time to count how many reads aligned to the genomic regions belonging to each gene.

While cufflinks takes a very intense statistical approach to resolving different isoforms and calculating differential expression, gene level analysis is much more straightforward:
For each read that was produced by the HiSeq machine for this experiment, we align it to the genome. (Done!)
Then we ask if that place where it aligned is an exon of a gene (this is our current step). If so, we give a 'point' to that gene.
We repeat this counting process for all of our tens of millions of reads for each sample. At the end, for each sample, we have a list of genes and their corresponding number of read counts (points) for this sample.
Then we use our replicate samples of our control and experimental conditions to calculate the average number of read counts for each gene for the control and experimental groups (and the variance of the counts around that average). Then we are basically just performing a big statistical comparison (similar to an ANOVA) to ask if the number of read counts assigned to each particular gene changed significantly between the control and experimental groups (controlling for the different total number of reads in each sample). 

Ok, so performing this counting of the reads could be really difficult and annoying to do, but fortunately someone (named Simon Anders) has written an awesome program that does this all for us. That program is called HTSeq-count.

HTSeq-count is a python script that is already available on Pegasus' installation of Python version 2.7.3.
Python is one of the only programming languages where nobody uses the latest version. There is a Python 3, but most scripts are still written for Python 2.7.x, including this one. So that is what we are going to use.

Get started by loading the python module:
module load python/2.7.3

Great! Now, take a look at the 'GetReadCounts_Example.job' and 'Explanation of GetReadCounts job file.txt' files. 
