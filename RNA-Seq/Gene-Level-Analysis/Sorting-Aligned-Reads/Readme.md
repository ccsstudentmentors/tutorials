This step is easy! 

All we need to do is take the aligned read data produced by tophat (or STAR) and sort it.
Right now, that data is just in whatever order the alignment program figured out its genomic site in. 
To efficiently count how many reads we have for each gene (which will be our next step), the program needs the reads to be sorted in a particular way.
We have two choices for how to sort the reads: by their position on the chromosome, or by the name of the read.

Intuitively, sorting by position may seem like the far more useful mode to choose, but this tutorial is actually going to use sorting by name instead.
The reason for this is because the program that does the counting over genomic intervals (HTSeq-count) handles paired ended reads much more memory efficiently when the data are name sorted than if they are position sorted.
If your data is single-ended, you are welcome to sort them either way (but name sorting will still work fine, and that is the way we will be going over here).

To sort the bam file data, we need to load a suite of tools called samtools using the following command:
module load samtools/0.1.19

You may have already loaded samtools to peek at the bam file using its 'view' command. But now we are going to use its 'sort' tool.

You can find the full manual on samtools' commands here: http://www.htslib.org/doc/samtools-1.1.html
See the file 'Explanation of BamSort_Example.txt' and the 'BamSort_Example.job' files for an example of how we will sort the data.
