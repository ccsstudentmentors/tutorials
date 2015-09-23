Now that you have your read data aligned to the genome, it is time to quantify the expression of the measured RNA and perform differential expression testing! 

At this point, you need to choose whether you want to analyze your data at the transcript level or at the gene level (although you probably made this decision long ago when you actually prepared the RNA for RNA-Seq). 

If you already know which route you want to take, head to that folder now and continue processing your data! If not, here is a brief list of the pros and cons of each approach to help you decide which better suits your scientific goals. 

# Gene Level Analysis:
Pros:
Easier and faster to perform the analysis.
More straightforward results.
Results more likely to be reproducible/validatable than transcript level results

Cons:
The most accurate gene level analysis programs only examine the expression of a gene across samples -- there is no way to compare the expression of one gene to another gene within the same sample.
The gene level analysis programs we will be covering require you to work in yet another programming environment -- R. 

# Transcript Level Analysis: 
Pros:
You get to see an entire additional level of regulation beyond the gene level.
You can still observe gene level changes in expression -- you just also have more information about what is happening to different forms of that gene.

Cons:
You get substantially less power in testing differential expression (partly because there will always be more transcripts than genes, and partly just due to a weakness in cufflinks itself).
Results are often difficult to interpret and difficult to validate in the lab.
