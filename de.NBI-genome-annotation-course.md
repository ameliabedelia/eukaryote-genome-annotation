![](images/ikmb_bfx_logo.png)
 
# Overview

The following exercices will introduce you to methods in eukaryotic genome annotation. The following topics will be discussed in order:

**A. Crash course/lecture**  
 * Insights into eukryotic genome annotation *(not included in this wiki)*  
  
  
**B. Practical course**  
 * **B1.** Genome Browsers  
 * **B2.** Understanding male morphs in the ruff (*Philomachus pugnax*) - A genomics approach  

Note: the exercises are split into two separate tasks (B1 and B2). We suggest that you take a look at B2 first, as it requires a bit more runtime on your computer. While B2 is running, you can get started on B1.  

For the exercises in B2, you will need the **IGV genome browser**. The program can be copied from this location to your home folder */mnt/lectures/biol258/Day7/IGV* using:  
`cp -R /mnt/lectures/biol258/Day07/IGV $HOME`  

*Please note! The commands given in these pages are examples and they assume that you make sure that the input and output files are where they should be.*  
*Some times you may have to copy files to your current directory or expand the example commands to include correct/full paths!*  
  
### Software containers  
Many of the tools used in this exercise are pretty complex to install. In order to avoid having to install them on each computer we are going to use software *containers*. These are *portable* environments in which you can install an operating system, a program or many programs and then copy the whole thing to whichever computer you like to run the software. Popular solutions are Docker and **Singularity** (which we will use here). The way to run a tool inside a container is pretty straight-forward:  

If you would run a command without container like this:  
`samtools sort -o rnaseq.sorted.bam rnaseq.bam`  

You would run the same command with a container like this:  
`singularity exec container.img samtools sort -o rnaseq.sorted.bam rnaseq.bam`  

Here `container.img` is a place holder, which you have to replace by the name of the actual container.  

Copy the singularity image for this exercise to your local hard drive:  
`cp -R /mnt/lectures/biol258/Day07/singularity.img $HOME`  

In order to be able to run the examples in the following exercises, you will have to add the following code in front of each command:  
`singularity exec $HOME/singularity.img <command execution here>` 

# B1. Genome Browsers  
Genome browsers are an important tool for anyone working with eukaryote genome sequences. They can help visualize rich data sets and are often critical to the biological interpretation of findings from e.g. sequencing experiments.  
It is important to familiarize yourself with these resources so you understand what kind of features can be annotated and where to find this type of information.  
Try out the following URLs to get an impression of different Genome-Browser flavours. What features/annotations do these browsers offer? How do they display them? If you need some inspiration regarding which genes to inspect in EnsEMBL and UCSC, you can start with these:  
BRCA1, BRCA2, HOXA1  
* Ensembl Genome Browser (http://www.ensembl.org/index.html)
* UCSC Genome Browser (http://genome.ucsc.edu)
* NCBI's Genome Data Viewer (https://www.ncbi.nlm.nih.gov/genome/gdv/)
* Metazome (http://www.metazome.net)  
Once you have done that, try completing the following tasks using the EnsEMBL genome portal liked above.  

**QB1.1:** How many protein-coding genes are annotated in the mouse (check the info page for mouse, accessible from the front page and then search for annotation statistics).  

**QB1.2:** How many non-coding genes are annotated in human? And in chicken?  

**QB1.3:** How many transcript variants are annotated in EnsEMBL for the human BRCA1 gene? How long is the longest protein product, how long the shortest?  
