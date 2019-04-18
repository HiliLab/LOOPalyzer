############################################################################
# LOOPalyzer (release 180924)
# Wayland Yeung
############################################################################

LOOPalayzer is a program for analyzing LOOPER-generated aptamer libraries.
It is written in Python 2.7 and has been tested working on Ubuntu 18.04 LTS.

How to run the example:

>> python loopalyzer.py r6_sample.fastq -op options.txt -b annnn.tsv

############################################################################

1) LOOPalyzer requires the NumPy, SciPy, and Biopython modules. If any of 
these are not installed, please install it with its respective commands:

>> pip install biopython
>> pip install numpy
>> pip install scipy


2) LOOPalyzer requires an options file which lists information about the 
sequence in plain text:

pr_l:   the left side primer (relative to the sequencing file) used in
        LOOPER polymerization.
pr_r:   the right side primer (relative to the sequencing file) used in 
        LOOPER polymerization. The reading region is defined as the sequence 
		between pr_l and pr_r
format: the format of the LOOPER codon where A/C/G/T denotes a constant 
        position. N denotes a random nucleotide ;X denotes a modification 
        defining base pair
        example: ANNXX was the original pentamer system used against thrombin
codons: number of codons in the reading region, the length of the reading 
        region is assumed to be the number of codons multiplied by the length 
		of the given format

		
3) This step is optional. LOOPalyzer also calculates the conservation of 
each codon position by Jensen-Shannon divergence. It compares the distribution
of codons under nonselective conditions against the sequences being analyzed. 
This distribution of codons under nonselective conditions can be given tab 
delimited format. The distribution is calculated by percentage of the whole, 
so absolute counts are acceptable.


4) LOOPalyzer accepts either fastq or fasta formatted reads as input. The 
program can be run with the following command:

>> python loopalyzer.py <input> -op <options_file> -b <codon_distribution>

All outputs are in plain text and are listed below:

<input>.count.tsv :
  Lists all unique sequences in the sequencing file as tab separated values
  field 1 is absolute counts, field 2 is proportion, field 3 is the unique sequence
	
<input>.reads.tsv :
  Lists all sequences in the same order as <input>.counts.tsv
  field 1 is absolute counts, field 2 is error analysis, field 3 is reading region
  
<input>.codons.txt :
  Shows the distribution of codons found in the sequences
  Also shows a breakdown of errors found in the sequences
  
<input>.divergence.txt : 
  This file is output if a reference distribution is provided
  Shows the divergence of each codon position against the nonselective set
