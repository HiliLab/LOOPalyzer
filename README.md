# LOOPalyzer (release 180924)
--------------------------------------------
LOOPalayzer is a Python 2.7 program for analyzing LOOPER-generated aptamer libraries. It was written by Wayland Yeung for the Hili Lab and has been tested working on Ubuntu 18.04 LTS.

## Publications
--------------------------------------------

Kong, D., Yeung, W., & Hili, R. (2017). In vitro selection of diversely functionalized aptamers. Journal of the American Chemical Society, 139(40), 13977-13980.

Kong, D., Yeung, W., Movahedi M., Chen D., & Hili, R. (2019). Evaluation of the evolutionary outcomes of LOOPER-derived aptamers isolated from in vitro evolution. Manuscript submitted.

## Prerequisites
--------------------------------------------
LOOPalyzer requires several Python libraries to run: NumPy, SciPy, and Biopython. These can be easily installed with pip by running the following commands:

    pip install biopython
    pip install numpy 
    pip install scipy

## Example
--------------------------------------------
If you're impatient, you can run the example with the following command:

    python loopalyzer.py r6_sample.fastq -op options.txt -b annnn.tsv

## Usage
--------------------------------------------
LOOPalyzer requires an options file which lists information about the 
sequence in plain text:

| Setting | Description |
| ------- | ----------- |
| pr_l    | LOOPER primer flanking the reading region on the **left** side, relative to the FASTQ file. |
| pr_r    | LOOPER primer flanking the reading region on the **right** side, relative to the FASTQ file. The reading region will be defined as the sequence in between pr_l and pr_r. |
| format  | LOOPER codon format. **A/C/G/T** denotes a constant position. **N** denotes a random nucleotide. **X** denotes a random nucleotide that also encodes the LOOPER modification. |
| codons  | Number of codons in the reading region. The length of the reading region will be assumed to be the number of codons multiplied by the length of the codon format |

If you have prior knowledge of your LOOPER library's codon distribution, LOOPalyzer can calculate the conservation of each codon position. Using Jensen-Shannon divergence, it compares the distribution of codons under nonselective conditions against the sequences being analyzed. This distribution of codons under nonselective conditions can be given tab 
delimited format. The distribution is calculated by percentage of the whole, 
so absolute counts are acceptable. **This step is optional.**

LOOPalyzer accepts either fastq or fasta formatted reads as input. The  program can be run with the following command:

    python loopalyzer.py <input> -op <options_file> -b <codon_distribution>

If you do not know your library's prior codon distribution, simply run the program without passing it:

    python loopalyzer.py <input> -op <options_file>


## Output
--------------------------------------------

All outputs will be in plain text format.

| Output File| Description |
| ---------- | ----------- |
| filename.count.tsv      | Lists all unique sequences in the sequencing file as tab separated values: field 1 is absolute counts, field 2 is proportion, field 3 is the unique sequence. |
| filename.reads.tsv      | Lists all sequences in the same order as <input>.counts.tsv: field 1 is absolute counts, field 2 is error analysis, field 3 is reading region. |
| filename.codons.txt     | Shows the distribution of codons found in the sequences and also shows a breakdown of errors found in the sequences. |
| filename.divergence.txt | This file is output if a reference distribution is provided. Shows the divergence of each codon position against the nonselective set |



