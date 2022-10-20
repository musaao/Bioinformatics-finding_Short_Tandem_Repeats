# Bioinformatics-finding_Short_Tandem_Repeats

This is a project for a graduate course (Bioinformatics)

>The fasta files contain the sequences considered for running the code on.

>The code for this project was written in jupyter notebook (.ipynb) and exported as a python file (.py)

>dictionaryforproject.txt shows all combination of di- and tri- nucleotides screened for in the sequences. 
These occur in repeats of at least 3 (e.g "ATATAT") to be considered as short tandem repeats (STRs)

### Output
e.g.

_4 STRs ending at pos. 819_

_3 STRs ending at pos. 1047_

_3 STRs ending at pos. 8688_

_TAA total STRs in the sequence : 10_

_TAA binomial distribution, Pb: 1.4803323638208016e-69_

>In a DNA sequence of 15000+ nucleotides, 4 repeats was found ending at position 819 in the sequence, 3 ..., 3... (all belonging to TAA); 
all these form a total of 10 STRs in the entire sequence. Pb shows the binomial distribution of such repeat in the sequences.
