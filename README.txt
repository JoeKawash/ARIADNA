ARIADNA: Machine Learning Method for Ancient DNA Variant Discovery

ARIADNA (ARtificial Intelligence for Ancient DNA) is a machine learning based approach for detecting single-nucleotide variants in read alignments of ancient DNA samples. ARIADNA utilizes a pre-built decision tree model for efficient and accurate variant calling. 

Source code can be found at:https://github.com/JoeKawash/ARIADNA or https://osf.io/5bph4/

Contact: jokawash@scarletmail.rutgers.edu

Requirements:

Linux operating system
python 2.7
python module scikit-learn up to version 0.18.2 (http://scikit-learn.org/stable/) 
GROM-ANC_v1_0_2 (https://osf.io/5bph4/)

All files from the ARIADNA package should be unpacked within the same directory. 

Input:

ARIADNA requires a reference genome, sorted BAM file, and associated BAM index file. 

Example of Usage: 

test files are in /test
validation file is in /validate

> ARIADNA.py test/test.bam test/test.fasta validate/trial_run.vcf

trial_run.vcf should list the same variants as ARIADNA_validate.vcf
ARIADNA works best in data with moderate to high read coverage and scaffolds of at least 1MB.

Please cite:
https://doi.org/10.1093/dnares/dsy029

Kawash, Joseph K., et al. "ARIADNA: machine learning method for ancient DNA variant discovery." DNA Research 25.6 (2018): 619-627.
