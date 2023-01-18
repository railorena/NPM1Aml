# Pipeline to investigate RNAseq data

A reference-free pipeline to investigate RNAseq data using k-mers and Machine Learning methods to better understand prognosis classification and key mechanisms of AML pathogenesis.


##Commands
- Kmtricks
It's needed:
	samples file
	fastq fil
 '''
~/kmtricks pipeline --file samples --run-dir temp --kmer-size 31 -t 4 --until merge
~/kmtricks aggregate --run-dir ./temp --matrix kmer --format text --sorted -t 4 > sorted_matrix.tsv
'''
