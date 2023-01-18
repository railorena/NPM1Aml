# Pipeline to investigate RNAseq data

A reference-free pipeline to investigate RNAseq data using k-mers and Machine Learning methods to better understand prognosis classification and key mechanisms of AML pathogenesis.

## Required
Kmtricks tool [git](<https://github.com/tlemane/kmtricks/>)

R version 4.1.2: Package seqinr [cran](<https://cran.r-project.org/web/packages/seqinr/index.html>)

C++: Boost String Algorithms Library (doc] (<https://www.boost.org/doc/libs/1_81_0/doc/html/string_algo/reference.html#header.boost.algorithm.string_hpp>)


## Running
**Kmtricks**

Necessary files:
- samples
- fastq files

```
~/kmtricks pipeline --file samples --run-dir temp --kmer-size 31 -t 4 --until merge
~/kmtricks aggregate --run-dir ./temp --matrix kmer --format text --sorted -t 4 > sorted_matrix.tsv
```

**Adapting matrix script**

Necessary files:
- samples
- sorted_matrix.tsv
```
./adapt_mat.sh
```

**Coefficient of variation code**

Necessary files:
- sorted_matrix.tsv.gzip
- samples_cond
```
./cvar
tsv_to_fasta.R
```

**STAR and Samtools**

Necessary files:
- output.fasta
```
STAR --genomeDir /data/indexes/STAR/2.7.6a/GRCh38.97 --readFilesType Fastx  --readFilesIn output.fa --outFileNamePrefix leucegene_ --outStd Log --runMode alignReads --runThreadN 10 --outSAMunmapped Within
samtools view -b output_Aligned.out.sam |samtools sort > output_Aligned.out.bam	
```

**Generating model**

Necessary files:
- train.tsv
```
./xgb_classification.py
```



Note: examples of the files can be found in the folder *examples*
