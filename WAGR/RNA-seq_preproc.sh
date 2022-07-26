#!bin/bash

conda activate rna-seq

##get data from SRA

echo "Starting prefetch"
seq -f %1.0f 3371765  3371768 | while read x; do echo GSM${x} ;done |  parallel -j0 --bar 'prefetch  {.}'
echo "DONE prefetch"

mkdir -p raw_data
echo "fastq-dump for paired-reads"
find . | grep [0-9].sra | while read x; do a=$(basename $x); sub=${a%%.*} ; echo fastq-dump  --defline-qual "+" -v --split-files  -O raw_data/ --gzip ${sub}  ;done
echo "DONE fastq-dump"


mkdir -p QC
echo "fastqc and multiqc"
ls *fastq.gz | while read x; do fastqc -t 6 -o QC/ ${x} ;done
cd QC
multiqc .
cd ..
echo "DONE QC"

echo "building STAR index"
STAR --runThreadN 48 --runMode genomeGenerate --genomeDir GRCh38.primary_assembly.index --genomeFastaFiles GRCh38.primary_assembly.genome.fa --sjdbGTFfile gencode.v26.primary_assembly.annotation.gtf --sjdbOverhang 99
echo "Done building STAR index"

mkdir -p results
echo "STAR mapping"
cd .. raw_data/
ls *fastq.gz | while read x; do gunzip ${x} & done
cd ..
ulimit -n 10000
find . | grep 1.fastq$ | while read x; do a=$(basename ${x}) ; sub=${a%%_*} ; STAR --genomeDir GRCh38.primary_assembly.index/ --runThreadN 48 --readFilesIn raw_data/${sub}_1.fastq raw_data/${sub}_2.fastq --outSAMtype BAM SortedByCoordinate --outSAMunmapped Within --outSAMattributes Standard --outFileNamePrefix ${sub} ;done
ulimit -n 1024
echo "Done STAR mapping"

echo "featureCounts"
cd results/
featureCounts 	-p --countReadPairs \
	     	-a ../gencode.v26.primary_assembly.annotation.gtf \
	     	-T 48 \
		-t exon \
		-g gene_id \
		-o counts.out  \
		*bam

cat counts.out | cut -f 1,7- > counts_extracted.out
echo "Done featureCounts"

cd ..
