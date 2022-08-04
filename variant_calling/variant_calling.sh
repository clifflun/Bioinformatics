#!usr/bin/bash

#data acquisition and preprocessing
prefetch -p -C yes -c SRR044234

fasterq-dump SRR044234/SRR044234.sra --qual-defline "+"

fastqc -t 16 *fastq

multiqc .

trimmomatic PE -threads 24 raw_data/SRR044234_1.fastq raw_data/SRR044234_2.fastq output_forward_paired.fq.gz output_forward_unpaired.fq.gz output_reverse_paired.fq.gz output_reverse_unpaired.fq.gz ILLUMINACLIP:TruSeq3-PE.fa:2:30:10:2:True SLIDINGWINDOW:4:20

#index building and alignment using Bowtie2

wget http://ftp.ensembl.org/pub/grch37/current/fasta/homo_sapiens/dna/Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz

bowtie2-build -f -q  --threads 48 Homo_sapiens.GRCh37.dna.chromosome.22.fa.gz hg19.chr22 ; echo "done"

bowtie2 -x ./index/hg19.chr22 -1 output_forward_paired.fq -2 output_reverse_paired.fq --fr -p 24 -S output_unsorted.sam --no-unal &

samtools view -h -b -o output_unsorted.bam output_unsorted.sam

sambamba sort -t 24 -o output_sorted.bam output_unsorted.bam

#variant calling using BCFtools

bcftools mpileup -f index/Homo_sapiens.GRCh37.dna.chromosome.22.fa output_sorted.bam | bcftools call -mv -Ob -o call.bcf

bcftools view -i '%QUAL>=30 && DP>10' call.bcf > call_filtered.bcf

bcftools norm -Ov -f index/Homo_sapiens.GRCh37.dna.chromosome.22.fa -d all -o call_norm.vcf call_filtered.bcf

bcftools stats call_norm.vcf
