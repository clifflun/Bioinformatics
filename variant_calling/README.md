![svc](flowchart.png)

This is a very brief outline of the pipeline for short variant calling.

First, I got one subject SRA data from GEO/SRA database. 

FASTQ data was obtained using sratool-kit prefetch and fasterq-dump/fastq-dump.

FASTQC is used to determine if trimming is needed.

Bowtie2 was used for alignment in this case. The sequence was aligned to chr22 of hg19, which is downlaoded from Ensembl.

Since this is somatic cells, bcftools is used. 

For germline variant calling, GATK pipelines could be used. (Haplotypecaller)

VCF normalization and filtering was done to filter out lower quality variant calls.

We can then visualize the vcf/bcf file in a genome viewer like IGV..

![igv](igv.png)

