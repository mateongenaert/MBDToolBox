# MBDToolBox
MBDToolBox - Docker container to support the data descriptor paper describing an MBD-seq dataset on neuroblastoma samples.

## Reference:
DNA-methylation profiling of primary neuroblastoma tumors using methyl-CpG-binding domain sequencing (MBD-seq)
Anneleen Decock, Maté Ongenaert, Wim Van Criekinge, Frank Speleman, Jo Vandesompele

## Tools:
- sratoolkit (v2.5.4-1) - NCBI toolkit to work with files from SRA (e.g. convert SRA to FASTQ)
- fastqc (v0.11.4) - QC on raw sequencing reads (FASTQ files)
- bowtie2 (v2.2.6) - maps reads to the reference genome (hg19 index files included as well)
- samtools (v0.1.19) - tools to manipulate SAM files (used to sort/index BAM files)
- picard (v1.140) - tools to manipulate next-generation sequencing formats (BAM files)(used to mark duplicate reads in the BAM file)
- samstat (v1.5.1) - QC statistics on mapped reads (SAM/BAM files)
- BamUtil (v1.0.13) - statistics and operations, QC on mapped reads (BAM files)
- MACS (v1.4.3) - peak caller for ChIP-seq or other capture methods (such as MBD-seq)

## Usage:

This toolbox works with all kinds of enrichment strategies, the example given here is based on a paired-end MBD-seq experiment on a primary neuroblastoma sample

### Steps:
1. Get the SRA-file (SRA repository, NCBI) and convert it to FASTQ files (paired-end) (SRAtoolkit)</BR>
2. Perform QC on the FASTQ raw reads (FastQC)</BR>
3. Map reads to the human hg19 reference genome taking into account the paired-end nature (bowtie2)</BR>
4. Sort and index the SAM/BAM files (samtools)</BR>
5. Mark duplicates (PCR duplicates during library prep) (picard)</BR>
6. Produce statistics and QC on BAM files (samstat, BamUtils)</BR>
7. Call peaks to identify enriched regions, covered by MBD (MACS)</BR>

### Complete script:

<b>1. Get the SRA-file (SRA repository, NCBI) and convert it to FASTQ files (paired-end) (SRAtoolkit)</b></BR>

curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX209%2FSRX209021/SRR629532/SRR629532.sra -o SRR629532.sra
fastq-dump --split-files SRR629532.sra --gzip

<b>2. Perform QC on the FASTQ raw reads (FastQC)</b></BR>

mkdir OUT_fastqc
fastqc SRR629532_1.fastq.gz --outdir OUT_fastqc/

<b>3. Map reads to the human hg19 reference genome taking into account the paired-end nature (bowtie2)</b></BR>

bowtie2 -q --phred33 --very-sensitive -X 500 --threads 8 -t -x hg19 -1 SRR629532_1.fastq.gz -2 SRR629532_2.fastq.gz | samtools view -bhSo SRR629532_unsorted.bam

<b>4. Sort and index the SAM/BAM files (samtools)</b></BR>

samtools sort -o SRR629532.sorted.bam SRR629532_unsorted.bam

<b>5. Mark duplicates (PCR duplicates during library prep) (picard)</b></BR>

java -jar picard/picard.jar MarkDuplicates INPUT=SRR629532.sorted.bam OUTPUT=SRR629532.sorted_nodups.bam ASSUME_SORTED=true METRICS_FILE=Picard_SRR629532_duplicates.txt
samtools index SRR629532.sorted_nodups.bam

<b>6. Produce statistics and QC on BAM files (samstat, BamUtils)</b></BR>

mkdir OUT_samstat

samstat SRR629532.sorted_nodups.bam > OUT_samstat/SRR629532_samstats.txt

mkdir OUT_bamUtils

bam stats --in SRR629532.sorted_nodups.bam --basic --bamIndex SRR629532.sorted_nodups.bam.bai --pBaseQC OUT_bamUtils/SRR629532_pbaseQC.txt > OUT_bamUtils/SRR629532_stats_basic.txt
bam stats --in SRR629532.sorted_nodups.bam --phred --bamIndex SRR629532.sorted_nodups.bam.bai > OUT_bamUtils/SRR629532_stats_phred.txt

<b>7. Call peaks to identify enriched regions, covered by MBD (MACS)</b></BR>

mkdir OUT_macs

macs -t SRR629532.sorted_nodups.bam --outdir=OUT_macs/ --name=SRR629532_macspeaks -f BAM --petdist=200 -g hs --wig --single-profile 




