# MBDToolBox
MBDToolBox - Docker container to support the data descriptor paper describing an MBD-seq dataset on neuroblastoma samples.

## Reference:
DNA-methylation profiling of primary neuroblastoma tumors using methyl-CpG-binding domain sequencing (MBD-seq)
Anneleen Decock, Maté Ongenaert, Wim Van Criekinge, Frank Speleman, Jo Vandesompele

## Tools:
- sratoolkit (v2.5.4-1) - NCBI toolkit to work with files from SRA (e.g. convert SRA to FASTQ)
- fastqc (v0.11.4) - QC on raw sequencing reads (FASTQ files)
- bowtie2 (v2.2.6) - maps reads to the reference genome (hg19 index files included as well)
- samtools - tools to manipulate SAM files (used to sort/index BAM files)
- picard (v1.140) - tools to manipulate next-generation sequencing formats (BAM files)(used to mark duplicate reads in the BAM file)
- samstat (v1.5.1) - QC statistics on mapped reads (SAM/BAM files)
- BamUtil (v1.0.13) - statistics and operations, QC on mapped reads (BAM files)
- MACS (v1.4.3) - peak caller for ChIP-seq or other capture methods (such as MBD-seq)

## Usage:

This toolbox works with all kinds of enrichment strategies, the example given here is based on a paired-end MBD-seq experiment on a primary neuroblastoma sample

Steps:
1. Get the SRA-file (SRA repository, NCBI) and convert it to FASTQ files (paired-end) (SRAtoolkit)
2. Perform QC on the FASTQ raw reads (FastQC)
3. Map reads to the human hg19 reference genome taking into account the paired-end nature (bowtie2)
4. Sort and index the SAM/BAM files (samtools)
5. Mark duplicates (PCR duplicates during library prep) (picard)
6. Produce statistics and QC on BAM files (samstat, BamUtils)
7. Call peaks to idnetify enriched regions, covered by MBD (MACS)




