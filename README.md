# MBDToolBox
MBDToolBox - Docker container to support the data descriptor paper describing an MBD-seq dataset on neuroblastoma samples.

### References:
<b>DNA-methylation profiling of primary neuroblastoma tumors using methyl-CpG-binding domain sequencing (MBD-seq)</b>

Anneleen Decock, Maté Ongenaert, Wim Van Criekinge, Frank Speleman, Jo Vandesompele

Scientific Data, 2016</BR>
doi:10.1038/sdata.2016.4</BR>
http://www.ncbi.nlm.nih.gov/pubmed/26836295</BR>

<b>Methyl-CpG-binding domain sequencing reveals a prognostic methylation signature in neuroblastoma</b>

Anneleen Decock, Maté Ongenaert, Robrecht Cannoodt, Kimberly Verniers, Bram De Wilde, Geneviève Laureys, Nadine Van Roy, Ana P. Berbegall, Julie Bienertova-Vasku, Nick Bown, Nathalie Clément, Valérie Combaret, Michelle Haber, Claire Hoyoux, Jayne Murray, Rosa Noguera, Gaelle Pierron, Gudrun Schleiermacher, Johannes H. Schulte, Ray L. Stallings, Deborah A. Tweddle for the Children’s Cancer and Leukaemia Group (CCLG), Katleen De Preter, Frank Speleman, Jo Vandesompele.

Oncotarget, 2015</BR>
doi: 10.18632/oncotarget.6477</BR>
http://www.ncbi.nlm.nih.gov/pubmed/26646589</BR>

# Tools:
- sratoolkit (v2.5.4-1) - NCBI toolkit to work with files from SRA (e.g. convert SRA to FASTQ)
- fastqc (v0.11.4) - QC on raw sequencing reads (FASTQ files)
- bowtie2 (v2.2.6) - maps reads to the reference genome (hg19 index files included as well)
- samtools (v0.1.19) - tools to manipulate SAM files (used to sort/index BAM files)
- picard (v1.140) - tools to manipulate next-generation sequencing formats (BAM files)(used to mark duplicate reads in the BAM file)
- samstat (v1.5.1) - QC statistics on mapped reads (SAM/BAM files)
- BamUtil (v1.0.13) - statistics and operations, QC on mapped reads (BAM files)
- MACS (v1.4.3) - peak caller for ChIP-seq or other capture methods (such as MBD-seq)

# Usage:

This toolbox works with all kinds of enrichment strategies, the example given here is based on a paired-end MBD-seq experiment on a primary neuroblastoma sample

## Steps:

1. Get the SRA-file (SRA repository, NCBI) and convert it to FASTQ files (paired-end) (SRAtoolkit)
2. Perform QC on the FASTQ raw reads (FastQC)
3. Map reads to the human hg19 reference genome taking into account the paired-end nature (bowtie2)
4. Sort and index the SAM/BAM files (samtools)
5. Mark duplicates (PCR duplicates during library prep) (picard)
6. Produce statistics and QC on BAM files (samstat, BamUtils)
7. Call peaks to identify enriched regions, covered by MBD (MACS)

## Complete script:

**1. Get the SRA-file (SRA repository, NCBI) and convert it to FASTQ files (paired-end) (SRAtoolkit)**
> set some local options

```sudo locale-gen en_US.UTF-8```

> get the SRA file from SRA archives through FTP

```curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX103%2FSRX1039196/SRR2040898/SRR2040898.sra -o 1509_MBD.sra```

```curl ftp://ftp-trace.ncbi.nlm.nih.gov/sra/sra-instant/reads/ByExp/sra/SRX%2FSRX103%2FSRX1039241/SRR2040943/SRR2040943.sra -o 1509_input.sra```

> output: downloaded SRA file

> SRA to FASTQ reads, paired-end (--split-files used) and GZipped (giving two fastq.gz files: one with _1 and the other with _2, representing the paired reads)

```fastq-dump --split-files 1509_MBD.sra --gzip```
```fastq-dump --split-files 1509_input.sra --gzip```

> output: 2 fastq files, ending with _1 and _2, containing the paired reads

**2. Perform QC on the FASTQ raw reads (FastQC)**

> make output directory

```mkdir OUT_fastqc```

> FastQC analysis on one of the pairs, results in the created output directory

```fastqc 1509_MBD_1.fastq.gz --outdir OUT_fastqc/```

> output: in *OUT_fastqc* directory: a HTML file and a ZIP file, containing the QC report and the images

**3. Map reads to the human hg19 reference genome taking into account the paired-end nature (bowtie2)**

> read mapping, qualities from SRA are always phred+33 (--phred33), paired-end with max. 500 bp between pairs (-X 500)
> 8 parallel processes (--threads 8), reference genome hg19 (-x hg19; bowtie2 index files included in container)
> output passed to samtools to create BAM file instead of SAMfile

```bowtie2 -q --phred33 --sensitive -X 500 --threads 8 -t -x bowtie2/hg19 -1 1509_MBD_1.fastq.gz -2 1509_MBD_2.fastq.gz | samtools view -bhSo 1509_MBD_unsorted.bam -```

```bowtie2 -q --phred33 --sensitive -X 500 --threads 8 -t -x bowtie2/hg19 -1 1509_input_1.fastq.gz -2 1509_input_2.fastq.gz | samtools view -bhSo 1509_input_unsorted.bam -```

> output: unsorted BAM file (Bowtie2 output SAM format, which is directly, on the fly converted to BAM using *samtools view* command; don't forget the "-" at the end of the command, it is required for this passing (pipe) to samtools to work

> if you want to, you could remove all FASTQ files and the SRA

```rm 1509*.sra```

```rm 1509*.fastq.gz```

**4. Sort and index the SAM/BAM files (samtools)**

> sort BAM file using samtools

```samtools sort 1509_MBD_unsorted.bam 1509_MBD.sorted```
```samtools sort 1509_input_unsorted.bam 1509_input.sorted```

> output: sorted BAM file, ready for processing with Picard tools

> if you want to, you could remove the unsorted BAM file

```rm 1509_MBD_unsorted.bam```
```rm 1509_input_unsorted.bam```

**5. Mark duplicates (PCR duplicates during library prep) (picard)**

> mark duplicate reads

```java -Xmx8g -jar picard/picard.jar MarkDuplicates INPUT=1509_MBD.sorted.bam OUTPUT=1509_MBD.sorted_nodups.bam ASSUME_SORTED=true METRICS_FILE=Picard_1509_MBD_duplicates.txt```

```java -Xmx8g -jar picard/picard.jar MarkDuplicates INPUT=1509_input.sorted.bam OUTPUT=1509_input.sorted_nodups.bam ASSUME_SORTED=true METRICS_FILE=Picard_1509_input_duplicates.txt```

> output: sorted BAM file, with duplicate reads marked; run summary by Picard (including the read numbers, paired reads, reads marked duplicate)

> if you want to, you could remove the BAM file without duplicates marked

```rm 1509_MBD.sorted.bam```
```rm 1509_input.sorted.bam```

> index sorted BAM file with duplicates flagged

```samtools index 1509_MBD.sorted_nodups.bam```
```samtools index 1509_input.sorted_nodups.bam```

> output: BAM indes, ending with .bam.bai

**6. Produce statistics and QC on BAM files (samstat, BamUtils)**

> perform samstat on sorted BAM file

```samstat 1509_MBD.sorted_nodups.bam```

> output: HTML file with samstat quality statistics

> perform bamUtils on sorted BAM file

```bam stats --in 1509_MBD.sorted_nodups.bam --basic --bamIndex 1509_MBD.sorted_nodups.bam.bai```

```bam stats --in 1509_MBD.sorted_nodups.bam --phred --bamIndex 1509_MBD.sorted_nodups.bam.bai```

> output: BamUtil statistics (outputted on screen)

**7. Call peaks to identify enriched regions, covered by MBD (MACS)**

> make output directory

```mkdir OUT_macs```

> call peaks and write a single WIG file for visualization in IGV etc. (--wig --single-profile)

```macs -t 1509_MBD.sorted_nodups.bam -c 1509_input.sorted_nodups.bam --outdir=OUT_macs/ --name=1509_MBD_macspeaks -f BAM --petdist=200 -g hs --wig --single-profile ```

> output: MACS peak files (Excel-file and BED files); WIG files for visualization, all in the *OUT_macs* directory


# Alternative scripts & visualization

## Alternative scripts

Other included scripts are examples of further analysis.

> Example how to run BALM - also 'peak / enriched regions' based

```script_BALM.sh ```

> R scripts and utilities to create a count matrix around transcription start sites (TSS) or in 5 kb genomic windows

```
script_count2000.R
TSS2000.bed
```

> R script demonstrating MEDIPS, an R/BioConductor package, calculating absolute methylation scores for all CG sites

```script_MEDIPS.R ```

## Visualization

**1. Get the TDF and BED.gz files from the GEO repositories (Cohort I:GSE69224; Cohort II: GSE69243)**
</BR>These files are in the "Supplementary file" section of GEO; you can download all these files at the same time.
</BR>For the supplied IGV session example file, you need the files for cohort II, which can be all downloaded from this location: http://www.ncbi.nlm.nih.gov/geo/download/?acc=GSE69243&format=file
</BR>Caution: file size: 32 GB

**2. Decompress the combined archive if needed (if downloading multiple files) and place them in a directory with the supplied IGV session file (example for cohort II, file:"MBD cohort II.xml") and the sample annotations file ("annotations.txt")**
</BR>It is crucial that the files are in the same directory as the session-file as the paths are defined relative

**3. Download and open IGV (https://www.broadinstitute.org/igv/) and load the session file.**
</BR>Alternatively, individual tracks (TDF files and BED files) can be loaded into IGV.
