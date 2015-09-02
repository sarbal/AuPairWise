#################################################################################
#################################################################################
#################################################################################

The sample_brainspan.Rdata contains the RNA-seq expression estimates for the 578 samples from brainspan.
There are 3 variables in the file.
- exprs
- rows
- cols



The sample_ENCODE.Rdata contains the RNA-seq expression estimates for the 42 samples from the ENCODE dataset (GSE ).
There are three expression estimates, one from RSEQtools (RPKM), one from Cufflinks (FPKM), and one from HTseq (counts).
See below for more details and parameters.
Each row is a gene, labelled by their human entrez gene ID.
Each column is a samplem labeled by the SRR ID (SRA run).

There are 4 variables in the file.
- exprs.rseq
- exprs.htseq
- exprs.cuff
- cols

The SRA files for each samples were downloaded and converted to FASTQ files.
Each sample was QC'd using the FASTX toolkit. Bowtie2 was used to align the reads to the genome.

# Versions and filenames
sratoolkit.2.3.2-5-centos_linux64
gencode.v3c.annotation.GRCh37.gtf
cufflinks-2.1.1
hg19
# from rseqtools
knownGene_2x20_spliceJunctions
knownGene_2x70_spliceJunctions
knownGene_composite.interval
knownToEnsemblToEntrez

# FASTX QC
fastx_artifacts_filter -v -Q 33 -i $fastq -o $file.fastq.artifacts 2>&1
fastq_quality_filter -v -Q 33 -q $Q -p $P -i $file.fastq.artifacts -o $file.fastq.filtered 2>&1
fastx_quality_stats -Q 33 -i $file.fastq.filtered -o $file.stats.txt
fastq_quality_boxplot_graph.sh -i $file.stats.txt -o $file.quality.png

# Bowtie2 alignments
bowtie2 -N 1 -p 16 -3 $t hg19 $file.fastq.filtered -S $file.sam
samtools view -bS  $file.sam | samtools sort - $file.sorted

# RSEQtools
sam2mrf < $file.sam >  $file.bowtie2.mrf

# Cufflinks
cufflinks -p 16 $file.sorted.bam -G gencode.v3c.annotation.GRCh37.gtf -o $out

# HTseq
htseq-count $file.sam gencode.v3c.annotation.GRCh37.gtf > $file.htseq.out

