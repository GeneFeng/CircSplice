# CircSplice
Pipeline for CircSplice

 

Detection of circ-AS in one sample:

 

1. Quality Control
Librarys of total RNA and rRNA depleted and treated by RNase R are suggested in detecting circRNA and circ-AS. Quality control for the sequencing reads is suggested before mapping to genome.

 

1.1 Run FastQC
fastqc Sample.R1.fq.gz Sample.R2.fq.gz

 

1.2 Remove adapter and low quality bases
trim_galore -q 20 --phred33 --stringency 3 --length 20 -e 0.1 --paired Sample.R1.fq Sample.R2.fq --gzip -o Sample

 

1.3 Check the quality again.
fastqc Sample.R1_trimmed.fq.gz Sample.R2_trimmed.fq.gz

 

2. Map to genome
 

2.1 Generate genome indices
STAR --runThreadN NumberOfThreads --runMode genomeGenerate --genomeDir /path/to/genomeDir --genomeFastaFiles /path/to/genome/fasta1 /path/to/genome/fasta2 --sjdbGTFfile /path/to/annotations.gtf

 

2.2 Mapping
STAR --genomeDir /path/to/genomeDir --readFilesIn Sample.R1_trimmed.fq.gz Sample.R2_trimmed.fq.gz --readFilesCommand zcat --runThreadN 10 --chimSegmentMin 20 --chimScoreMin 1 --alignIntronMax 100000 --outFilterMismatchNmax 4 --alignTranscriptsPerReadNmax 100000 --outFilterMultimapNmax 2 --outFileNamePrefix Sample

 

3. Run CircSplice
CircSplice.pl Chimeric.out.sam hg38.genome.fa bed_refFlat_hg38.txt

 

Input files:

Chimeric.out.sam is in the output of mapping results by STAR.

 

hg38.genome.fa is the genome fasta file.

 

bed-refFlat_hg38.txt is a bed format file generated from UCSC refFlat gene annotation file. User can download it from UCSC genome browser and re-construct themselves. We provide the human hg38 version here. For other species or genome versions, user can re-construct their own file according to the format.

 

The column of bed-refFlat_hg38.txt is as following:

1. Chromosome
2. Start coordinates -2
3. End coordinates +2
4. Transcript ID
5. Number of exons
6. Strand
7. Gene symbol
8. Transcript ID
9. Chromsome
10. Strand
11. Start coordinates
12. end coordinates
13. Start coordinates of CDS
14. End coordinates of CDS
15. Number of exons
16. Start coordinates of each exon
17. End coordinates of each exon
18. Gene type (lncRNA or mRNA)

 

The code for transforming the bed-refFlat format is providing in here. User can run it as:

reftobed.pl refFlat.txt

** Bedtools 2.26 or latest version and Samtools (flagstat) is required to be installed for running CircSplice.

 

Output results:

CircSplice generates two output files: Chimeric.out.sam.result.as and Chimeric.out.sam.result.circ

 

Each row of Chimeric.out.sam.result.as represents one circ-AS events. And each column is explained as following order:

1. The AS type and genomic coordinates. The coordinates are reported according to the position of AS event, which is marked by red line in Figure 1.
2. Number of reads supporting this AS event.
3. The normalized number of reads of this AS event, which is calculated by (Number of reads in this event)/(Number of total chimeric reads)*10^6.
4. Read ID supporting this event.
5. Annotation for this event: gene symbol, transcript and gene type (lncRNA or mRNA).
6. Strand.
7. Chromosome.
8. Genomic coordinates of reads supporting this event.
9. Genomic coordinates of annotated exons related to this event.

 

Each row of Chimeric.out.sam.result.circ represents one circRNA. And each column is explained as following order:

1. Genomic coordinates of donor and acceptor sites of this circRNA.
2. Number of back-splicing reads supporting this circRNA.
3. The normalized number of reads of this circRNA, which is calculated by (Number of reads in this circRNA)/(Number of total chimeric reads)*10^6.
4. Read ID for detecting this circRNA.
5. Annotation for this circRNA: gene symbol, transcript and gene type (lncRNA or mRNA).
6. Strand.
7. Chromosome.

 

Comparison of circ-AS between samples
CircSplice-merge provides merging circRNA or circ-AS results from multiple samples, which let user to compare the circRNA or circ-AS events in different diseases and tissues. CircSplice-merge merges results according to the genomic coordinates with 2bp mismatch toleration.

 

Running CircSplice-merge:
Name the circRNA and circ-AS results from CircSplice as samplename-type-as or samplename-type-circ. For example, s60N-normal-as or s60T-cancer-circ.

Then put circRNA and circ-AS results from different samples separately into two folders named as dir-as and dir-circ. And then run as following:

CircSplice-merge.pl dir-as dir-circ

 

Output files:

This code generates two results: AS_result and Circ_result, which contain the merging results across all input samples.

 

Each row of AS_result represents one merged AS events across all input samples. Each column is explained as following order:

1. The AS type and genomic coordinates. The coordinates is reported according to the position of AS event, which is marked by red line in Figure 1.
2. Genomic coordinates of donor and acceptor sites of the circRNA including this AS event.
3. Sample type
4. Sample name and number of reads supporting this AS event.
5. Sample name and normalized number of reads of this AS event.
6. Annotation for this event: gene symbol, transcript and gene type (lncRNA or mRNA).
7. Strand.
8. Chromosome.
9. Genomic coordinates of reads supporting this event.
10. Genomic coordinates of annotated exons related to this event.

 

Each row of Circ_result represents one merged circRNA across all input samples. And each column is explained as following order:

1. Genomic coordinates of donor and acceptor sites of this circRNA.
2. Sample type.
3. Sample name and number of reads supporting this AS event.
4. Sample name and normalized number of reads of this AS event.
5. Annotation for this event: gene symbol, transcript and gene type (lncRNA or mRNA).
6. Strand.
7. Chromosome.

 

CircSplice-merge keeps all circRNA or circ-AS from all input samples in the result, which let user can select and compare the counts and abundance of different events according to their own requirement.

 

All problems and suggestions could be reported to che@whu.edu.cn or gfeng@whu.edu.cn.
