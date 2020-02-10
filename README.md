
# Towards Unravelling Transcriptomic Changes In Human Airway Smooth Muscles 

## Abstract: 

Asthma used to be a chronic inflammatory airway disease which can be treated by and glucocorticosteroids and β2-agonists. Airway smooth muscle is one of the targeted primary tissues for treatment . An Experiment was performed for retrieving RNA-Seq in order to characterize the human airway smooth muscle (HASM) transcriptome at baseline and under three asthma treatment conditions [2]. In our analysis, Two conditions only were used Differential expression analysis for Airway smooth muscle tissues first at baseline and secondly at treatment with glucocorticoids will be performed to identify the differentially expressed genes in order to characterize transcriptomic changes in human ASM cell lines that were treated with dexamethasone a potent synthetic glucocorticoid and without treatment.

# Methods and Materials:

## Datasets and pre-processing: 
The analysis adopted a data set obtained from The Illumina TruSeq assay for project (PRJNA229998) under SRA (SRP033351) [3]. Short read data of Airway smooth muscle human samples were collected from the NCBI  database. There were four samples in total, these samples were classified into two conditions, each has three Runs: an untreated samples (at baseline) with Run accession number; (SRR1039512,SRR1039520,SRR1039516) and treated samples with Dexamethasone as Glucocorticoid; (SRR1039513,SRR1039517,SRR1039521). NCBI’s fastq-dump from sra-toolkit was used to download the short reads for NCBI short read archive (SRA) 
The reads, in a single file, were paired-ends so it needs to be split for the downstream analysis

The following command can download the data directly from NCBI, without prior prefetch command, split it into two reads and subset the reads according to our parameters. Troubleshooting downloading large dataset size is explained in detalis at the following issue 
[Issue 1](https://github.com/hagarelsayed/Ngs_2nd_Abstract/issues/1)
```
cd ~/workdir/sample_data
fastq-dump --outdir fastq --gzip --skip-technical --readids --read-filter pass --dumpbase --split-3 -N 10000 -X 4010000 --clip SRR1039520
# all samples downloaded by the same command 

SRR1039508 , SRR1039509 , SRR1039512 ,SRR1039513 ,SRR1039516, SRR1039517 , SRR1039520 , SRR1039521
```
–split-3 was used to separate the read into left and right ends each in a single file R1& R2, the third file will be for left ends without a matching right end and for a right end without a matching left ends  they will be put in a single file.

--skip-technical to remove technical reads produced from the “Illumina multiplex library construction protocol” 

–clip to Apply left and right clips to remove sequences that include tags which were used for amplification.

 -N 10000 -X 4010000 It was used to subset the data to 4Million reads, It was started from position 10000 in order to skip the unusual sequences at the beginning of the fastq files.

Previously the data was downloaded like this, which consumes lots of time as it caused elongation of time needed to download the whole set first then subset it .

```
 prefetch SRR1039508
 fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039508
 prefetch SRR1039509
 fastq-dump --outdir fastq --gzip --skip-technical  --readids --read-filter pass --dumpbase --split-3 --clip SRR1039509
```
Other previous attempts for downloading the data was explained in details via [Issue1](https://github.com/hagarelsayed/Ngs_2nd_Abstract/issues/1
)
# Reference Genome

There were multiple Reference genomes that have been used , Each at a time for the troubleshooting matters and to enhance the alignment rates. All the references were downloaded from NCBI Assemly, Ensembl and Gencodes websites, respectively. (Homo_sapiens.GRCh37.75.gtf.gz, Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz, Homo_sapiens.GRCh38.86.gtf.gz, Homo_sapiens.GRCh38.dna.chromosome.11.fa, gencode.v33.annotation.gtf, gencode.v33.transcripts.fa. 


```
wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
wget https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
gunzip -k Homosapiens_GRCh37.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.chr.gtf.gz

```

Downloading one chromosome only to shorten the time needed for the other steps;

```
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
gunzip -k Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz 

wc Homo_sapiens.GRCh38.dna.chromosome.11.fa
head -n 425000 Homo_sapiens.GRCh38.dna.chromosome.11.fa | tail
cat Homo_sapiens.GRCh38.dna.chromosome.11.fa |  grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'
```

## Quality control of reads was performed vi FastQC tool

Results: the quality for the 12 reads for the 6 runs can be obtained in details from this link https://github.com/Saraelattar/fastqc/issues/2 [](https://github.com/Saraelattar/fastqc/issues/2
) . 
The following plot indicates a good quality of the samples as it all lies under 30%.

![per_base_quality](https://user-images.githubusercontent.com/60422836/74183300-5f1ab500-4c4d-11ea-8cdf-0b5ba62dd582.png)

The percentage of sequence remaining if deduplicated  is 66.6% which indicates that 66.6% of the reads are duplicated between 10-50 times as shown in the following plot
![duplication_levels](https://user-images.githubusercontent.com/60422836/74184271-2a0f6200-4c4f-11ea-91f2-5a6b98ba1a53.png)

The N content is almost zero across all bases which is good as te read were subsetted from the spot 10000 to skip the ambiguity of the Ns that is usually present at the begining of the fastq files. 
![per_base_n_content](https://user-images.githubusercontent.com/60422836/74184954-65f6f700-4c50-11ea-9293-e7f72736879a.png)

The Sequence content across all bases is quite good as it shows that. The proportion of each base position in the reads are with equal proportion .
![per_base_sequence_content](https://user-images.githubusercontent.com/60422836/74185175-ca19bb00-4c50-11ea-979f-dc8d18f1e8b5.png)

The GC content across the whole read length is very good as it resembles the modelled normal distribution of GC content as shown below.
![per_sequence_gc_content](https://user-images.githubusercontent.com/60422836/74185656-bf135a80-4c51-11ea-835a-939610b8291f.png)



## Indexing by Hisat2
In order to do the RNA-Seq analysis , we need to align the data against the transcriptome. 

```
INDEX=./Homo_sapiens.GRCh37
REF= ./Homo_sapiens.GRCh37.75.fa

hisat2-build $REF $INDEX
```

## Alignment by Hisat2: 
```
INDEX=~/Diff_proj/index/Homo_sapiens.GRCh37.75.gtf
RUNLOG=runlog.txt
READS_DIR=~/Diff_proj/sample_data/
mkdir bam

### Names of the samples were edited manually to begin with TTT for treated samples with Glucocorticoid and UNT for untreated samples 

for SAMPLE in UNT TTT ; do     for REPLICATE in 08 09 12 13 16 17 20 21;     do         R1=$READS_DIR/${SAMPLE}_SRR5${REPLICATE}*r1.fastq;         R2=$READS_DIR/${SAMPLE}_SRR5${REPLICATE}*r2.fastq;         BAM=bam/${SAMPLE}_${REPLICATE}.bam;          hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM;         samtools index $BAM;     done; done
```
Alignment rate was very low, around 3%

## Re_Indexing and Re-Alignment by Hisat2 
 The alignment was redone against transcriptome instead of the whole genome.


```
INDEX=./gencode.v33.transcripts
REF= ./gencode.v33.transcripts.fa

hisat2-build $REF $INDEX
```
```

INDEX=~/Diff_proj/index/gencode.v33.transcripts
RUNLOG=runlog.txt
READS_DIR=~/Downloads/fastqq/fastq/fastg/
mkdir bam


for SAMPLE in UNT;
do
    for REPLICATE in 12 16 20;
    do
        R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_1.fastq.gz
        R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_2.fastq.gz
        BAM=bam/${SAMPLE}_${REPLICATE}.bam

        hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM
        samtools index $BAM
    done
done

for SAMPLE in TTT;
       do 
         for REPLICATE in 13 17 21; 
         do 
           R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_1.fastq.gz  
           R2= $READS_DIR/${SAMPLE}Rep${REPLICATE}*pass_2.fastq.gz   BAM=bam/${SAMPLE}${REPLICATE}.bam
hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM
        samtools index $BAM
        done
done
```
Alignment rate was around  86% 

Alignment Summary;
```
4000001 reads; of these:
  4000001 (100.00%) were paired; of these: 
     656023(16.40%) aligned concordantly 0 times
     894885 (22.37%) aligned concordantly exactly 1 time 
     2449093 (61.23%) aligned concordantly >1 times
       ---- 
      656023 pairs aligned concordantly 0 times; of these: 
        6338 (0.97%) aligned discordantly 1 time
       ---- 
      649685 pairs aligned 0 times concordantly or discordantly; of these: 
            1299370 mates make up the pairs; of these:
            1107153 (85.21%) aligned 0 times 
            61740 (4.75%) aligned exactly 1 time
            130477 (10.04%) aligned >1 times
 86.16% overall alignment rate

```

## Quantification:

```
GTF=~/Diff_proj/gencode.v33.transcripts.gtf 
featureCounts -a $GTF -g gene_name -o counts.txt  Bam_org/UNT*.bam  Bam_org/TTT*.bam
```
The following error popped up; 

` The chromosome name of "ENST00000456328.2|ENSG00000223972.5|OTTHUMG00000000961.2|OTTHUMT00000362751.1|DDX11L1-202|DDX11L1|1657|processed_transcript|" contains 125 characters, longer than the upper limit of 99 featureCounts has to stop running`

Troublshooting of this error explained in details in [Issue 5](https://github.com/hagarelsayed/Ngs_2nd_Abstract/issues/5
)
and  [issue 3](https://github.com/hagarelsayed/Ngs_2nd_Abstract/issues/3
)
An even smaller subset was retrieved form the fastq data for testing and checking the other workflow.

Subset small number ; For test only for other workflow ;
`for file in ./*.fastq.gz ; do     echo $file ;     seqtk sample -s100 $file 500 > ${file/.fastq.gz/.fastq}; done`

## Alignment
As previous alignment trials was consuming too much time and computation power so It was an idea to work on the ERCC that was gotten from the course material to work on to check for the feature count error. 
```
INDEX=chr22_with_ERCC92
RUNLOG=runlog.txt
READS_DIR=~/workdir/sample_data/renamed/
mkdir bam


for SAMPLE in UNT;
do
    for REPLICATE in 12 16 20;
    do
        R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_1.fastq
        R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_2.fastq
        BAM=bam/${SAMPLE}_${REPLICATE}.bam

        hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM
        samtools index $BAM
    done
done

```
### Results of the first set; 
This is one result summary , the remaining is all included in the issue 
500 reads; of these:
  500 (100.00%) were paired; of these:
    478 (95.60%) aligned concordantly 0 times
    22 (4.40%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    478 pairs aligned concordantly 0 times; of these:
      0 (0.00%) aligned discordantly 1 time
    ----
    478 pairs aligned 0 times concordantly or discordantly; of these:
      956 mates make up the pairs; of these:
        951 (99.48%) aligned 0 times
        5 (0.52%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
4.90% overall alignment rate


![dd](https://user-images.githubusercontent.com/60422836/74152521-265fe900-4c17-11ea-9683-e6bfcb2aa090.png)

Code for second set;

```
for SAMPLE in TTT;
do
    for REPLICATE in 13 17 21;
    do
        R1=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_1.fastq
        R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_2.fastq
        BAM=bam/${SAMPLE}_${REPLICATE}.bam

        hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM
        samtools index $BAM
    done
done
```



## Results of second set : 

500 reads; of these:
  500 (100.00%) were paired; of these:
    486 (97.20%) aligned concordantly 0 times
    14 (2.80%) aligned concordantly exactly 1 time
    0 (0.00%) aligned concordantly >1 times
    ----
    486 pairs aligned concordantly 0 times; of these:
      0 (0.00%) aligned discordantly 1 time
    ----
    486 pairs aligned 0 times concordantly or discordantly; of these:
      972 mates make up the pairs; of these:
        967 (99.49%) aligned 0 times
        5 (0.51%) aligned exactly 1 time
        0 (0.00%) aligned >1 times
3.30% overall alignment rate


## Quantification

```
GTF=~/workdir/diff_exp/ref/ERCC92.gtf 
featureCounts -a $GTF -g gene_name -o counts.txt  bam/UNT*.bam  bam/TTT*.bam
```

The following error came up : 

`Failed to open the annotation file /home/ngs/workdir/diff_exp/ref/ERCC92.gtf, or its format is incorrect, or it contains no 'exon' features
`

![crop](https://user-images.githubusercontent.com/60422836/74151840-b309a780-4c15-11ea-9e28-c5789f79fffb.png)

The Reference genome changed to 
`GTF=~/workdir/sample_data/gencode.v29.annotation.gtf
`
The feature count worked smoothly and this is the out put results 
![22](https://user-images.githubusercontent.com/60422836/74152073-31fee000-4c16-11ea-928c-43f2d1997710.png)


```
cat counts.txt | cut -f 1,7-12 > simple_counts.txt 
less simple_counts.txt
```

## Results of Quantification:

The results could not be uploaded to git but found at this link 

[Results of Quantification](https://drive.google.com/open?id=14K528G5B0N6N54OKRwXT5pZlUbBnYUWU
)


## Aligning to our original index (not a subset like in the past example) :
There was an attempt  towards changing the reference to our original reference that gave high alignment rate with the dataset so that the feature count will result in a count matrix containing different versatile values for deseq to work on.

```
INDEX=gencode.v33.transcripts

###Other variables are the same and defined previously 
RUNLOG=runlog.txt
READS_DIR=~/workdir/sample_data/renamed/
mkdir bam
 
for SAMPLE in UNT; do     for REPLICATE in 12 16 20;     do         R1=$READS
_DIR/${SAMPLE}_Rep${REPLICATE}*pass_1.fastq;         R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_2.fastq;         BAM=bam/${SAMPLE}_${REPLICATE
}.bam;          hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM;         samtools index $BAM;     done; done
```
Results of the Alignment score for the untreated condition
500 reads; of these:
  500 (100.00%) were paired; of these:
    71 (14.20%) aligned concordantly 0 times
    108 (21.60%) aligned concordantly exactly 1 time
    321 (64.20%) aligned concordantly >1 times
    ----
    71 pairs aligned concordantly 0 times; of these:
      4 (5.63%) aligned discordantly 1 time
    ----
    67 pairs aligned 0 times concordantly or discordantly; of these:
      134 mates make up the pairs; of these:
        112 (83.58%) aligned 0 times
        9 (6.72%) aligned exactly 1 time
        13 (9.70%) aligned >1 times
88.80% overall alignment rate


Same step for the other group
```
for SAMPLE in UNT; do     for REPLICATE in 12 16 20;     do         R1=$READS
_DIR/${SAMPLE}_Rep${REPLICATE}*pass_1.fastq;         R2=$READS_DIR/${SAMPLE}_Rep${REPLICATE}*pass_2.fastq;         BAM=bam/${SAMPLE}_${REPLICATE
}.bam;          hisat2 $INDEX -1 $R1 -2 $R2 | samtools sort > $BAM;         samtools index $BAM;     done; done
```
Alignment summary results :+1: 
500 reads; of these:
  500 (100.00%) were paired; of these:
    87 (17.40%) aligned concordantly 0 times
    111 (22.20%) aligned concordantly exactly 1 time
    302 (60.40%) aligned concordantly >1 times
    ----
    87 pairs aligned concordantly 0 times; of these:
      0 (0.00%) aligned discordantly 1 time
    ----
    87 pairs aligned 0 times concordantly or discordantly; of these:
      174 mates make up the pairs; of these:
        143 (82.18%) aligned 0 times
        6 (3.45%) aligned exactly 1 time
        25 (14.37%) aligned >1 times
85.70% overall alignment rate

![err2](https://user-images.githubusercontent.com/60422836/74164856-00ddda00-4c2d-11ea-97a0-8548c55c1d0c.png)


## Quantification

```
GTF=~/workdir/hisat_align/hisatIndex/proj_index/gencode.v33.annotation.gtf
featureCounts -a $GTF -g gene_name -o counts.txt  bam/UNT*.bam  bam/TTT*.bam
```
Result of Quantification 


Differential Expression:

`cat simple_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.tsv
`
Another Error popped up ;
![33](https://user-images.githubusercontent.com/60422836/74167861-924f4b00-4c31-11ea-99cd-0d34080e36ac.png)

The table of the count is availabe at this link 

https://github.com/hagarelsayed/Ngs_2nd_Abstract/issues/6




# Differential Expression 


`cat counts.txt | cut -f 1,7-12 > simple_counts.txt`

Simple counts produced 
[simple_counts.txt](https://github.com/hagarelsayed/Ngs_2nd_Abstract/files/4181182/simple_counts.txt)

All the matrix is Zero may be because of the low alignment rate for the small genome

while trying to do the differential analysis by Deseq, 
`cat simple_counts.txt | Rscript deseq1.r 3x3 > results_deseq1.tsv
`
The following error came up :
Error in parametricDispersionFit(means, disps) :
Parametric dispersion fit failed. Try a local fit and/or a pooled estimation.

![Deseq Error](https://user-images.githubusercontent.com/60422836/74156848-d9ccdb80-4c1f-11ea-9147-c0582d770393.png)

The error may be due to values on the matrix is zero which was a reult of aligning to a small portion of the reference genome 
Now will try to get back to the alignment step to index another genome
Details in the following link [issue 6](https://github.com/hagarelsayed/Ngs_2nd_Abstract/issues/6)

## Quality checking for the fastq files after processing and calculating GC content

Details of Fastqc files at the following issue [](https://github.com/Saraelattar/fastqc/issues/1)

convert bam file to fastq file

```
samtools bam2fq TTT_13.bam > TTT_13.fastq
 samtools bam2fq TTT_17.bam > TTT_17.fastq
 samtools bam2fq TTT_21.bam > TTT_21.fastq
 
samtools bam2fq UNT_12.bam > UNT_12.fastq
samtools bam2fq UNT_20.bam > UNT_20.fastq
 samtools bam2fq UNT_16.bam > UNT_16.fastq
```
![image](https://user-images.githubusercontent.com/60422836/74192075-7a8dbc00-4c5d-11ea-9f38-0901d507e659.png)

