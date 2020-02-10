
# Towards Unravelling Transcriptomic Changes In Human Airway Smooth Muscles 

## Abstract: 

Asthma used to be a chronic inflammatory airway disease which can be treated by and glucocorticosteroids and β2-agonists. Airway smooth muscle is one of the targeted primary tissues for treatment . An Experiment was performed for retrieving RNA-Seq in order to characterize the human airway smooth muscle (HASM) transcriptome at baseline and under three asthma treatment conditions [2]. In our analysis, Two conditions only were used Differential expression analysis for Airway smooth muscle tissues first at baseline and secondly at treatment with glucocorticoids will be performed to identify the differentially expressed genes in order to characterize transcriptomic changes in human ASM cell lines that were treated with dexamethasone a potent synthetic glucocorticoid and without treatment.

# Methods and Materials:

## Datasets and pre-processing: 
The analysis adopted a data set obtained from The Illumina TruSeq assay for project (PRJNA229998) under SRA (SRP033351) [3]. The Experiment used high through-put (to be completed)


NCBI’s fastq-dump from sra-toolkit was used to download the short reads for NCBI short read archive (SRA).
The reads, in a single file, were paired-ends so it needs to be split for the downstream analysis.
The following command can download the data directly from NCBI, without prior prefetch command, split it into two reads and subset the reads according to our parameters. 
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

Getting the Reference Genome in gtf Format :



```
wget http://ftp.ensembl.org/pub/release-75/gtf/homo_sapiens/Homo_sapiens.GRCh37.75.gtf.gz
wget https://www.ncbi.nlm.nih.gov/assembly/GCF_000001405.13
gunzip -k Homosapiens_GRCh37.fa.gz 
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.primary_assembly.fa.gz
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.gtf.gz
wget ftp://ftp.ensembl.org/pub/release-86/gtf/homo_sapiens/Homo_sapiens.GRCh38.86.chr.gtf.gz

```

Downloading one chromosome only to shorten the time needed for the other steps , 

```
wget ftp://ftp.ensembl.org/pub/release-86/fasta/homo_sapiens/dna/Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz
gunzip -k Homo_sapiens.GRCh38.dna.chromosome.11.fa.gz 

wc Homo_sapiens.GRCh38.dna.chromosome.11.fa
head -n 425000 Homo_sapiens.GRCh38.dna.chromosome.11.fa | tail
cat Homo_sapiens.GRCh38.dna.chromosome.11.fa |  grep -v ">" | perl -ne 'chomp $_; $bases{$_}++ for split //; if (eof){print "$_ $bases{$_}\n" for sort keys %bases}'
```


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

Troublshooting of this error explained in detalis in [Issue 5](https://github.com/hagarelsayed/Ngs_2nd_Abstract/issues/5
)
and  [issue 3](https://github.com/hagarelsayed/Ngs_2nd_Abstract/issues/3
)
It was solved in issue 5

