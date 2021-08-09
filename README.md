# GBS introgression line analysis

This repository contains source code for analysis of introgression lines. 

This is in context of a soon to be submitted manuscript (todo: update when paper is online)


## Context

Identify position of introgression in material based on GBS sequencing data from donor, background and introgression line

## Preprocessing

### 1. Mapping

Map GBS sequencing reads to a reference genome. At time this analysis was done, large chromosomes could not be directly processed. They needed to be split into parts. The scripts in this repository expect a BED file detailing the split positions.

We used [bwa](http://bio-bwa.sourceforge.net/) for mapping and processed mappings with [samtools] (http://www.htslib.org/).

Process every individual gbs data set (including donor, introgression lines and background) as below.

```
bwa index reference_parts.fasta
samtools faidx reference_parts.fasta
bwa mem reference_parts.fasta gbs.fastq > alignment.sam
samtools sort -o alignment.bam alignment.sam
samtools index alignment.bam
samtools mpileup -BQ0 -f reference_parts.fasta alignment.bam > alignment.pileup

```

### 2. Integrate donor, background and introgression line variations

The java program will load every position in the reference that has coverage from background line mapping. It records positions where background allele differs from reference allele. These position will be excluded from further processing. Next, positions are loaded where donor has a different allele than reference. Only positions are regarded that have coverage from background. Finally, the introgression line data is loaded. For each mega-base of the reference, the donor alleles present in the introgression line are counted.

Executing the scripts requires [jre 1.6 or higher](https://www.oracle.com/uk/java/technologies/javase-java-archive-javase6-downloads.html)


```
java -jar GBSSNPs.jar -b <background.pileup> -d <donor.pileup> -i <introgression.pileup> -o <output1.txt> -p <parts2chrs.bed> -l <interval_length>
```

An example for the parts2chrs.bed file is [here](https://github.com/steuernb/GBS_introgression_line_analysis/blob/main/161010_Chinese_Spring_v1.0_pseudomolecules_parts_to_chr.bed). For our analysis, we used an interval_length of 1000000 (1 mega-base). 


### 3. Plot output

The output file from step 2, <output1.txt>, can be used to plot percentages per MB. To accurately plot chromosome lengths, a fasta index of the (not split) reference file is used: `samtools faidx <reference.fasta>` will create an "fai"-file used in script below. 

```{R}
lengths <- read.table("genome.fasta.fai")

chrs <- c("chr1A", "chr1B", "chr1D",
          "chr2A", "chr2B", "chr2D",
          "chr3A", "chr3B", "chr3D",
          "chr4A", "chr4B", "chr4D",
          "chr5A", "chr5B", "chr5D",
          "chr6A", "chr6B", "chr6D",
          "chr7A", "chr7B", "chr7D")



file <- read.table("output1.txt" , sep = "\t", stringsAsFactors = FALSE, header = FALSE)
par(mfrow=c(7,3) ,mar = c(3, 4, 4, 2) + 0.1, mgp = c(2,1,0))
x<- file$V2
y<- (file$V4)*100/(file$V3+0.0000001)
for( i in chrs){
  plot(    x[file$V1==i],y[file$V1==i],
              type = "h", xlab = i,  ylab ="", 
              ylim = c(0, max(y)), 
              xlim = c(0,lengths$V2[lengths$V1==i])
      )
}

```


