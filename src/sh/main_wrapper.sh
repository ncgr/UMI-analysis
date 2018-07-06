#!/bin/bash

#
# Using bioconda to manage bioinformatics tools. This makes sure they
# are in the path (e.g. STAR)
#
# Do this first:
# source activate your_conda_enviroment

#==================================================================
# YOU MUST edit $CBIN/umi_analyzer_uthash/exe.c to make the width of the 
# chromosome position variable and the width of the UMI match your data, and
# recompile

#==================================================================
# EDIT HERE
LOCALBIN="/your/local/build/dir" # If you want to hack scripts, put them here and use instead of PERLBIN or CBIN
UMI_EDIT_DISTANCE=1              # UMI differ by <= this, collapse them
MAX_SHIFT=3                      # The biggest mapping shift for adjacent posiions to be in the same cluster
UMI_LENGTH=10                    # 
GGG=3                            # How many GGG to remove after UMI
QUAL_LIMIT=30                    # No reads with Q less than this in the UMI
OFFSET=0                         # Where to begin the UMI substring (UMI)

#==================================================================
# CHOOSE GFF FILE AND STAR INDEX FOLDER
# files. Edit these for your choice of organism and genome
star_genomedir=/path/to/physco/genome/ref_3.3/star_index
gff_file=/path/to/physco/genome/ref_3.3/Ppatens_318_v3.3.gene_exons.gff 

# Locations of executables. In the Git project dir
PERLBIN=/path/to/umi_pipeline/src/perl
CBIN=/path/to/umi_pipeline/src/c

# The executables
filter_exe=$CBIN/fastq_qual_filter/fastq_qual_filter
clipper_exe=$CBIN/fastq_umi_clipper/fastq_umi_clipper
counterizer_exe=$PERLBIN/umi_counterizer.pl
analyzer_exe=$CBIN/umi_analyzer_uthash/exe
clusterizer_exe=$PERLBIN/umi_clusterizer.pl
gene_umi_counter_exe=$PERLBIN/gene_umi_counter.pl
sam_to_bed_exe=$PERLBIN/sam_to_bed.pl
gene_read_counter_exe=$PERLBIN/gene_read_counter.pl

# comment things out to run one at a time
$filter_exe seq.fq good.fq bad.fq log $QUAL_LIMIT $OFFSET $UMI_LENGTH
$clipper_exe good.fq $UMI_LENGTH $GGG > clip.fq 

# If you use an aligner other than STAR, it needs to specify unique aligners with this: "NH:i:1" or next step will break
STAR  --runThreadN 16 --readFilesIn clip.fq --genomeDir $star_genomedir --outFilterMismatchNmax 2 --alignIntronMax 1 --alignEndsType EndToEnd --outSAMattributes All
cat Aligned.out.sam | pcregrep '^@' > tmp.sam
cat Aligned.out.sam | pcregrep '\tNH:i:1\t' >> tmp.sam
mv tmp.sam Aligned.out.uniq.sam
samtools view -b Aligned.out.uniq.sam -o Aligned.out.uniq.bam
samtools sort -o Aligned.sorted.out.uniq.bam Aligned.out.uniq.bam
samtools view -h Aligned.sorted.out.uniq.bam > Aligned.sorted.out.uniq.sam
$counterizer_exe Aligned.sorted.out.uniq.sam $UMI_EDIT_DISTANCE > umi_counts.txt
$analyzer_exe  umi_counts.txt  > umi_positions.txt
$clusterizer_exe umi_positions.txt $MAX_SHIFT > position_umi_counts.bed
bedtools intersect -a position_umi_counts.bed -b $gff_file -wa -wb | pcregrep '\tgene\t' > umi_bed_intersection_gff.txt
$gene_umi_counter_exe umi_bed_intersection_gff.txt > gene_umi_counts.txt
$sam_to_bed_exe Aligned.sorted.out.uniq.sam > position_read_counts.bed
bedtools intersect -a position_read_counts.bed -b $gff_file -wa -wb | pcregrep '\tgene\t' > read_bed_intersection_gff.txt
$gene_read_counter_exe  read_bed_intersection_gff.txt > gene_read_counts.txt

