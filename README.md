# UMI (Unique Molecular Identifier) Processing Pipeline

[NCGR] and [WPI] are supported by the [Gordon and Betty Moore Foundation][moore] grant number 4823 to develop single cell RNA seq in plants. This research includes preparation and sequencing of Illumina RNA-seq libraries having Unique Molecular Identifiers (UMI) to enable counting of mRNA molecules. We developed a series of programs in Perl and C to process sequence data containing UMI. Programmers and analysts familiar with handling sequence data with these languages will be able to use and adapt these programs and we hope you will find them useful. However, this is not finished software. Its usefulness will be limited unless you can read and edit C and Perl, and unless you are familiar with sequence files, formats, and common practices. You will need to install libraries, edit makefiles to match your environment, and set paths to executables. Our use case only involves single ended sequencing at this stage. Steps in the workflow and the scripts to use are as follows:

#### Quality filter your sequence data

Use the C project src/c/fastq\_qual\_filter. This requires the [Better String Library](https://github.com/msteinert/bstring) shared object to be in your LD\_LIBRARY\_PATH and bstring.h should be somewhere sensible. This library adds some desireable features to string manipulations in C. Processing large FASTQ files in Perl can be slow which is why C is used here. Go into src/c/fastq\_qual\_filter and edit the makefile. The compiled program takes a FASTQ file and filters on any substring of your choice. This is important because you may wish to apply a stricter quality threshold to the UMI itself than the rest of your sequence, or your UMI position in the read may not be the same as ours. You can filter the entire read on one quality then take the output and then filter just on the UMI. This assumes your FASTQ file has been demultiplexed already, i.e. contains only one index, or that you don't care what the indices are.

``fastq_qual_filter infile goodfile badfile logfile min_ok offset length``

* infile: path to a fastq file
* goodfile: path to output fastq file for passing reads
* badfile: path to output fastq file for failing reads
* logfile: path to file for logging (not used)
* min_ok: pass/fail mean quality threshold
* offset: 0-based position to begin the substring to examine
* length: length of substring to examine
   
   Like:
   ``fastq_qual_filter input.fq good.fq bad.fq log 30 0 10``

   Warning: not much arg checking is done so unexpected behavior may be encountered. For example if you enter a substring longer than the read length, or negative values. Output files will be overwritten without warning.

   #### Remove the UMI from each sequence read
   Use the C project src/c/fastq\_umi\_clipper. This requires the Better String Library (https://github.com/msteinert/bstring) shared object to be in your LD\_LIBRARY\_PATH and bstring.h should be somewhere sensible. Edit the makefile to match your environment. You tell the program how long the UMI is at the start of each read. It replaces the index at the end of the FASTQ header with the UMI. Optionally, it erases more bases at the 5' end of each read. This is used when you make your RNA-seq library using the [Clontech SMRTER](http://www.clontech.com/US/Products/cDNA_Synthesis_and_Library_Construction/cDNA_Synthesis_Kits/SMARTer_Kits#) template switching system, which leaves GGG or GGGG at the 5' end of the read and you don't want these bases  in your alignments. 

  ``fastq_umi_clipper fastqfilename clip_length g_length``

   * fastqfilename: path to a fastq file
   * clip\_length: number of nucleotides to move from the start of the read and put in header 1
   * g\_length: number of nucleotides (Gs) to delete after the UMI (can be 0)

   Warning: not much arg checking is done so unexpected behavior may be encountered. For example if you enter a clip\_length longer than the read length, or negative values. Output files will be overwritten without warning. 

   #### Align your reads to a genomic reference sequence
   Do this however you want. We use [STAR](https://github.com/alexdobin/STAR) specifying no introns. To get to the next step you need a sorted SAM file.

   #### Eliminate duplicate UMI at the the same genomic position
   Run the Perl script umi\_counterizer.pl on your sorted SAM file. This requires the packages [Text::Levenshtein::XS](http://search.cpan.org/~ugexe/Text-Levenshtein-XS-0.503/lib/Text/Levenshtein/XS.pm) and [Try::Tiny](http://search.cpan.org/~ether/Try-Tiny-0.28/lib/Try/Tiny.pm). 

   ``umi_counterizer.pl sorted_samfile max_edit_distance``
   At the same mapping position on the genome, UMI greater than max\_edit\_distance from each other will be considered different UMI and not collapsed into one.

   What it does:
   * Open a sorted SAM file. 
   * Identify UMI duplicates at the same position. 
   * Collapse duplicates pairwise into most abundant UMI
   * Count the UMIs at each position
   * Print to STDOUT

   You end up with the UMI associated with reads mapping to positions on the genome: position, number of UMI, the UMIs (space delimited):
   ```
   Chr01:219194 1 CAAGGTATAG
   Chr01:219196 6 TGATGGAGGG TCAACCGGGT TGGGGAAAGG AAGTAGAAGG TATGTAACAG GCTTCCAGGG
   Chr01:219204 1 GTGGGGTCCT
   Chr01:219221 1 TCCTTTGGTT
   Chr01:219226 2 GGTCAAACCA TGGTCATCCC
   Chr01:219243 2 TGGTTTGTGA CAGTTTGTGC
   Chr01:219245 3 GGGGCCCACT TCTGTATGTA CGCGGGTCCC
   etc.
   ```
Some logging is done and sent to all\_sam\_umis.txt (every UMI observed in the SAM file) and collapsed\_sam\_umis.txt (non-redundant list of UMI found in the SAM file). Both are sorted.
   
   
   #### Find the genomic positions of every UMI
   Next step is to find all of the positions where a UMI maps in the genome. You do with with the C project in src/c/umi\_analyzer\_uthash. Keeping track of positions associated with a UMI requires a hash or map structure, which C lacks. Unfortunately using a scripting language for this task is far too slow. Fortunately there are C implementations of hashes. One of these is [uthash](http://troydhanson.github.io/uthash/userguide.html), which is implemented entirely as preprocessor macros. So put uthash.h somewhere sensible and there is no need to link to a library. The input to this program is the output of umi\_counterizer.pl. There is no flexibility here, and not much error checking. The emphasis here is solely on speed, making important assumptions about the structure of the input. This makes it possible to use fgets() instead of getline() to peel lines off a filehandle, and strtok\_r() expects space delimited input lines. You need to edit exe.c before compiling. Specify the length of the UMI and the maximum length of the position string in bytes. Also specify the maximum expected length of a line in the umi counts file. 

Output is each UMI and where it maps in the genome. This file is usually called umi\_positions.txt.
   ```
   TGTAAATTGG Chr01:76001 Chr01:76002 
   GTTGCGGGCT Chr01:76002 
   ACAGGTTAGC Chr01:76007 
   GGGCAGTGTT Chr01:76059 
   GTTTGGGGGG Chr01:76065 Chr01:6054609 
   TAAATGAGGG Chr01:76075 
   GCGTTGTCTT Chr01:76079 
   CGCGTTGTCT Chr01:76080 
   TCTGTGGCGC Chr01:76087 
   GCAGTGAGGG Chr01:76103 Chr01:76104 
   GGTTTTGTGG Chr01:76122 Chr01:76123 
   TGGGTTGGGG Chr01:76122 Chr01:76123 Chr01:76124 Chr01:10482858 
   etc.
   ```
   #### Count UMI at each genomic position
   Next step is to count unique UMI at any position where a UMI is found. This is one more step towards what we really want, which is a count of UMI for every gene, to be used in differential gene expression analysis. Run the Perl script umi\_clusterizer.pl using the umi\_positions.txt file from the previous step as input. The output will be a BED format file (usually position\_umi\_counts.bed) consisting of positions on the genome and the number of unique UMI found in reads mapping to those positions. BED format is used to enable merging with GFF in the next step:

   ```
Chr01	10070102	10070102	1
Chr01	10070498	10070498	6
Chr01	10070499	10070499	8
Chr01	10070500	10070500	4
Chr01	10071648	10071648	1
Chr01	10076521	10076521	1
Chr01	10076540	10076540	2
Chr01	10076541	10076541	4
Chr01	10076542	10076542	1
Chr01	10076543	10076543	3
   etc.
   ```

   #### Find the genes in a GFF file matching the UMI positions
   This is done by taking the positions of the UMI in the previous step and finding the corresponding annotations in a GFF file. Run bedtools thus:

   ``bedtools intersect -a position_umi_counts.bed -b file.gff -wa -wb | pcregrep '\tgene\t' > umi_bed_intersection_gff.txt``

   * gff\_file: GFF fprmatted annotations for your organism. Positions are grouped on the feature "gene".
   * position_umi_counts.bed: output file from previous step

The output is a file with lines that connects the UMI coordinates to their cognate genes (umi\_bed\_intersection_gff.txt):

```
Chr01	10070102	10070102	1	Chr01	phytozomev10	gene	10067999	1007
0452	.	+	.	ID=Pp3c1_13800;gene_id=Pp3c1_13800;ancestorIdentifier=Pp3c1_13800.v3
.1
Chr01	10076521	10076521	1	Chr01	phytozomev10	gene	10076060	1007
8901	.	+	.	ID=Pp3c1_13830;gene_id=Pp3c1_13830;ancestorIdentifier=Pp3c1_13830.v3
.1
Chr01	10076540	10076540	2	Chr01	phytozomev10	gene	10076060	1007
8901	.	+	.	ID=Pp3c1_13830;gene_id=Pp3c1_13830;ancestorIdentifier=Pp3c1_13830.v3
.1
Chr01	10076541	10076541	4	Chr01	phytozomev10	gene	10076060	1007
8901	.	+	.	ID=Pp3c1_13830;gene_id=Pp3c1_13830;ancestorIdentifier=Pp3c1_13830.v3
.1
Chr01	10076542	10076542	1	Chr01	phytozomev10	gene	10076060	1007
8901	.	+	.	ID=Pp3c1_13830;gene_id=Pp3c1_13830;ancestorIdentifier=Pp3c1_13830.v3
.1
Chr01	10076543	10076543	3	Chr01	phytozomev10	gene	10076060	1007
8901	.	+	.	ID=Pp3c1_13830;gene_id=Pp3c1_13830;ancestorIdentifier=Pp3c1_13830.v3
.1

```
   #### Extract the UMI counts per gene
Execute the perl script gene\_umi\_counter.pl on the above output:

``gene_umi_counter.pl umi_bed_intersection_gff.txt > gene_umi_counts.txt``

   The output is a gene and UMI count table like this, which is the start of a differential gene expression analysis experiment:
   
   ```
Pp3c5_13519     11047
Pp3c18_8100	    9694
Pp3c22_5610	    9088
Pp3c21_3950	    7194
Pp3c3_18750	    6002
Pp3c19_21160	5662
Pp3c7_14450 	4566
Pp3c13_7900	    4109
Pp3c19_20900	4088
Pp3c13_5930	    3325
Pp3c21_9980 	3172
Pp3c3_10740	    3055
Pp3c3_10750	    3054
Pp3c9_5310	    2812
Pp3c9_3440	    2588
Pp3c6_12510	    2520
   etc.
   ```

 #### Extract the read counts per gene
 
 To make a comparison between UMI counts and read counts per gene you will need an equivalent table mapping genes to their read counts.  The first step is to run sam_\to\_bed.pl in order to arrive at a BED file that indicates the number of mapping to each coordinate in the genome:
 
``sam_to_bed.pl Aligned.sorted.out.uniq.sam > position_read_counts.bed``

Then merge this file with your GFF file to identify genes in the GFF file corresponding to the mapped reads:

``bedtools intersect -a position_read_counts.bed -b $gff_file -wa -wb | pcregrep '\tgene\t' > read_bed_intersection_gff.txt``

Finally, count up the reads mapping to each gene:

``gene_read_counter.pl  read_bed_intersection_gff.txt > gene_read_counts.txt``

The output is a table mapping genes to reads. This can be all be executed with the bash wrapper main_wrapper.sh. It is advisable to use commenting to execute the steps one at a time until your are confident they will run without errors. This set of scripts automates the basic pipeline, but will require customization and additional code depending on a variety of factors such as the header format of your fastq files, where in the reads UMI are to be found and the protocol used for RNA-seq library preparation. The gene-UMI and gene-read tables are starting points to demonstrate the effect of correcting PCR amplification bias with UMI, and offer a starting point for gene expression analysis. Depending on your needs you may want to filter the GFF file differently than simply on the feature gene.
 
 #### Dependencies
 
C programs require the [Better String Library](https://github.com/msteinert/bstring) and [uthash.h](http://troydhanson.github.io/uthash/userguide.html). Perl modules are noted in use statements in each program. Bioinformatics tools are managed in a [Conda environment](bioconda.github.io).
 



---




   [ncgr]: <http://www.ncgr.org>
   [moore]: <http://www.moore.org>
   [wpi]: <https://www.wpi.edu>
