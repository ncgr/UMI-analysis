#!/home/linux/perl/bin/perl

use strict;

#
# Feed on merged BED and GFF files. The file is read_bed_intersection_gff.txt made
# by Bedtools joining a BED file: position_read_counts.bed (this BED file was made by 
# sam_to_bed.pl)...
# 
# NC_000067.6     100448946       100448946       1
# NC_000067.6     100616625       100616625       1
# NC_000067.6     100635787       100635787       1
# NC_000067.6     10064082        10064082        1
#
# with a GFF file. The intersection file has a line for each overlapping segment.
# Since the BED file has segments of length 1, this means every line in the intersection
# file describes one nucleotide in a gene, having a certain number of reads mapping there.
# These reads counts are added up for each gene.
#


unless ($ARGV[0])
{
   print "\nUsage: $0 read_bed_intersection_gff\n\n";
   exit 1;
}

my $file = shift;

open (UMI, "<$file")|| die "Couldn't open the file $!";

my $h; 

while (<UMI>)
{
   my ($chr1,$start,$stop,$read_count,$chr2,$method,$feature,$start2,$stop2,$score,$strand,$phase,$attributes) = split /\s+/;  
   my ($gene) = $attributes =~ /ID=(.*?);/;
   
   if($gene ne "")
   {
      $h->{$gene} += $read_count;
   }  
  
}


for my $gene (sort {$h->{$b} <=> $h->{$a}} keys %$h)
{
   print "$gene\t$h->{$gene}\n";
}
