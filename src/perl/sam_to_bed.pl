#!/home/linux/perl/bin/perl

use strict;

# Feed on a SAM file. The goal here is to make a BED file that specifies the number of reads 
# at each read start position, like this:
#
# NC_000067.6     106063295       106063295       2
# NC_000067.6     106936257       106936257       1
# NC_000067.6     107102708       107102708       1
# NC_000067.6     107202703       107202703       2
# 
# It is BED format but has nothing to do with showing features in a browser. It is only used 
# for compatibility with BEDtools so that we can intersect these read counts with a GFF file
# and count the number of reads in each gene.
#
# The reason NOT to simply intersect the BAM and GFF files is that the UMI counts per gene 
# calculator only considers the starting position of the read mapping, not the whole read.
# We want to compare UMI and read counts per gene so need to count them the same way.
#


unless ($ARGV[0])
{
   print "\nUsage: $0 uniquely_aligning_samfile\n\n";
   exit 1;
}

my $file = shift;

open (SAM, "<$file")|| die "Couldn't open the file $!";

my $h; 

while (<SAM>)
{
   next if(/^@/);
   
   my ($qname,$flag,$rname,$pos,$mapq,$cigar,$rnext,$pnext,$tlan,$seq,$qual,$tag1) = split /\s+/;
      
   $h->{$rname}->{$pos}++; 
  
}


for my $contig (sort keys %$h)
{
   for my $pos (sort keys %{$h->{$contig}})
   {
      print "$contig\t$pos\t$pos\t$h->{$contig}->{$pos}\n";
   }
}
