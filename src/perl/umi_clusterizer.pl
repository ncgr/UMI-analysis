#!/home/linux/perl/bin/perl

use strict;

# Feed on a list of UMIs and their positions
# umi_positions.txt
# TGCAGCGTTG Chr01:629718 Chr06:10727258 Chr15:12754424 Chr16:7912214
#
# These positions must be sorted in ascending order
#
# Group identical UMIs into clusters if they are within N nucleotides of each other. A UMI may have 
# multiple positions that occur on the same chromosome or transcript that are legitimate, i.e. occur from when
# the same transctipt chooses the same UMI. However, when positions are clustered close to one another, something
# fishy is happening. How close is close enough to be suspicious? That's the second argument. 1 or 2 or 3.
#
# Output a BED file:
#
# Chr01   10069826        10069826        1
# Chr01   10070497        10070497        2
# Chr01   10070498        10070498        3
# Chr01   10070499        10070499        1
# Chr01   10076520        10076520        1
# 
# The laast column is the number of UMI mapping to that position


unless ($ARGV[1])
{
   print "\nUsage: $0 umi_positions_file max_spacing\n\n";
   exit 1;
}

# umi positions file and maximum spacing between positions that we want to consider interesting.
my ($umifile,$max_spacing) = @ARGV;

open (UMI, "<$umifile")|| die "Couldn't open the file $!";

my $umi_to_position; #hash reference. ref->umi->{cluster_coords}

# open each line at a time
while (<UMI>)
{
   my ($umi,@positions) = split /\s+/;  
   my $umi_data; #hash reference for each UMI
   
   for my $position (@positions)
   {
      my ($chr,$nuc) = split /:/,$position;   
      push @{$umi_data->{$chr}},$nuc;   # $umi_data->{$chr} is an array reference
   }
   
   # Now the data are packed by chromosome. Each chromosome or transcript or whatever unit
   # has an arrayref of positions.
   
   # Go through positions on the current chromosome for the current UMI  
   CHR: while ( my ($chr,$nuc_positions_arrayref) = each %$umi_data)
   {
      # This will be the coord on the chromosome that will be recorded as a UMI cluster.
      # Set it to the smallest or the only one in the array.
      my $cluster_position = $chr . ":" . $nuc_positions_arrayref->[0];
      
      ###print "$umi   $chr   @$nuc_positions_arrayref\n";
      
      # If there is only one position on the chromosome there are no actual spacings.Spacing is an array of one. 
      if (scalar @$nuc_positions_arrayref == 1)
      {
         push@{$umi_to_position->{$umi}},$cluster_position;         
         ###print "   $cluster_position\n\n";         
         next CHR;
      }
      
      # Else there is more than one position on the chromosome. Make the clusters.
      COORD: for ( my $i = 0 ; $i <= $#$nuc_positions_arrayref-1 ; $i++) 
      { 
         my $spacing = $nuc_positions_arrayref->[$i+1] - $nuc_positions_arrayref->[$i];
         
         # If the spacing is <= threshold, cluster position stays the same as 
         # already set
         if ($spacing <= $max_spacing)
         {
            # If $i+1 is the last one, then record the cluster and go to the next chromosome.
            if($i+1 == $#$nuc_positions_arrayref)
            {
               push @{$umi_to_position->{$umi}},$cluster_position;         
               ###print "   $cluster_position\n\n";  
               next CHR;                     
            }
         }
                 
         # If the spacing is greater than the threshold, then record the current cluster
         # position for the UMI 
         if($spacing > $max_spacing)
         {
            push @{$umi_to_position->{$umi}},$cluster_position;         
            ###print "   $cluster_position\n\n";
            
            # If $i+1 is the last one, then the last one counts as a cluster.
            # Record it as a cluster position then stop here and go to the next chromosome
            if($i+1 == $#$nuc_positions_arrayref)
            {
               $cluster_position = $chr . ":" . $nuc_positions_arrayref->[$i+1];
               push @{$umi_to_position->{$umi}},$cluster_position;         
               ###print "   $cluster_position\n\n";  
               next CHR;                     
            }
            # If $i+1 is not the last one, then set the current cluster position to $i+1
            # and go to the next iteration. This is another cluster for this UMI on the 
            # same chromosome
            else
            {
               $cluster_position = $chr . ":" . $nuc_positions_arrayref->[$i+1];
               next COORD;
            }
            
         }
         
      }
      
   }   
   
}



# Turn the umi_to_position hash into position and UMI counts

my %position_umi_counts;

while ( my($umi,$position_array_ref) = each %$umi_to_position )
{
   for my $position(@$position_array_ref)
   {
      $position_umi_counts{$position}++;
   }  
}


# for my $position (sort keys %position_umi_counts)
# {
#   print "$position  $position_umi_counts{$position}\n";
# }


for my $position (sort keys %position_umi_counts)
{
  my ($chr,$coord) = split /:/,$position;
  print "$chr\t$coord\t$coord\t$position_umi_counts{$position}\n";
}



# output:
# Pp3c10_10290V3.1:827  1
# Pp3c10_10290V3.1:846  1
# Pp3c10_10290V3.1:87  2
# Pp3c10_10290V3.1:88  4
# Pp3c10_10290V3.1:89  19
# Pp3c10_10290V3.1:90  8
# Pp3c10_10290V3.1:91  1
# Pp3c10_10290V3.1:939  1
# Pp3c10_10290V3.1:969  1
# Pp3c10_10490V3.3:2223  1
# Pp3c10_10710V3.2:2  1
# Pp3c10_10710V3.2:4  1
# Pp3c10_10710V3.2:8  1
# Pp3c10_11000V3.3:435  1
#etc.
