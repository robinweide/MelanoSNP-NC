#!/bin/env perl
 
use strict;
use warnings;
 
my $usage = "Usage: $0 <infile.jaspar>\n";
my $infile = shift or die $usage;
 
open(IN,'<',$infile) || die "Could not open $infile: $!\n";
while(<IN>){
   chomp;
#>MA0004.1 Arnt
#4       19      0       0       0       0
#16      0       20      0       0       0
#0       1       0       20      0       20
#0       0       0       0       20      0
   if (/^>/){
      s/\s+/_/g;
      print "$_\n";
      my $freq = [];
      for(1..4){
         chomp(my $line = <IN>);
         my @line = split(/\s+/, $line);
         push($freq,\@line);
      }
      for(my $i=0; $i < scalar(@{$freq->[0]}); ++$i){
         my $line_to_print = '';
         for(my $j=0; $j < scalar(@{$freq}); ++$j){
            $line_to_print .= "$freq->[$j][$i]\t";
         }
         $line_to_print =~ s/\t$//;
         print "$line_to_print\n";
      }
      print "<\n";
   }
}
close(IN);
 
exit(0);