#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;


#I/O files
my $positionIn = $ARGV[0];
my $chromFile = $ARGV[1];
my $vcfIn = $ARGV[2];
my $vcfOut = $ARGV[3];

#Check if 3 arguments
if (scalar(@ARGV) != 4)
{
    print "Usage: perl regionRemover.pl <FilterToAll.txt> <chroms.txt> <input.vcf> <filtered.vcf>\n";
    exit;
}


###################
#                 #
#   Chromosomes   #
#                 #
###################


#Open the input position file, in read-only mode
open(my $chromFileFH, "<", $chromFile) or die "Error opening input vcf file $chromFile: $!\n";

#store chromosomes in array
my @chroms;

#Read position file
while (my $line = <$chromFileFH>)
{
    chomp($line); #strip trailing carriage return
    next if $line eq ''; #skip if line is empty.
    push(@chroms, $line);
}

#Close position file
close($chromFileFH);


#################
#               #
#   Positions   #
#               #
#################


#Open the input position file, in read-only mode
open(my $positionInFH, "<", $positionIn) or die "Error opening input vcf file $positionIn: $!\n";

#Store positions in array
my %positions;

#Read position file
while (my $line = <$positionInFH>)
{
    chomp($line); #strip trailing carriage return
    next if $line eq ''; #skip if line is empty.
    
    #If only one chromosome (Mycobaterium)
    #filter files has only one column
    if ( scalar(@chroms) == 1)
    {
        my $chrom = $chroms[0];
        my $pos = $line;
        push (@{ $positions{$chrom}{$pos} }, "1");
    }
    else # if (scalar(@chroms) > 1)
    {
        my ($chrom, $pos) = split(/\t/, $line);
        push (@{ $positions{$chrom}{$pos} }, "1");
    }
    #If no chrom is present, this error is handled in the parent bash script
}

#Close position file
close($positionInFH);


###########
#         #
#   VCF   #
#         #
###########


#Open output vcf file
open(my $vcfOutFH, ">", $vcfOut) or die "Error write to output vcf file $vcfOut: $!\n";

#Open input VCF file
open(my $vcfInFH, "<", $vcfIn) or die "Error writing to output vcf file $vcfIn: $!\n";


while (my $line = <$vcfInFH>)
{
    chomp($line); #remove carriage return
    next if $line eq ''; #skip if empty
    
    if ($line =~ /^#/) #header line
    {
        print ($vcfOutFH $line . "\n");# print it to output file
    }
    else #data
    {
        my ($chrom, $pos) = (split(/\t/, $line))[0..1];
        
        if ( exists( $positions{$chrom}{$pos} ) ) #in the list of positions to filter
        {
            next;
        }
        else
        {
            print ($vcfOutFH $line . "\n");
        }
    }
    
}

#Close input VCF file 
close($vcfInFH);

#Close ouput VCF file 
close($vcfOutFH);
