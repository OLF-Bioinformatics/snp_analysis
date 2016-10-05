#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;


#I/O files
my $positionIn = $ARGV[0];
my $vcfIn = $ARGV[1];
my $vcfOut = $ARGV[2];

#Check if 3 arguments
if (scalar(@ARGV) != 3)
{
    print "Usage: perl regionRemover.pl <FilterToAll.txt> <input.vcf> <filtered.vcf>\n";
    exit;
}


#Open the input position file, in read-only mode
open(my $positionInFH, "<", $positionIn) or die "Error opening input vcf file $positionIn: $!\n";

#Store positions in array
my @positions;

#Read position file
while (my $line = <$positionInFH>)
{
    chomp($line); #strip trailing carriage return
    next if $line eq ''; #skip if line is empty.
    push(@positions, $line);
}

#Close position file
close($positionInFH);


#Open input VCF file
open(my $vcfInFH, "<", $vcfIn) or die "Error write to output vcf file $vcfIn: $!\n";

#To store header info
my @header;

#To store VCF fields info
my %vcf;

while (my $line = <$vcfInFH>)
{
    chomp($line); #remove carriage return
    next if $line eq ''; #skip if empty
    
    
    if ($line =~ /^#/) #header line
    {
        push (@header, $line);
        next;
    }
    else #data
    {
        #put VCF line into array
        my @fields = split(/\t/, $line);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SAMPLE) = @fields[0..9];
        my $genotype = (split(/:/, $SAMPLE))[0];

        #Put data in hash of hashes of array
        push(@{ $vcf{$CHROM}{$POS} }, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SAMPLE) unless $genotype eq './.';
    }
    
}

#Close ouput VCF file 
close($vcfInFH);


#Open output vcf file
open(my $vcfOutFH, ">", $vcfOut) or die "Error write to output vcf file $vcfOut: $!\n";

print ($vcfOutFH join("\n", @header), "\n");

foreach my $chromosome (sort keys %vcf)
{
    foreach my $position (sort keys %{ $vcf{$chromosome} })
    {
        print ($vcfOutFH "$chromosome\t$position\t", join("\t", @{ $vcf{$chromosome}{$position} }), "\n");
    }
}


#Close ouput VCF file 
close($vcfOutFH);
