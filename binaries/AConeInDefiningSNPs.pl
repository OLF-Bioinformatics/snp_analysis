#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Array::Utils qw(:all);

#I/O
my $definingSNPs = $ARGV[0];
my $vcfFolder = $ARGV[1];
my $report = $ARGV[2];


#Check if 3 arguments
if (scalar(@ARGV) != 3)
{
    print "Usage: perl AConeInDefiningSNPs.pl <definingSnps.tsv> <vcfFolder> <section2.txt>\n";
    exit;
}


=pod
Designation Chromosome  Position    ALT Label
1   AF2122_NC002945 2134486 T   Group
2   AF2122_NC002945 389472  T   Group
=cut


#####################
#                   #
#   Defining SNPs   #
#                   #
#####################


#Open the input position file, in read-only mode
open(my $positionInFH, "<", $definingSNPs) or die "Error opening input vcf file $definingSNPs: $!\n";

#Store positions in array
my %defining;

#skip header
my $header = <$positionInFH>; #read one line

#Read position file
while (my $line = <$positionInFH>)
{
    chomp($line); #strip trailing carriage return
    next if $line eq ''; #skip if line is empty.
    my ($chrom, $pos) = (split(/\t/, $line))[1..2];
    push( @{ $defining{$chrom} }, $pos) ;
}

#Close position file
close($positionInFH);


#################
#               #
#   VCF files   #
#               #
#################


#open directory with VCF files to look into
opendir(my $vcfDirFH, $vcfFolder) or die "Could not open $vcfFolder: $!\n";

#Create a list of all the VCF files (array)
my @vcfFiles = grep { -f } glob("$vcfFolder/*.vcf");

#Close directory handle
closedir($vcfDirFH);

#create an array for the file handles
my @vcfFilesFH;
foreach my $file (@vcfFiles)
{
    open(my $fh, "<", $file) or die "Can't open < $file: $!\n";
    push(@vcfFilesFH, $fh);
}

#process every VCF file in a row and store content in array
my  %ac1InDefining;

foreach my $handle (@vcfFilesFH)
{
    my $sampleName;
    
    while (my $line = <$handle>)
    {
        chomp($line); #remove carriage return
        next if $line eq ''; #skip if empty
        next if ($line =~ /^##/); #skip if VCF header line
        
        if ($line =~ /^#CHROM/)
        {
            $sampleName = (split(/\t/, $line))[9];
            next;
        }
        
        #put VCF line into array
        my @fields = split(/\t/, $line);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SAMPLE) = @fields[0..9];
        my $AC = (split(/;/, $INFO))[0];
        
        
        #AC=1 postions which are Defining SNPs
        if ($AC eq 'AC=1')
        {
            if (grep { $POS eq $_ } @{ $defining{$CHROM} })
            {
                #Put data in hash of hashes of array
                push(@{ $ac1InDefining{$sampleName}{$CHROM} }, $POS);
            }
        }
    }
}


#close VCF files handles
foreach my $handle (@vcfFilesFH)
{
    close($handle);
}


##############
#            #
#   Report   #
#            #
##############

#Create output file
open (my $reportFH, '>', $report)  or die "Error writing to output file $report: $!\n";

#loop though the hash of hashes of array
foreach my $sample (sort keys %ac1InDefining)
{
    foreach my $chrom (sort keys %{ $ac1InDefining{$sample} })
    {
        if (@{ $ac1InDefining{$sample}{$chrom} }) #if any
        {
            print ($reportFH "AC=1 found in defining SNPs for $sample in chromosome $chrom at position(s): join(", ", @{ $ac1InDefining{$sample}{$chrom} }).\n");
        }
    }
}

close ($reportFH);
