#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Data::Dumper;
use List::Util qw(sum);
use Math::Round;

#I/O files
my $positionIn = $ARGV[0];
my $vcfFolder = $ARGV[1];
my $qualityOut = $ARGV[2];

#Check if 3 arguments
if (scalar(@ARGV) != 3)
{
    print "Usage: perl qualityExtractor.pl <position.txt> <vcf_folder> <quality.txt>\n";
    exit;
}

#open directory with VCF files to look into
opendir(my $vcfDirFH, $vcfFolder) or die "Could not open $vcfFolder: $!\n";

#Create a list of all the VCF files (array)
my @vcfFiles = grep { -f } glob("$vcfFolder/*.vcf");

#Close directory handle
closedir($vcfDirFH);

#Put all VCF files in a hash hashes of array
#Chromosome -> position -> mapping scores

#create an array for the file handles
my @vcfFilesFH;
foreach my $file (@vcfFiles)
{
    open(my $fh, "<", $file) or die "Can't open < $file: $!\n";
    push(@vcfFilesFH, $fh);
}

my %qual;

foreach my $handle (@vcfFilesFH)
{
    while (my $line = <$handle>)
    {
        chomp($line); #remove carriage return
        next if $line eq ''; #skip if empty
        next if ($line =~ /^#/); #skip if VCF header line
        
        #put VCF line into array
        my @fields = split(/\t/, $line);
        my ($c, $p, $info) = @fields[0..1,7]; #chromosome, position, info field
        
        my $MQ;

        # AF2122_NC002945 306698  .   N   .   .   .   .   GT  ./.
        if ($info eq ".")
        {
            $MQ = 0;
        }
        else
        {
            #extract MQ from the info field
            my @inf = split(/;/, $info);
            
            my @info = grep(/^MQ/, @inf);
            $MQ = $info[0];
            $MQ =~ s/MQ=//;
        }
        
        #http://docstore.mik.ua/orelly/perl2/prog/ch09_02.htm
        push @{ $qual{$c}{$p} }, $MQ;
    }
}

#close VCF files handles
foreach my $handle (@vcfFilesFH)
{
    close($handle);
}

#Create a new hash of hashes of array with the mean mapping score for each chromosome/position
my %averageQual;

#loop through the hash with all the quality values
foreach my $chromKey (sort keys %qual)
{
    foreach my $posKey (sort keys %{ $qual{$chromKey} }) #dereference they key
    {
        my @qualities = @{ $qual{$chromKey}{$posKey} };
        my $meanMQ = sum(@qualities)/scalar(@qualities); #average quality
        $meanMQ = round($meanMQ); #round quality
        $averageQual{$chromKey}{$posKey} = $meanMQ; #save to new hash of hashes
    }
}

#open output vcf file
open(my $qualityOutFH, ">>", $qualityOut) or die "Error write to output vcf file $qualityOut: $!\n";

#write header to output file
my $header = "reference_pos\tmap-quality\n";
print $qualityOutFH "$header";

#Open the input position file, in read-only mode
open(my $positionInFH, "<", $positionIn) or die "Error opening input vcf file $positionIn: $!\n";

#store position in hash
my %positions;

#Read position file
while (my $line = <$positionInFH>)
{
    chomp($line); #strip trailing carriage return
    next if $line eq ''; #skip if line is empty.
    
    my $info = (split(/ /, $line))[1];
    my ($chrom, $pos) = split(/-/, $info); #split fields by space
    $positions{$chrom}{$pos} = $averageQual{$chrom}{$pos}; #no value added to the hash of hashes
}

close($positionInFH);

foreach my $chromKey (sort keys %positions)
{
    foreach my $posKey (sort keys %{ $positions{$chromKey} })
    {
        
        print $qualityOutFH "$chromKey-$posKey\t$positions{$chromKey}{$posKey}\n";
    }
}

close($qualityOutFH);
