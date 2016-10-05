#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Number::Range;

#I/O files
my $tsvIn = $ARGV[0];
my $outDir = $ARGV[1];

#Check if 2 arguments
if (scalar(@ARGV) != 2)
{
    print "Usage: perl rangeExpander.pl <filterFile.txt> <output_folder>\n";
    exit;
}

#For debug
#my $tsvIn = "/home/bioinfo/analyses/mbovis_script2/2016-09-20at10h41m32s-FilterFiles/filterFile.txt";
#my $outDir = "/home/bioinfo/analyses/mbovis_script2/2016-09-20at10h41m32s-FilterFiles";

#Open tsv file
open(my $inFH, '<', $tsvIn) or die "Could not open input file $tsvIn: $!\n";

#Read fisrt line (header)
chomp(my $fistLine = <$inFH>);   #First line is read here

#create array with header
my @header = split(/\t/, $fistLine);
my $headerLastIndex = $#header;

#hash to store column info
my %groups;

while (my $line = <$inFH>)
{
    chomp($line); #remove carriage return
    next if $line eq ''; #skip if empty
    
    my @items = split(/\t/, $line);
    
    for (my $i = 0; $i <= $headerLastIndex; $i++)
    {
        if($items[$i])
        {
            push @{ $groups{$header[$i]} }, $items[$i];
        }       
    }
    
}

#expand the ranges
foreach my $key (sort keys %groups)
{
    for my $i (0 .. $#{$groups{$key}})
    {
        my $value = $groups{$key}[$i];
        if ($value =~ m/-/)
        {
            $value =~ s/-/../;
            my @newRange = Number::Range -> new($value) -> range;
            $groups{$key}[$i] = \@newRange;
        }
    }
}

#output the hash, one file per key

foreach my $key (sort keys %groups)
{
    my $outFile = $outDir . "/" . $key . ".txt";
    
    #Create output file (handle)
    open(my $outFH, '>', $outFile) or die "Could write to output file $outFile: $!\n";
    
    for my $i (0 .. $#{$groups{$key}})
    {
        if(ref($groups{$key}[$i]) eq 'ARRAY') #if an array (previously a range)
        {
            foreach my $item (@{ $groups{$key}[$i] })
            {
                print $outFH $item, "\n";
            }       
        }
        else
        {
            print $outFH $groups{$key}[$i], "\n";
        }
    }
    
    #Close output file handle
    close($outFH);
}
