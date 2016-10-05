#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Spreadsheet::XLSX;

#I/O files
my $excelIn = $ARGV[0];
my $tsvOut = $ARGV[1];

#Check if 2 arguments
if (scalar(@ARGV) != 2)
{
    print "Usage: perl parseXLSX.pl <input.xlsx> <output.tsv>\n";
    exit;
}

#Open excel file
my $excel = Spreadsheet::XLSX -> new($excelIn);

#Open desired worksheet
my $sheet = $excel -> {Worksheet}[0]; #['New groupings']};

#Open output file for writing
open(my $tsvOutFH, '>', $tsvOut) or die "Could not write to $tsvOut: $!\n";

#Set row and column ranges to convert
#$sheet -> {MaxRow} ||= $sheet -> {MinRow};
my $row_min = $sheet -> {MinRow};
my $row_max = $sheet -> {MaxRow};



foreach my $row ($row_min..$row_max)
{
#   $sheet -> {MaxCol} ||= $sheet -> {MinCol};
    my $col_min = $sheet -> {MinCol};
    my $col_max = $sheet -> {MaxCol};
    
    #store data in array
    my @data;
    foreach my $col ($col_min..$col_max)
    {
        my $cell = $sheet -> {Cells} [$row] [$col];
        push @data, defined($cell -> {Val}) ? $cell -> {Val} : ''; #accounts for empty cells
    }
    
    #expand the ranges, if any, in the array
    
    
    #print array content to output file
    print $tsvOutFH join("\t", @data), "\n";
}

# Close output file
close($tsvOutFH);
