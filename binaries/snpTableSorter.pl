#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use List::MoreUtils qw(firstidx);
use Data::Dumper qw(Dumper);


#I/O
my $tableIN = $ARGV[0]; #.table.txt
my $sampleOrderIN = $ARGV[1]; #cleanedAlignment.txt
my $organizedOUT = $ARGV[2]; #.organized_table.txt


#Check if 3 arguments
if (scalar(@ARGV) != 3)
{
   print "Usage: perl snpTableSorter.pl <table.txt> <cleanedAlignment.txt>  <organized_table.txt>\n";
   exit;
}


########################
#                      #
#   Parse table file   #
#                      #
########################


#open table.txt
open(my $tableFH, '<', $tableIN) or die "Can't read $tableIN: $!\n";

my %table;

chomp(my $header = <$tableFH>); #first line -> reference_pos
my @refPosList = split(/\t/, $header);
#@refPos = $refPos[1..$#refPos]; #"delete" fist elment of array

while (my $line = <$tableFH>)
{
    chomp($line); #remove carriage return
    next if $line eq ''; #skip if empty
    
    my @fields = split(/\t/, $line);
    my $sample = $fields[0];
    for my $i (1..$#fields)
    {
        push( @{ $table{$sample}{$refPosList[$i]} }, $fields[$i] );
    }
}

#print Dumper \%table;


###################################
#                                 #
#   Order positions and samples   #
#                                 #
###################################


my %countPos;
#my %countSam; #Commented because order comes from the RAxML output

#number of REF alleles
foreach my $sample(sort keys %table)
{
    my $cntREF = 0;
    foreach my $pos (sort keys %{ $table{$sample} })
    {
        my $REF = @{ $table{'reference_call'}{$pos} }[0];
        my $allele = @{ $table{$sample}{$pos} }[0];
        
        $countPos{$pos} += 1 if ( $allele eq $REF );
#        $countSam{$sample} += 1 if ( $allele eq $REF); #Commented because order comes from the RAxML output
    }
}

#print Dumper \%counts;

#Order positions based on REF allele counts and save to array
my @orderedPos;
foreach my $pos (sort { $countPos{$a} <=> $countPos{$b} } keys %countPos)
{
#	print $counts{$pos}; #count values
#	print "\n";
	push(@orderedPos, $pos); #position names
} 

#Commented because order comes from the RAxML output
#Order samples based on REF allele counts and save to array
#my @orderedSam;
#foreach my $pos (sort { $countSam{$b} <=> $countSam{$a} } keys %countSam)
#{
##	print $counts{$pos}; #count values
##	print "\n";
#	push(@orderedSam, $pos); #position names
#}



#open cleanedAlignment.txt
open(my $samFH, '<', $sampleOrderIN) or die "Can't read $sampleOrderIN: $!\n";

#skip fisrt line
chomp(my $header1 = <$samFH>); #first line -> reference_pos

my @orderedSam;

#read the rest of the file
while (my $line = <$samFH>)
{
    chomp($line); #remove carriage return
    next if $line eq ''; #skip if empty
    
    push(@orderedSam, $line);
}

#Commented because order comes from the RAxML output
#"reference_call" has to be the first
#my $index = firstidx { $_ eq "reference_call" } @orderedSam; #find it's index
#splice(@orderedSam, $index, 1); #Remove "reference_call" from array
#unshift(@orderedSam, "reference_call"); #Add back "reference_call" at the beginning of the array

#refine position sorting by determining the row number of the first ALT allele
my %rank;
my $cnt = 0;

foreach my $sample (@orderedSam)
{
	$cnt += 1;
	foreach my $pos (@orderedPos)
	{
		my $REF = @{ $table{'reference_call'}{$pos} }[0];
		my $allele = @{ $table{$sample}{$pos} }[0];
		
		if ( $allele ne $REF && !exists $rank{$pos} && !defined $rank{$pos} )
		{
			$rank{$pos} = $cnt;
		}
	}
}


#Order positions based on REF allele counts and save to array
my @orderedPosTuned;
#foreach my $pos (sort { $countPos{$a} <=> $countPos{$b} || $rank{$a} <=> $rank{$b} } keys %countPos)
foreach my $pos (sort { $rank{$a} <=> $rank{$b} || $countPos{$a} <=> $countPos{$b} } keys %countPos)
{
#	print $counts{$pos}; #count values -> debug
#	print "\n";
	push(@orderedPosTuned, $pos); #position names
} 


###########################
#                         #
#   Output sorted table   #
#                         #
###########################


#Write to file
open(my $tmpFH, '>', $organizedOUT) or die "Could not write to $organizedOUT: $!\n";

print { $tmpFH } ("reference_pos\t");

my $positions = join("\t", @orderedPosTuned);
print { $tmpFH } ("$positions\n");


foreach my $sample (@orderedSam)
{
	print { $tmpFH } ("$sample");
	
	foreach my $pos (@orderedPosTuned)
	{
		my $allele = @{ $table{$sample}{$pos} }[0];
		print { $tmpFH } ("\t$allele");
	}
	print { $tmpFH } ("\n");
}
