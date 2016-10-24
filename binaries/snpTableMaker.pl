#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use Array::Utils qw(:all);
use Bio::SeqIO;


#I/O files
my $genomeFile = $ARGV[0];
my $vcfFolder = $ARGV[1];
my $minQual = $ARGV[2]; #QUAL=150 # Minimum quality for calling a SNP
my $minAltQual = $ARGV[3]; #highEnd=200 # QUAL range to change ALT to N
my $ac1Out = $ARGV[4];
my $countOut = $ARGV[5];
my $fastaOutFolder = $ARGV[6];
my $fastaTable = $ARGV[7];

#Check if 7 arguments
if (scalar(@ARGV) != 8)
{
   print "Usage: perl snpTableMaker.pl <ref.fasta> <vcfFolder> <minQual> <minAltQual> <AC1Report.txt> <section4.txt> <fastaOutFolder> <fastaTable.tsv>\n";
   exit;
}


#################
#               #
#   Functions   #
#               #
#################


sub uniq {
    my %h;
    return grep { !$h{$_}++ } @_
}


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

#AC=1 report
open(my $ac1OutFH, ">", $ac1Out) or die "Error writing to output file $ac1Out: $!\n";

#Need hashes in case multiple chromosomes
my %vcfs;
my %AC1;
my %AC2;

foreach my $handle (@vcfFilesFH)
{
    my $sample;
    
    while (my $line = <$handle>)
    {
        chomp($line); #remove carriage return
        next if $line eq ''; #skip if empty
        next if ($line =~ /^##/); #skip if VCF header line
        
        if ($line =~ /^#CHROM/)
        {
            $sample = (split(/\t/, $line))[9];
            next;
        }
        
        #put VCF line into array
        my @fields = split(/\t/, $line);
        my ($CHROM, $POS, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SAMPLE) = @fields[0..9];
        my $AC = (split(/;/, $INFO))[0];

        #AC = Alternative allele count
        #The variant caller has been used in diploid mode, even though we're working on bacteria
        #AC=2 means that the sample is homozygote for the ALT allele (alternative allele is found) 
        #AC=1 means that the sample has both REF and ALT alleles.
        #Since Brucella spp. and Mycobacterium spp. are bacteria, they should should only have one allele present at the time.
        #Thus, AC=1 suggest that sample might contain multiple isolates.

        #AC=1
        push (@{ $AC1{$sample}{$CHROM} }, $POS) if ($AC eq 'AC=1' && $QUAL > 0);
        
        #AC=2
        push (@{ $AC2{$sample}{$CHROM} }, $POS) if ($AC eq 'AC=2' && $QUAL > $minQual);
                
        #Put data in hash of hashes of array
        push(@{ $vcfs{$sample}{$CHROM}{$POS} }, $ID, $REF, $ALT, $QUAL, $FILTER, $INFO, $FORMAT, $SAMPLE);
    }
}

#close VCF files handles
foreach my $handle (@vcfFilesFH)
{
    close($handle);
}

#Put all AC=2 positions of all VCF files in a single array
my %allAC2;
foreach my $sample ( sort keys %AC2)
{
    foreach my $chrom (sort keys %{ $AC2{$sample} } )
    {
        push (@{ $allAC2{$chrom} }, @{ $AC2{$sample}{$chrom} });
    }
}

#Join, remove duplicates and sort allAC2 hash of arrays
#The has is always needed to keep the chromosome information
foreach my $chrom (sort keys %allAC2)
{
    @{ $allAC2{$chrom} } = sort { $a <=> $b } (uniq (@{ $allAC2{$chrom} }));
}


#Find AC=1 also in AC=2 (intersection)
my %finalAC1;

foreach my $sample ( sort keys %AC1)
{
    foreach my $chrom (sort keys %{ $AC1{$sample} } )
    {
        if ($AC2{$sample}{$chrom})
        {
            my %isect;
            my @test = intersect(@{ $AC1{$sample}{$chrom} }, @{ $allAC2{$chrom} });
            #AC=1 also in AC=2 (intersection)
            push ( @{ $isect{$chrom} },  uniq(sort { $a <=> $b } (intersect(@{ $AC1{$sample}{$chrom} }, @{ $allAC2{$chrom} }))) ); #sort numerically and keep unique
            
            #Store AC=1 found in AC=2 in a hash of hashes of array for each sample (VCF file)
            push ( @{ $finalAC1{$sample}{$chrom} }, @{ $isect{$chrom} });
        }
    }
}

#Report AC=1 also in AC=2 per sample
foreach my $sample (sort keys %finalAC1)
{
    foreach my $chrom (sort keys %{ $finalAC1{$sample} })
    {
        print ($ac1OutFH "$sample\nAC=1 is also found in AC=2 in chromosome $chrom at position(s): " . join(', ', @{ $finalAC1{$sample}{$chrom} }) . "\n\n");
    }
}

#Close AC1 report file
close($ac1OutFH);


##################
#                #
#   Ref genome   #
#                #
##################


#Parse reference genome fasta file into an hash (works if multiple chromosomes)
my %ref;

#use Bio::SeqI
my $seqIO = Bio::SeqIO -> new( -format => "fasta",
                            -file => "$genomeFile");

#process each chromosome
while (my $seq = $seqIO -> next_seq)
{
    # process each seq
    #push value to hash of arrays
    push ( @{ $ref{$seq->id} }, $seq->seq() );
}


################
#              #
#   All SNPs   #
#              #
################


#Filtered SNP positions
#Replace missing SNP with REF
#Only incluse AC=1 positions if found as AC=2 in other samples
my (%fastas, %counts);
my $allele;

#variable to store the last value of the next loop
#my %hash;
my $s;
#my $c;
#my $p;

foreach my $sample ( sort keys %AC2)
{
    $s = $sample;
    
    foreach my $chrom (sort keys %{ $AC2{$sample} } )
    {
        foreach my $pos ( sort { $a <=> $b } (@{ $allAC2{$chrom} }) )
        {
            
            #position was genotyped in sample
            # or is AC=1, but was also found in AC=2
            if( grep(/\b$pos\b/, @{ $AC2{$sample}{$chrom} }) || grep(/\b$pos\b/, @{ $finalAC1{$sample}{$chrom} }) ) #"\b is for word boundary -> exact word match"
            {
                $allele = @{ $vcfs{$sample}{$chrom}{$pos} }[2]; #ALT allele
            }
            #Make sure all SNP positions are in all samples
            #Fill with reference genome allele information
            else
            {
                #Fill with reference genome allele information
                #print ("$sample, $chrom, $pos\n");
                $allele = substr( @{ $ref{$chrom} }[0], $pos-1, 1); #or die "$sample, $chrom, $pos";
            }
            #print "$sample, $chrom, $pos\n";
            push ( @{ $fastas{$sample}{$chrom}{$pos} }, $allele);
            push ( @{ $counts{$chrom}{$pos} }, $allele) unless (grep {$_ eq $allele} @{ $counts{$chrom}{$pos} } );
        }
    }
}

my %informativePos;

foreach my $sample ( sort keys %fastas)
{
    foreach my $chrom (sort keys %{ $fastas{$sample} } )
    {
        foreach my $pos ( sort { $fastas{$sample}{$chrom}{$a} <=> $fastas{$sample}{$chrom}{$b} } keys %{ $fastas{$sample}{$chrom} } ) #sort numerically by position
        {
            if ( scalar( @{ $counts{$chrom}{$pos} } ) > 1 )
            {
                push ( @{ $informativePos{$sample}{$chrom}{$pos} },  @{ $fastas{$sample}{$chrom}{$pos} });
            }
        }
    }
}


#################
#               #
#   SNP count   #
#               #
#################


#Count report, append mode because file already exists
open(my $countOutFH, ">>", $countOut) or die "Error writing to output file $countOut: $!\n";

my $filteredCount = 0; #= keys %{ $fastas{$s}{$c} };
my $informativeCount = 0; # = keys %{ $informativePos{$s}{$c} };


foreach my $chrom ( sort keys %{ $fastas{$s} } )
{
    $filteredCount += keys %{ $fastas{$s}{$chrom} };
    $informativeCount += keys %{ $informativePos{$s}{$chrom} };
}


#print to file
print ($countOutFH "Total filtered SNPs: $filteredCount\n");
print ($countOutFH "Total informative SNPs: $informativeCount\n\n");

#print to screen
print ("Total filtered SNPs: $filteredCount\n");
print ("Total informative SNPs: $informativeCount\n\n");

#Close count report file
close($countOutFH);


#################
#               #
#   Fasta out   #
#               #
#################

my %posList;
foreach my $chrom ( sort keys %{ $fastas{$s} } )
{
    foreach my $pos (sort keys %{ $informativePos{$s}{$chrom} })
    {
        push ( @{ $posList{$chrom} }, $pos);
    }
}

my %posListSorted;
foreach my $chrom ( sort keys %{ $fastas{$s} } )
{
    push ( @{ $posListSorted{$chrom} }, sort { $a <=> $b } @{ $posList{$chrom} });
}



#Output parsimony SNPs concatenated in table and fasta format


#Create output folder
mkdir ($fastaOutFolder) or die "Cannot create directory $fastaOutFolder: $!\n" unless -d $fastaOutFolder;

#print fasta files
my $id;
foreach my $sample (sort keys %informativePos)
{
    # print "$sample:\n";

    $id = $sample; #to save one last id to get the "root sequence". Which one doesn't matter because they all have the same positions.
    my $outName = "$fastaOutFolder/$sample.fas";
    
    open (my $fastaFH, ">", $outName)  or die "Error writing to output file $outName: $!\n";
    print ($fastaFH ">" . $sample . "\n"); #fasta header
    
    foreach my $chrom (sort keys %{ $informativePos{$sample} })
    {
        # print "\t$chrom:\n";
        foreach my $pos ( @{ $posListSorted{$chrom} } )
        # foreach my $pos (sort { $informativePos{$sample}{$chrom}{$a} <=> $informativePos{$sample}{$chrom}{$b} } keys %{ $informativePos{$sample}{$chrom} }) #sort numerically by position for output
        {
            print ($fastaFH join("", @{ $informativePos{$sample}{$chrom}{$pos} }));
            # print "\t\t$pos\n";
        }
    }
    print ($fastaFH "\n");
    
    close ($fastaFH);
}


################
#              #
#   Root seq   #
#              #
################


#Create "root" sequence
#Made from REF for all informative positions

my $rootName = "$fastaOutFolder/root.fas";
open (my $fastaFH, ">", $rootName)  or die "Error writing to output file $rootName: $!\n";
print ($fastaFH ">root" . "\n"); #fasta header

my %finalPos = %{ $informativePos{$id} };


# print "Root sequence:\n";
foreach my $chrom (sort keys %finalPos)
{
    # print "\t$chrom:\n";
    # foreach my $pos (sort { $finalPos{$chrom}{$a} <=> $finalPos{$chrom}{$b} } keys %{ $finalPos{$chrom} }) #sort numerically by position for output
    foreach my $pos ( @{ $posListSorted{$chrom} } )
    {
        # print "\t\t$pos\n";
        my $allele = substr( @{ $ref{$chrom} }[0], $pos-1, 1);
        print ($fastaFH $allele);
    }
}

print ($fastaFH "\n");
close ($fastaFH);


#################
#               #
#   SNP table   #
#               #
#################


#Create SNP table output file
open (my $tableFH, ">", $fastaTable)  or die "Error writing to output file $fastaTable: $!\n";

#referecen_pos
print ($tableFH "reference_pos");

foreach my $chrom (sort keys %finalPos)
{
    # foreach my $pos (sort { $finalPos{$chrom}{$a} <=> $finalPos{$chrom}{$b} } keys %{ $finalPos{$chrom} }) #sort numerically by position for output
    foreach my $pos ( @{ $posListSorted{$chrom} } )
    {
        print ($tableFH "\t$chrom-$pos");
    }
}

#reference_call
print ($tableFH "\nreference_call");

foreach my $chrom (sort keys %finalPos)
{
    # foreach my $pos (sort { $finalPos{$chrom}{$a} <=> $finalPos{$chrom}{$b} } keys %{ $finalPos{$chrom} }) #sort numerically by position for output
    foreach my $pos ( @{ $posListSorted{$chrom} } )
    {
        my $allele = substr( @{ $ref{$chrom} }[0], $pos-1, 1);
        print ($tableFH "\t$allele");
    }
}

#samples
foreach my $sample (sort keys %informativePos)
{
    print ($tableFH "\n$sample");
    foreach my $chrom (sort keys %{ $informativePos{$sample} })
    {
        # foreach my $pos (sort { $informativePos{$sample}{$chrom}{$a} <=> $informativePos{$sample}{$chrom}{$b} } keys %{ $informativePos{$sample}{$chrom} }) #sort numerically by position for output
        foreach my $pos ( @{ $posListSorted{$chrom} } )
        {
            print ($tableFH "\t@{ $informativePos{$sample}{$chrom}{$pos} }");
        }
    }
}

print ($tableFH "\n");
close ($tableFH);
