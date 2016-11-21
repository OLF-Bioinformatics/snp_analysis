#!/usr/bin/env perl

use strict;
use warnings;
use diagnostics;
use File::Copy;
use Parallel::ForkManager;
use Sys::CPU;
use Data::Dumper;

#I/O
my $definingSNPs = $ARGV[0];  # /home/bioinfo/prog/snp_analysis/script_dependents/Mycobacterium_bovis/DefiningSNPsGroupDesignations.txt
my $vcfFolder = $ARGV[1];  # /home/bioinfo/analyses/mbovis_script2_All/
my $report = $ARGV[2];  # /home/bioinfo/Desktop/section3.txt
my $minQual = $ARGV[3];  # 300


#Check if 4 arguments
if (scalar(@ARGV) != 4)
{
   print "Usage: perl groupFinder.pl <definingSnps.tsv> <vcfFolder> <section3.txt> <minQual>\n";
   exit;
}


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
    my ($group, $chrom, $pos, $type) = (split(/\t/, $line))[0..2, 4];
    push( @{ $defining{$type}{$chrom}{$pos} }, $group);
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

#process every VCF file in a row and store content in hash
my  %ac2InDefining;

#setting up the forking process
my $nCPU = Sys::CPU::cpu_count();
my $pm1 = Parallel::ForkManager -> new($nCPU);


$pm1->run_on_finish(sub {
    my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;

    %ac2InDefining = (%ac2InDefining, %{ $data_structure_reference->{ac2} });
});

foreach my $handle (@vcfFilesFH)
{
    my $pid = $pm1->start and next;

    my ($sample, %a2);

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
        my ($CHROM, $POS, $QUAL, $INFO) = @fields[0,1,5,7];
        my $AC = (split(/;/, $INFO))[0]; #AC is always the first field


        #AC=2 postions which are Defining SNPs and found in VCF file
        #push( @{ $defining{$type}{$chrom}{$pos} }, $group);
        foreach my $type (sort keys %defining)
        {
            #if (grep { $POS eq $_ } @{ $defining{$type}{$CHROM}{$POS} } && $AC eq 'AC=2')
            if ($defining{$type}{$CHROM}{$POS} && $AC eq 'AC=2' && $QUAL > $minQual)
            {
                #Put data in hash of hashes of array
                push(@{ $a2{$sample}{$type} }, @{ $defining{$type}{$CHROM}{$POS}});
            }
        }
    }
    $pm1 -> finish(0, { ac2 => \%a2 } );
}

$pm1 -> wait_all_children();


#close VCF files handles
foreach my $handle (@vcfFilesFH)
{
    close($handle);
}


################################
#                              #
#   Report and Copy to folder  #
#                              #
################################


#Create output file
open (my $reportFH, '>', $report)  or die "Error writing to output file $report: $!\n";

#Write header to report file
print ($reportFH "Sample\tGroup\tSubgroup\tClade\n");

#Array to keep print order for report
my @order = qw( Group Subgroup Clade ); #qw means "quote on whitespace into a list"

#Carefull, a sample can belong to more than one group or can belong to no group!

#loop though the hash of hashes of array
#push(@{ $ac2InDefining{$sample}{$type} }, @{ $defining{$type}{$CHROM}{$POS}});
foreach my $sample (sort keys %ac2InDefining)
{
    # print to report
    print ($reportFH "$sample");

    foreach my $type (@order)
    {
        if ($ac2InDefining{$sample}{$type}) #Check if it exists first
        {
            foreach my $id (@{ $ac2InDefining{$sample}{$type} })
            {
                my $folder = "$vcfFolder/$type-$id";
                mkdir ($folder) or die "Cannot create directory $folder : $!\n" unless -d $folder;

                #copy("sourcefile","destinationfile") or die "Copy failed: $!";
                my $vcf = "$vcfFolder/$sample.SNPsZeroCoverage.vcf";
                copy($vcf, $folder) or die "Copy of $sample failed: $!";

                # print to report
                print ($reportFH "\t$id");
            }
        }
        else #if not
        {
            #Print empty field
            print ($reportFH "\t");
        }
    }

    print ($reportFH "\n");
}

close ($reportFH);
