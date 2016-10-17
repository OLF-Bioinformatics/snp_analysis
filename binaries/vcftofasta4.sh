#!/bin/bash

: <<'END'
This script is the second script in a two script workflow.
Script 2 genotypes Mycobacterium tuberculosis complex and Brucella species from SNP data contained in VCFs. 
It operates on VCFs generated with the same reference output from script 1.
VCFs are collected into a single working directory.
Comparisons are output as SNP tables and alignment FASTA files to view as trees in your program of choice.

Script 2 will run and output tables and alignments from just the data generated from script 1,
however the data will be more informative if additional files are provide.

Those files are:
1) A file that contains positions to cluster individual isolates/VCFs into groups, subgroups and clades.
2) Files that contain positions to remove from the analysis.

Paradigm:
1) Once a SNP occurs and establishes in a population it does not revert back
2) Observations of homoplasy are rare
3) Group, subgroup and clade clusters only show parsimony informative SNPs for the isolates within that cluster
4) SNPs observed in a single isolate are less informative than SNPs seen in multiple isolates and therefore established in a population

Workflow summary -->

Available options
    -c with look for positions to filter.  By default, with no -c, this will not be done.
    -m will email just "M"e
    -e will run the bovis "E"lite representative samples
    -a get "a"ll_vcf alignment table

Based on VCF reference set parameters and link file dependencies
    set file to change sample names
    set defining SNPs
    turn on or off filtering
    if reference has only one chromosome filter from Excel file
    if multiple chromosomes filter from text file
    set mininum QUAL value for selecting SNPs
    set high/low QUAL value for calling a SNP "N"
    set copy location
    set email list

File checks
    convert dos files to unix
    remove special characters to renaming samples
    test for duplicate files

Count chromosome number
    1 chromosome filter from Excel worksheet
    >2 chromosomes filter from text file

Rename files to improve tree and table readability -> #####removed, done with outside scritp#####

Look for AC=1 calls (mixed SNPs) at defining SNP positions

Remove very low QUAL SNPs

Change low QUAL SNPs to "N"

Change mix SNP calls to IUPAC nomenclature

Filter positions

Group VCF files based on defining SNPs

Select SNPs with >150/300 QUAL and AC=2 call (VCF created with ploidy set to 2)

Prevent defaulting back to reference if low quality, deletion or AC=1 call present

Make aligned FASTA and alignment table files for each group

Make trees using RAxML

Organize the SNP tables

Add Map Quality averages to SNP tables.

END


######################
#                    #
#    User Defined    #
#                    #
######################


#Where analysis will take place
# baseDir=""${HOME}"/analyses/mbovis_script2"
baseDir=""${HOME}"/analyses/script2v4_test"

#Where the VCF files are
vcfPath=""${HOME}"/Desktop/vcf_mbovisCAN"

#script dependencies (the "script_dependents" folder in the "snp_analysis" folder)
dependents=""${HOME}"/prog/snp_analysis/script_dependents"


#####################
#                   #
#   Data Stucture   #
#                   #
#####################


#time stamp of analysis
uniqdate="$(date "+%Y-%m-%dat%Hh%Mm%Ss")"
# uniqdate="2016-09-15at14h39m04s"

#Files containing positions to filter
filterdir=""${baseDir}"/"${uniqdate}"-FilterFiles"

#create folder if doesn't exist
[ -d "$filterdir" ] || mkdir -p "$filterdir"


######################
#                    #
#     Resources      #
#                    #
######################


# Computer cores to use when analyzing
cpu=$(nproc) 


###############
#             #
#    Debug    #
#             #
###############


#for debug
# alias pause='read -p ""$LINENO" Enter"'


######################
#                    #
#     Text subs      #
#                    #
######################


# Sed searches put into variables
tbNumberV="s/_.*//" #Remove all charaters at and beyond "_"
tbNumberW="s/\..*//" #Remove all charaters at and beyond "."
tbNumberOnly="s/.*\([0-9]\{2\}-[0-9,FM]\{4,6\}\).*/\1/" #Only tb Number, *laboratory specific*
dropEXT="s/\(.*\)\..*/\1/" #Just drop the extention from the file


##########################
#                        #
#   Initiating reports   #
#                        #
##########################


echo -e "****************************** START ******************************\n"

echo "Start Time: "$(date "+%F %A %H:%M:%S")"" | tee "${baseDir}"/sectiontime.txt
starttime=$(date +%s)
argUsed="$1"


######################
#                    #
#     VCF files      #
#                    #
######################


#populate baseDir folder (symbolic links)
for i in $(find "$vcfPath" -type f | grep -F "SNPsZeroCoverage.vcf"); do
    name=$(basename "$i")
    ln -s "$i" "${baseDir}"/"$name"
done


#################
#               #
#    Options    #
#               #
#################


# Set flags
# flag -c with look for positions to filter.  By default, with no -c, this will not be done.
# flag -m will email just "M"e
# flag -e will run the bovis "E"lite representative samples
# flag -a get "a"ll_vcf alignment table

cflag=
mflag=
eflag=
aflag=
while getopts 'cmea' OPTION; do
    case "$OPTION" in
        c) cflag=1
        ;;
        m) mflag=1
        ;;
        e) eflag=1
        ;;
        a) aflag=1
        ;;
        ?) echo "Invalid option: -"$OPTARG"" >&2
        ;;
    esac
done
shift $(($OPTIND - 1))


#if use the "elite" flag not with bovis
if [ "$eflag" ] && [ "$1" != "bovis" ]; then 
    echo "The \"-e\" option can only be used with \"bovis\""
    exit 1
fi


# Environment controls:

if [ "$1" = "ab1" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt" #to change the name of vcf files

    #reference genome
    genome=""${dependents}"/Brucella_abortus/NC_00693c.fasta"

    # When more than one chromosome
    # Genbank files must have "NC" file names that match NC numbers in VCF chrom identification in column 1 of vcf
    # Example: File name: NC_017250.gbk and "gi|384222553|ref|NC_017250.1|" listed in vcf
    gbk_file1=""${dependents}"/Brucella_abortus/NC_006932.gbk"
    gbk_file2=""${dependents}"/Brucella_abortus/NC_006933.gbk"
    gbk_files=("$gbk_file1" "$gbk_file2")
    # echo "${gbk_files[@]}"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_abortus/Abortus1_Defining_SNPs.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_abortus/FilterFiles" #Files containing positions to filter
    
    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using Brucella abortus bv 1, 2 or 4 variables" | tee "${baseDir}"/section5.txt

    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "mel" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_melitensis/BmelitensisM5-90.fasta"

    gbk_file1=""${dependents}"/Brucella_melitensis/NC_017246.gbk"
    gbk_file2=""${dependents}"/Brucella_melitensis/NC_017247.gbk"
    gbk_files=("$gbk_file1" "$gbk_file2")

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_melitensis/Mel_Defining_SNPs.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_melitensis/FilterFiles" #Files containing positions to filter
    
    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. melitensis variables" | tee "${baseDir}"/section5.txt

    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "suis1" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_suis_bv1/NC_01725c.fasta"

    gbk_file1=""${dependents}"/Brucella_suis_bv1/NC_017250.gbk"
    gbk_file2=""${dependents}"/Brucella_suis_bv1/NC_017251.gbk"
    gbk_files=("$gbk_file1" "$gbk_file2")

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_suis_bv1/Suis1_Defining_SNPs.txt"

    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_suis_bv1/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. suis bv1 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "suis2" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_suis_bv2/Bsuisbv2-94-11.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_suis_bv2/suis2_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_suis_bv2/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. suis bv2 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "suis3" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_suis_bv3/B-REF-BS3-686.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_suis_bv3/Suis3_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_suis_bv3/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. suis bv3 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "suis4" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_suis_bv4/B-REF-BS4-40.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_suis_bv4/Suis4_Defining_SNPs.txt"

    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_suis_bv/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. suis bv4 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "canis" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_canis/BcanisATCC23365.fasta"

    gbk_file1=""${dependents}"/Brucella_canis/NC_010103.gbk"
    gbk_file2=""${dependents}"/Brucella_canis/NC_010104.gbk"
    gbk_files=("$gbk_file1" "$gbk_file2")

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_canis/Canis_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_canis/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "vcftofasta.sh ran as B. canis"
    echo "Script vcftofasta.sh ran using B. canis variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "ceti1" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_ceti-grp1/Bceti1Cudo.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_ceti-grp1/Ceti1_Defining_SNPs.txt"

    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_ceti-grp1/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B ceti group 1 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "ceti2" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_ceti-grp2/BBceti2-TE10759.fasta"

    gbk_file1=""${dependents}"/Brucella_abortus/NC_022905.gbk"
    gbk_file2=""${dependents}"/Brucella_abortus/NC_022906.gbk"
    gbk_files=("$gbk_file1" "$gbk_file2")

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_ceti-grp2/Ceti2_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_ceti-grp2/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B ceti group 2 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

elif [ "$1" = "ovis" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"

    #reference genome
    genome=""${dependents}"/Brucella_ovis/BovisATCC25840.fasta"

    gbk_file1=""${dependents}"/Brucella_ovis/NC_009505.gbk"
    gbk_file2=""${dependents}"/Brucella_ovis/NC_009504.gbk"
    gbk_files=("$gbk_file1" "$gbk_file2")

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_ovis/Ovis_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    filterdir=""${dependents}"/Brucella_ovis/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. ovis variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"

### TB ###

elif [ "$1" = "bovis" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/Mycobacterium_bovis/NC_002945.fasta"

    gbk_files=(""${dependents}"/Mycobacterium_bovis/NC_002945.gbk")

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Mycobacterium_bovis/DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using M. bovis variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    if [ "$eflag" ]; then
        echo "Only the "elite" bovis isolates are being ran"
    else
        echo "All bovis are being ran"
    fi

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/Mycobacterium_bovis/Filtered_Regions.xlsx"

    #Usage: perl parseXL.pl <input.xlsx> <output.tsv>
    #Usage: perl rangeExpander.pl <filterFile.txt> <output_folder>
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

elif [ "$1" = "tb1" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/TB1/NC_017528.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB1/tb1DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB1/tb1Filtered_Regions.xlsx"
    
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

elif [ "$1" = "tb2" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/TB2/NC_021251.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB2/tb2DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB2/tb2Filtered_Regions.xlsx"
    
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

elif [ "$1" = "tb3" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/TB3/NC_021193.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB3/tb3DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB3/tb3Filtered_Regions.xlsx"
    
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

elif [ "$1" = "tb4a" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/TB4a/NC002755.fasta"

    # gbk_file=""${dependents}"/NC_018143.gbk"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB4a/tb4aDefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB4a/tb4aFiltered_Regions.xlsx"
    
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

elif [ "$1" = "tb4b" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/TB4b/NC018143.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB4b/tb4bDefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB4b/tb4bFiltered_Regions.xlsx"
    
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

elif [ "$1" = "tb5" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/TB5/APKD01000001.fasta"

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB5/tb5DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB5/tb5Filtered_Regions.xlsx"
    
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

elif [ "$1" = "tb6" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/TB6/NC_015758.fasta"

    gbk_files=(""${dependents}"/TB6/NC_015758.gbk")

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB6/tb6DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB6/tb6Filtered_Regions.xlsx"
    
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

elif [ "$1" = "para" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"

    #reference genome
    genome=""${dependents}"/paraTB/NC_002944.fasta"

    gbk_files=(""${dependents}"/paraTB/NC_002944.gbk")

    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/paraTB/DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using para variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/paraTB/Filtered_Regions.xlsx"
    
    parseXLSX.pl \
        "$excelinfile" \
        /dev/stdout \
        | rangeExpander.pl \
            /dev/stdin \
            "$filterdir"

else

    echo -e "\
Usage: vcftofasta.sh [flag] <species>

Madatory argument (species; chose only one):
    ab1    Brucella abortus bv 1
    mel    Brucella melitensis
    suis1  Brucella suis bv 1
    suis2  Brucella suis bv 2
    suis3  Brucella suis bv 3
    suis4  Brucella suis bv 4
    canis  Brucella canis
    ceti1  Brucella ceti bv 1
    ceti2  Brucella ceti bv 2
    ovis   Brucella ovis
    bovis  Mycobacterium bovis (TBBOV)
    tb1    Mycobacterium bovis group 1
    tb2    Mycobacterium bovis group 2
    tb3    Mycobacterium bovis group 3
    tb4a   Mycobacterium bovis group 4a
    tb4b   Mycobacterium bovis group 4b
    tb5    Mycobacterium bovis group 5
    tb6    Mycobacterium bovis group 6
    para   Mycobacterium avium subsp. paratuberculosis

Optional flags:

    -c     look for positions to filter.  By default, with no -c, this will not be done.
    -m     email just "M"e
    -e     run the bovis "E"lite representative samples
    -a     get "a"ll_vcf alignment table
"

    rm -rf "$baseDir" #remove analysis folder
    exit 1 #stop script execution
fi



function removeIsolates ()
{
    if [ -n "$RemoveFromAnalysis" ]; then
        echo "Removing unwanted isolates"
        cat "$RemoveFromAnalysis" \
            | tr '\r' '\n' \
            | awk '{print $1}' \
            > "${baseDir}"/RemoveFromAnalysisUnixReady.txt

        removeList=$(cat "${baseDir}"/RemoveFromAnalysisUnixReady.txt)

        for i in "$removeList"; do
            rm *"${i}"* > /dev/null 2>&1
        done

        rm "${baseDir}"/RemoveFromAnalysisUnixReady.txt
    fi
}


# If there are 2 vcf files with the same name one of the files might unknowingly
# get cut out of the analysis and keep the undesired vcf instead.  This will
# alert if 2 vcf with the same TB number are present.
# The regular expression used in sed should be changed based on vcf naming convention
function testDuplicates ()
{
    echo -e "\nChecking for empty or duplicated VCF files."

    directorytest="${baseDir##*/}" #name of directory where script was launched (or name of "$baseDir")
    if [ "$directorytest" = "VCF_Source_All" ]; then # if where all the reference vcf files are stored
        echo "Change directory name and restart"
        exit 1
    fi

    for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
        if [ ! -s "$i" ]; then
            echo ""$i" is empty.  Fix and restart script"
            # exit 1
        fi

        getbase=$(basename "$i")
        number=$(echo "$getbase" | sed "$tbNumberV" | sed "$tbNumberW")
        echo "$number" >> "${baseDir}"/list.txt
    done

    duplist=$(cat "${baseDir}"/list.txt | sort | uniq -d)
    dupNumberSize=$(echo "$duplist" | wc | awk '{print $3}')

    rm "${baseDir}"/list.txt

    if [ "$dupNumberSize" -gt 4 ]; then
        echo "There are duplicated VCF files."
        echo "Please remove the duplicates and restart script."
        echo "$duplist"
        exit 1 # Error status
    else
        echo "Good! No duplicated VCF files present"
    fi
}


#################################################################################

# This function prepares the filter files.
# awk needs to see a number in the file, so if the file is blank 2 matching numbers are added.  2 numbers because duplicates are kept therefore allowing a number to be pasting into awk when comparing files.

function filterFilespreparation ()
{
    # For tb, inputXLS.py creates text files with positions to be filetered, and places them in filterdir
    echo "Waiting for filter file creation to complete"

    echo "Preparing Filter Files"
    for i in $(find -L "$filterdir" -type f | grep -F ".txt"); do
        getbase=$(basename "$i")
        number=$(echo "$getbase" | sed 's/\(.*\)\..*/\1/')

        cat "$i" \
            | sort | uniq \
            > "${baseDir}"/"${number}".num

        if [ "$chromCount" -eq 1 ]; then
            echo "100000000" >> "${baseDir}"/"${number}".num
            echo "100000000" >> "${baseDir}"/"${number}".num
        elif [ "$chromCount" -eq 2 ]; then
            echo "chrom1    100000000" >> "${baseDir}"/"${number}".num
            echo "chrom1    100000000" >> "${baseDir}"/"${number}".num
            echo "chrom2    100000000" >> "${baseDir}"/"${number}".num
            echo "chrom2    100000000" >> "${baseDir}"/"${number}".num
        else
            echo "Greater than 2 chromosomes present."
        fi

        rm "$i"
        mv "${baseDir}"/"${number}".num "${baseDir}"/"${number}".txt
    done

    echo "Finished preparing filter files"
}


#################################################################################

# Change SNPs with low QUAL values to N, based on parameter set above in variable settings
function changeLowCalls ()
{
    echo -e "\nChanging low calls, started --> "$(date "+%F %A %H:%M:%S")""
    for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
        base=$(basename "$i")
        baseNoExt="${base%.*}"

        cat "$i" \
            | awk -v x="$QUAL" -v y="$highEnd" 'BEGIN {OFS="\t"} { if ($6 >= x && $6 <= y) print $1, $2, $3, $4, "N", $6, $7, $8, $9, $10; else print $0 }' \
            > "${baseDir}"/"${baseNoExt}".txt

        rm "$i"
        mv "${baseDir}"/"${baseNoExt}".txt "${baseDir}"/"${baseNoExt}".vcf
    done
}

#################################################################################

function findpositionstofilter ()
{
    echo "Finding positions to filter --> "$(date "+%F %A %H:%M:%S")""

    d="$1"

    # positions have already been filtered via cutting specific positions.
    cp "${d}"/filtered_total_ref_pos.txt "${d}"/total_pos.txt
    cat "${d}"/total_pos.txt | awk '{print $1}' > "${d}"/prepositionlist.txt

    for n in $(cat "${d}"/prepositionlist.txt); do
        front=$(echo "$n" | sed 's/\(.*\)-\([0-9]*\)/\1/')
        back=$(echo "$n" | sed 's/\(.*\)-\([0-9]*\)/\2/')
        # echo "front: $front"
        # echo "back: $back"

        positioncount=$(cat $(find "${d}" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" ' $1 = f && $2 = b {count++} END {print count}')

        #echo "position count: $positioncount"
        if [ "$positioncount" -gt 2 ]; then
            #printf "%s\t%s\n" "$front" "$back"
            echo "$n" >> "${d}"/positionlist.txt
        else
            echo "$n" >> "${d}"-DONOT_filtertheseposition.txt # "$d" from fasta_table ()
        fi
    done

    echo "Filtering --> "$(date "+%F %A %H:%M:%S")""

    for p in $(cat "${d}"/positionlist.txt); do
        front=$(echo "$p" | sed 's/\(.*\)-\([0-9]*\)/\1/')
        back=$(echo "$p" | sed 's/\(.*\)-\([0-9]*\)/\2/')
        #echo "front: $front"
        #echo "back: $back"

        maxqual=$(cat $(find "$d" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" 'BEGIN{max=0} $1 = f && $2 = b {if ($6>max) max=$6} END {print max}' \
            | sed 's/\..*//')

        avequal=$(cat $(find "${d}" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" '$6 != "." && $1 = f && $2 = b {print $6}' \
            | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' \
            | sed 's/\..*//')

        maxmap=$(cat $(find "$d" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" ' $1 = f && $2 = b {print $8}' \
            | sed 's/.*MQ=\(.....\).*/\1/' | awk 'BEGIN{max=0}{if ($1>max) max=$1} END {print max}' \
            | sed 's/\..*//')

        avemap=$(cat $(find "$d" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" '$6 != "." && $1 = f && $2 = b {print $8}' \
            | sed 's/.*MQ=\(.....\).*/\1/' \
            | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' \
            | sed 's/\..*//')

        #change maxmap from 52 to 56 2015-09-18
        if [ "$maxqual" -lt 1300  ] || [ "$avequal" -lt 800 ] || [ "$maxmap" -lt 58  ] || [ "$avemap" -lt 57 ]; then
            echo "maxqual "$maxqual"" >> "${d}"/filterpositiondetail.txt
            echo "avequal "$avequal"" >> "${d}"/filterpositiondetail.txt
            echo "maxmap "$maxmap"" >> "${d}"/filterpositiondetail.txt
            echo "avemap "$avemap"" >> "${d}"/filterpositiondetail.txt
            echo "position "$p"" >> "${d}"/filterpositiondetail.txt
            echo ""  >> "${d}"/filterpositiondetail.txt
            echo "$p" >> "${d}"-filtertheseposition.txt
        else
            echo "$p" >> "${d}"-DONOT_filtertheseposition.txt
            #echo "maxqual $maxqual"
            #echo "maxmap $maxmap"
            #echo "avemap $avemap"
            #echo "position $p"
            #echo ""
        fi
    done

    #cleanup
    rm "${d}"/positionlist.txt
    rm "${d}"/prepositionlist.txt
    rm "${d}"/total_pos.txt
}

#################################################################################

#   Function: fasta and table creation
function fasta_table ()
{
    # Loop through the directories
    if [ -n "$eflag" ] || [ -n "$aflag" ] && [ "$1" = ""${baseDir}"/all_vcfs" ]; then
        directories=("$1")
    else
        # directories=($(find ""${baseDir}"/all_subgroups" -maxdepth 1 -type d | awk 'NR > 1'))
        directories=($(find "$1" -maxdepth 1 -type d | awk 'NR > 1'))
    fi

    wait

    # echo "${directories[@]}" | tr " " "\n"

    for d in "${directories[@]}"; do
        # d=""${baseDir}"/all_groups/Group-9"
        dName=$(basename "$d")

        #backup vcf files
        [ -d "${d}"/starting_files ] || mkdir "${d}"/starting_files
        cp "${d}"/*.vcf "${d}"/starting_files

        if [ "$d" != ""${baseDir}"/all_vcfs" ]; then #skip is in all_vcfs folder, the filtering has been done already
            #Filterout group/subgroup/clade specific positions
            if [ "$FilterGroups" = "yes" ]; then

                echo -e "\nFiltering positions for "$dName":"

                #Mark vcf allowing areas of the genome to be removed from the SNP analysis
                for i in $(find "$d" -maxdepth 1 -type f | grep -F ".vcf"); do #don't look at the backed up vcf folder (starting_files)
                    # i=""${d}"/MBWGS083.SNPsZeroCoverage.vcf"
                    m=$(basename "$i")
                    # n=$(echo "$m" | sed "$dropEXT")
                    n="${m%%.*}"

                    echo -e "Working on "$n"..."

                    #Usage: perl regionRemover.pl <FilterToAll.txt> <chroms.txt> <input.vcf> <filtered.vcf>
                    regionRemover.pl \
                        "${filterdir}/FilterToAll.txt" \
                        "${baseDir}"/chroms.txt \
                        "$i" \
                        "${i}".filtered
                    wait

                    #Replace original vcf by the filtered one
                    mv "${i}".filtered "$i"
                done
            fi
        fi

        echo ""${dName}":" >> "${baseDir}"/section4.txt

        #Usage: perl snpTableMaker.pl <ref.fasta> <vcfFolder> <minQual> <minAltQual> <AC1Report.txt> <section4.txt> <fastaOutFolder> <fastaTable.tsv>
        snpTableMaker.pl \
            "$genome" \
            "$d" \
            "$QUAL" \
            "$highEnd" \
            "${d}"/"${dName}"_AC1Report.txt \
            "${baseDir}"/section4.txt \
            "${d}"/fasta \
            "${d}"/"${dName}".table.txt #> "${baseDir}"/"${dName}"_outputPostions.txt
        wait

        # Make a file containing all fasta files. Used awk instead of cat to insure newline between files
        cat "${d}"/fasta/*.fas > "${d}"/"${dName}"_alignment.fasta
    done
}


function alignTable ()
{
    d="$1"
    # d=""${baseDir}"/all_groups/Group-9/fasta"
    # d=""${baseDir}"/all_vcfs/fasta"

    #Group/Subgroup/Clade name
    parent=$(dirname "$d")
    dName="$(basename "$parent")"
    

    ###############
    #             #
    #    RAxML    #
    #             #
    ###############


    # Beginning in fasta folder
    echo -e "\nRAxML started on "$dName" --> "$(date "+%F %A %H:%M:%S")""

    cat "${d}"/*.fas \
        | awk '{print $0}'  \
        | sed '/root/{N;d;}' \
        >> "${d}"/fastaGroup.txt

    cat "${d}"/*.fas \
        | awk '{print $0}' \
        >> "${d}"/RAxMLfastaGroup.txt

    #Make the tree
    raxmlHPC-PTHREADS-AVX2 \
        -T "$cpu" \
        -s "${d}"/RAxMLfastaGroup.txt \
        -w "$d" \
        -n "$dName" \
        -m GTRCAT \
        -p "$RANDOM" \
        &>/dev/null
    wait

    # reroot tree
    # don't save to file because the trees are cut to the right
    # nw_reroot \
    #     "${d}"/RAxML_bestTree."$dName" \
    #     root \
    #     | tee >(tee "${d}"/tableinput."$dName" "${d}"/rooted_RAxML_bestTree."$dName" 2&> /dev/null) \
    #     | tee >(nw_display \
    #             -s \
    #             -w 1500 \
    #             -v 20 \
    #             -b 'opacity:0' \
    #             -i 'font-size:8' \
    #             -l 'font-family:serif;font-style:italic' \
    #             -d 'stroke-width:2;stroke:blue' \
    #             /dev/stdin \
    #                 | tee >(tee "${d}"/"${dName}"-tree.svg 2&> /dev/null) \
    #                 | tee >(rsvg-convert \
    #                         -f pdf \
    #                         "${d}"/"${dName}"-tree.svg \
    
    # reroot tree
    nw_reroot \
        "${d}"/RAxML_bestTree."$dName" \
        root \
        | tee "${d}"/tableinput."$dName" "${d}"/rooted_RAxML_bestTree."$dName" &>/dev/null
    wait

    mv "${d}"/rooted_RAxML_bestTree."$dName" "${d}"/RAxML_bestTree."$dName"

    #cleanup if exsists
    [ -f "${d}"/RAxML_parsimonyTree* ] && rm "${d}"/RAxML_parsimonyTree*

    cat "${d}"/tableinput."$dName" \
        | tr ":" "\n" \
        | tr "," "\n" \
        | sed -e 's/(//g' -e 's/)//g' \
        | grep -v "\.[0-9]*" \
        | grep -v "[0-9]e-[0-9]" \
        | grep -v "root" \
        > "${d}"/cleanedAlignment.txt


    # Place headers onto aligned file
    (echo "reference_call"; cat "${d}"/cleanedAlignment.txt) \
        > "${d}"/cleanedAlignment.txt.temp
    mv "${d}"/cleanedAlignment.txt{.temp,}
    wait

    (echo "reference_pos"; cat "${d}"/cleanedAlignment.txt) \
        > "${d}"/cleanedAlignment.txt.temp
    mv "${d}"/cleanedAlignment.txt{.temp,}
    wait

    #Sort and organize SNP table
    sortOrganizeTable.py \
        "${parent}"/"${dName}".table.txt \
        "${d}"/cleanedAlignment.txt \
        "${d}"/"${dName}".sorted_table.txt \
        "${d}"/"${dName}".organized_table.txt
    wait

    #replace sorted table
    mv "${d}"/"${dName}".organized_table.txt "${d}"/"${dName}".sorted_table.txt



    # Add map qualities to sorted table
    echo "Adding map qualities..."

    # Get just the position from first line.  The chromosome must be removed
    cat "${d}"/"${dName}".sorted_table.txt \
        | head -n 1 \
        | tr "\t" "\n" \
        | sed "1d" \
        | awk '{print NR, $0}' \
        > "${d}"/"${dName}"-positions.txt


    echo -e "Sorting mapping quality table..."

    #Usage: perl qualityExtractor.pl <position.txt> <vcf_folder> <quality.txt>
    qualityExtractor.pl \
        "${d}"/"${dName}"-positions.txt \
        "${parent}"/starting_files \
        "${d}"/quality.txt
    wait

    mapvalues.py \
        "${d}"/"${dName}".sorted_table.txt \
        "${d}"/quality.txt \
        "${d}"/"${dName}".transposed_table.txt \
        "${d}"/"${dName}".finished_table.txt
    wait

    #Overwrite the old table
    mv "${d}"/"${dName}".finished_table.txt "${d}"/"${dName}".sorted_table.txt



    echo "Organizing mappping quality table..."

    #run python script
    mapvalues.py \
        "${d}"/"${dName}".sorted_table.txt \
        "${d}"/quality.txt \
        "${d}"/"${dName}".transposed_table.txt \
        "${d}"/"${dName}".finished_table.txt
    wait

    #Overwrite the old table
    mv "${d}"/"${dName}".finished_table.txt "${d}"/"${dName}".organized_table.txt

    #Cleanup
    rm "${d}"/quality.txt
    rm "${d}"/"${dName}".transposed_table.txt
    rm -r "${parent}"/starting_files
}





#################################################################################
#################################################################################
#################################################################################
###################################### START ####################################
#################################################################################
#################################################################################
#################################################################################


# Clean the tag file that has been exported to Desktop
#chmod 777 ${genotypingcodes}  
#cat ${genotypingcodes} | tr '\r' '\n' | awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' | sed 's/\"##/##/' | sed 's/MN_Wildlife_Deer_//' > preparedTags.txt

#clean_tag.sh $genotypingcodes
####################
# Clean the genotyping codes used for naming output
#sed 's/\*//g' < preparedTags.txt -e 's/(/_/g' -e 's/)/_/g' -e 's/ /_/g' -e 's/-_/_/g' -e 's/\?//g' -e 's/_-/_/g' -e 's/,/_/g' -e 's#/#_#g' \
#| sed 's#\\#_#g' | sed 's/__/_/g' | sed 's/__/_/g' | sed 's/__/_/g' | sed 's/-$//g' | sed 's/_$//g' |awk 'BEGIN {OFS="\t"}{gsub("_$","",$1)}1' > "${baseDir}"/outfile.txt
#rm preparedTags.txt

#cat ${genotypingcodes} | tr '\r' '\n' | grep "Yes" | sed 's/_.*//' >> elite
#echo "Only samples in this file will be ran when elite is used as the secound argument" >> elite

####################

# Test for duplicate VCFs
testDuplicates
wait

#Prepare Filter files.
# filterFilespreparation
# wait

#Test for match coverage file
# checkMatchingCoverageFile

# if [ "$eflag" ]; then
#     echo "Only analyzing elite files"

#     for i in $(cat elite); do
#         name=$(ls starting_files | grep "$i")
#         cp ./starting_files/$name ./
#     done

#     for i in $(find ./starting_files/ -mtime -15); do
#         cp "$i" ./
#     done

# else
#     echo "all samples will be ran"
#     cp ./starting_files/* ./
# fi

# rm elite


#################################################################################


# Count the number of chromosomes used in the reference when VCFs were made.
#singleFile=`ls *.vcf | head -1`
echo -e "\nCounting the number of chromosomes in first 100 samples, started -->  "$(date "+%F %A %H:%M:%S")""
chromCount=$(cat $(find -L "$baseDir" -type f | grep -F ".vcf" | head -100) \
                | awk ' $0 !~ /^#/ {print $1}' \
                | sort | uniq -d \
                | awk 'END {print NR}')

echo "The number of chromosomes/segments seen in VCF files: "$chromCount""
cat $(find -L "$baseDir" -type f | grep -F ".vcf" | head -100) \
    | awk ' $0 !~ /^#/ {print $1}' \
    | sort | uniq -d \
    > "${baseDir}"/chroms.txt

#Chromosomes/segments found
echo -e "$(cat "${baseDir}"/chroms.txt)"


#################################################################################

# Remove selected isolates from comparison
# This is optional, and should be turned on or off based on laboratories preference
# removeIsolates

#################################################################################


#######################
#                     #
#    File renaming    #
#                     #
#######################


#This is optional

for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
    base=$(basename "$i")
    searchName=$(cut -d "." -f 1 <<< "$base")
    # echo "Original File: "$base""
    # echo "searchName: "$searchName""

    # Direct script to text file containing a "${baseDir}"/list.txt of the correct labels to use.
    # The file must be a txt file.
    if [ -f "${baseDir}"/outfile.txt ]; then
        p=$(cat "${baseDir}"/outfile.txt | grep "$searchName")
        echo "This is what was found in tag file: "$p""
        newName=$(echo "$p" | awk '{print $1}' | tr -d "[:space:]") # Captured the new name
        n=$(echo "$base" | sed "$tbNumberV" | sed "$tbNumberW")
        noExtention=$(echo "$base" | sed "$dropEXT")
        VALtest=$(echo "$i" | grep "VAL")
        # echo "VALtest: $VALtest"
        # h=`echo ${i%-AZ}`; g=`echo ${h%-Broad}`; echo $g

        #Check if a name was found in the tag file.  If no name was found, keep original name, make note in log and cp file to unnamed folder.
        if [ -z "$p" ]; then # new name was NOT found
            if [ -z "$VALtest" ]; then
                name="$searchName"
                echo "n is "$n""
                echo "$name" >> "${baseDir}"/section1.txt
                
                mkdir -p "${baseDir}"/FilesNotRenamed
                cp "$i" "${baseDir}"/FilesNotRenamed
                mv "$i" "${baseDir}"/"${name}".vcf
                # echo "A"
            else
                name="${searchName}"-Val
                mv "$i" "${baseDir}"/"${name}".vcf
                # echo "B"
            fi
        else # New name WAS found
            if [ -z "$VALtest" ]; then
                name="$newName"
                mv "$i" "${baseDir}"/"${name}".vcf
                # echo "C"
            else
                name="${newName}"-Val
                echo "newName is $name"
                mv "$i" "${baseDir}"/"${name}".vcf
                # echo "D"
            fi
        fi
    fi
done

[ -e "${baseDir}"/outfile.txt ] && rm "${baseDir}"/outfile.txt

##################### Start: Make Files Unix Compatiable #####################

#Fix validated (VAL) vcf files.  This is used in vcftofasta scripts to prepare validated vcf files opened and saved in Excel.
#Create "${baseDir}"/list.txt of isolates containing "VAL"
#Do NOT make this a child process.  It messes changing column 1 to chrom


#####################
#                   #
#    Dos to Unix    #
#                   #
#####################


#check if file is in dos format
isDos=()
for j in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
    dosLine=$(cat "$j" | grep -IU --color "^M")
    # [ -n "$dosLine" ] && isDos+=("$j")
done

wait

#if some files have windows-type carriage returns
if [ "${#isDos[@]}" -gt 0 ]; then
    echo "#############################"
    echo "Making Files Unix Compatiable"

    #Only change those files
    for v in "${isDos[@]}"; do
        dos2unix "$v" > /dev/null 2>&1 #Fixes files opened and saved in Excel
        cat "$v" | tr '\r' '\n' \
            | awk -F '\t' 'BEGIN{OFS="\t";} {gsub("\"","",$5);print;}' \
            | sed -e 's/\"##/##/' -e 's/\"AC=/AC=/' \
            > "${baseDir}"/"${v}".temp
        mv "${baseDir}"/"${v}".temp "${baseDir}"/"$v"
    done
fi


###############################
#                             #
#    AC1s in Defining SNPs    #
#                             #
###############################


# Looks for defining positions in VCF files.
# If an AC=1 is found at a defined position it is flagged as a posible mixed infection.
# These defining positions must be SNPs found cluster's main branch

#Usage: perl AConeInDefiningSNPs.pl <definingSnps.tsv> <vcfFolder> <section2.txt>
AConeInDefiningSNPs.pl \
    "$DefiningSNPs" \
    "$baseDir" \
    "${baseDir}"/section2.txt
wait




# Change low QUAL SNPs to N, see set variables
#gets rid of symbolic links and saves "real" files in "$baseDir"
changeLowCalls
wait


#######################
#                     #
#    AC1s to IUPAC    #
#                     #
#######################


echo -e "\nChanging AC=1 to IUPAC -->  "$(date "+%F %A %H:%M:%S")""

for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
    awk '
BEGIN { OFS = "\t"}

{ if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AG/ )
         print $1, $2, $3, $4, "R", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CT/ )
         print $1, $2, $3, $4, "Y", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GC/ )
         print $1, $2, $3, $4, "S", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AT/ )
         print $1, $2, $3, $4, "W", $6, $7, $8, $9, $10
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GT/ )
         print $1, $2, $3, $4, "K", $6, $7, $8, $9, $10         
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /AC/ )
         print $1, $2, $3, $4, "M", $6, $7, $8, $9, $10         
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /GA/ )
         print $1, $2, $3, $4, "R", $6, $7, $8, $9, $10         
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TC/ )
         print $1, $2, $3, $4, "Y", $6, $7, $8, $9, $10         
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CG/ )
         print $1, $2, $3, $4, "S", $6, $7, $8, $9, $10         
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TA/ )
         print $1, $2, $3, $4, "W", $6, $7, $8, $9, $10         
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /TG/ )
         print $1, $2, $3, $4, "K", $6, $7, $8, $9, $10         
else if ($6 > 300 && $8 ~ /^AC=1;/ && $4$5 ~ /CA/ )
         print $1, $2, $3, $4, "M", $6, $7, $8, $9, $10         
else
    print $0     
}' "$i" > "${i%vcf}temp"

    mv "${i%vcf}temp" "$i"
done
wait


############################################
#                                          #
#    Mark files / Remove marked regions    #
#                                          #
############################################


#Remove positions with missing genotypes and the ones to from the filter file

echo "Number of chromosomes:  "$chromCount""

if [ "$FilterAllVCFs" = "yes" ]; then

    echo -e "\nFiltering out positions from FilterToAll.txt --> "$(date "+%F %A %H:%M:%S")""

    #Label filter field for positions to be filtered in all VCFs
    if [ "$chromCount" -ge 1 ]; then

        for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
            m=$(basename "$i")
            # n=$(echo "$m" | sed "$dropEXT")
            n="${m%%.*}"

            echo -e "Working on "$n"..."

            #Usage: perl regionRemover.pl <FilterToAll.txt> <chroms.txt> <input.vcf> <filtered.vcf>
            regionRemover.pl \
                "${filterdir}/FilterToAll.txt" \
                "${baseDir}"/chroms.txt \
                "$i" \
                "${i}".filtered
            wait

            #Replace original vcf by the filtered one
            mv "${i}".filtered "$i"
            wait
        done
    else
        echo "No chromosome detected."
        exit 1
    fi
    wait

else
    echo "***VCF filtering not done. FilterAllVCFs was set to \"no\"." | tee -a "${baseDir}"/section5.txt
fi


###########################################################
#                                                         #
#    Categorize VCF files (Group, Subgroup and Clades)    #
#                                                         #
###########################################################


#Find group using DefiningSNP positions and copy files to group folders

echo -e "\nFind groups using DefiningSNP positions --> "$(date "+%F %A %H:%M:%S")""

#Usage: perl AConeInDefiningSNPs.pl <definingSnps.tsv> <vcfFolder> <section2.txt>
groupFinder.pl \
    "$DefiningSNPs" \
    "$baseDir" \
    "${baseDir}"/section3.txt \
    "$QUAL"


[ -d "${baseDir}"/all_vcfs ] || mkdir -p "${baseDir}"/all_vcfs #Make all_vcfs folder if one does not exist.
for i in $(find -L "$baseDir" -maxdepth 1 -type f | grep -F ".vcf"); do
    mv "$i" "${baseDir}"/all_vcfs
done


#########################################
#                                       #
#    Organize VCF files into folders    #
#                                       #
#########################################


[ -d "${baseDir}"/all_groups ] || mkdir -p "${baseDir}"/all_groups
mv "${baseDir}"/Group* "${baseDir}"/all_groups

[ -d "${baseDir}"/all_subgroups ] || mkdir -p "${baseDir}"/all_subgroups
mv "${baseDir}"/Subgroup*/ "${baseDir}"/all_subgroups/

[ -d "${baseDir}"/all_clades ] || mkdir -p "${baseDir}"/all_clades
mv "${baseDir}"/Clade*/ "${baseDir}"/all_clades/

wait


######################
#                    #
#    Fasta Tables    #
#                    #
######################


#all_vcfs
if [ -n "$eflag" ] || [ -n "$aflag" ]; then
    echo -e "\nCreating fasta table for all_vcfs --> "$(date "+%F %A %H:%M:%S")""
    fasta_table "${baseDir}"/all_vcfs
    wait
fi

#all_groups
if [ -n "$(ls -A "${baseDir}"/all_groups)" ]; then
    echo -e "\nCreating fasta table for all_groups --> "$(date "+%F %A %H:%M:%S")""
    fasta_table "${baseDir}"/all_groups
    wait
else
    echo "No Groups"
fi

#all_subgroups
if [ -n "$(ls -A "${baseDir}"/all_subgroups)" ]; then
    echo -e "\nCreating fasta table for all_subgroups --> "$(date "+%F %A %H:%M:%S")""
    fasta_table "${baseDir}"/all_subgroups
    wait
else
    echo "No SubGroup"
fi

#all_clades
if [ -n "$(ls -A "${baseDir}"/all_clades)" ]; then
    echo -e "\nCreating fasta table for all_clades --> "$(date "+%F %A %H:%M:%S")""
    fasta_table "${baseDir}"/all_clades
    wait
else
    echo "No Clades"
fi


#Cleanup
# rm -rf "$filterdir"

#For QA, keep a copy of these files to reproduce analysis
# cp "$DefiningSNPs" "$baseDir"
# cp "${dependents}"/Table_Template.xlsx "$baseDir"
# cp "$0" "$baseDir" #$0 is the name of the script itself


#########################
#                       #
#    Table alignment    #
#                       #
#########################


echo -e "\nOrganizing SNP tables --> "$(date "+%F %A %H:%M:%S")""

if [ -n "$eflag" ] || [ -n "$aflag" ]; then
    directories=($(find "$baseDir" -maxdepth 1 -type d \
        | awk 'NR > 1' \
        | grep -v "$uniqdate"))
else
    directories=($(find "$baseDir" -maxdepth 1 -type d \
    | awk 'NR > 1' \
    | grep -v "$uniqdate" \
    | grep -vF "all_vcfs"))
fi

wait

# echo "${directories[@]}" | tr " " "\n"

#all_groups/all_subgroups/all_clades
for d in "${directories[@]}"; do
    # echo "$d"
    if [ -n "$(find "$d" -maxdepth 1 -type d | grep -vF "starting_files\|fasta" | awk 'NR > 1')" ] && [ "$d" != ""${baseDir}"/all_vcfs" ]; then 
        subdirectories=($(find "$d" -maxdepth 1 -type d \
            | grep -vF "starting_files\|fasta" \
            | awk 'NR > 1'))

        for sd in "${subdirectories[@]}"; do
            # echo "$sd" 
            alignTable "${sd}"/fasta
            wait
        done

        wait
    elif [ "$d" = ""${baseDir}"/all_vcfs" ]; then #if all_vcfs folder, there is no subfolders for groups
        # echo "$d"
        alignTable "${d}"/fasta
        wait
    else
        echo -e "\nNothing to do for "$(basename "$d")""
    fi
done

wait


#########################
#                       #
#    Final reporting    #
#                       #
#########################


# echo -e "\nPreparing final report --> "$(date "+%F %A %H:%M:%S")""

[ -f "${baseDir}"/section1.txt ] && cat "${baseDir}"/section1.txt | column > "${baseDir}"/csection1.txt


#Script execution time
echo -e "\nEnd Time: "$(date "+%F %A %H:%M:%S")"\n" | tee -a "${baseDir}"/sectiontime.txt
endtime=$(date +%s)
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) | tee -a "${baseDir}"/sectiontime.txt


#####################
#                   #
#    Text report    #
#                   #
#####################


cat "${baseDir}"/sectiontime.txt >  "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt

cat "${baseDir}"/section5.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt

echo "These files did not get renamed:" >> "${baseDir}"/log.txt
[ -f "${baseDir}"/csection1.txt ] && cat "${baseDir}"/csection1.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt

echo "Possible mixed isolates (defining SNPs called AC=1)" >> "${baseDir}"/log.txt
cat "${baseDir}"/section2.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt
cat "${baseDir}"/section3.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt

echo "SNP counts:" >> "${baseDir}"/log.txt
cat "${baseDir}"/section4.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************" >> "${baseDir}"/log.txt

echo "AC1 called SNPs" >> "${baseDir}"/log.txt
cat "${baseDir}"/all_vcfs/all_vcfs_AC1Report.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************" >> "${baseDir}"/log.txt


#####################
#                   #
#    HTML report    #
#                   #
#####################


echo "<html>" > "${baseDir}"/email_log.html
echo "<Body>" >> "${baseDir}"/email_log.html

cat "${baseDir}"/sectiontime.txt \
    | awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' \
    >  "${baseDir}"/email_log.html
echo -e "\n****************************************************\n" >> "${baseDir}"/email_log.html

cat "${baseDir}"/section5.txt \
    | awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' \
    >> "${baseDir}"/email_log.html
echo -e "\n****************************************************\n" >> "${baseDir}"/email_log.html

echo "<p> These files did not get renamed: </p>" >> "${baseDir}"/email_log.html
[ -f "${baseDir}"/csection1.txt ] && cat "${baseDir}"/csection1.txt \
    | awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" "$i""</td>";print "</tr>"} END{print "</table>"}' \
    >> "${baseDir}"/email_log.html
echo -e "\n****************************************************\n" >> "${baseDir}"/email_log.html

echo "<p> Possible Mixed Isolates, Defining SNPs called AC=1 </p>" >> "${baseDir}"/email_log.html
cat "${baseDir}"/section2.txt \
    | awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" "$i""</td>";print "</tr>"} END{print "</table>"}' \
    >> "${baseDir}"/email_log.html
echo -e "\n****************************************************\n" >> "${baseDir}"/email_log.html

cat "${baseDir}"/section3.txt \
    | awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>"$i"</td>";print "</tr>"} END{print "</table>"}' \
    >> "${baseDir}"/email_log.html
echo -e "\n****************************************************\n" >> "${baseDir}"/email_log.html

echo "<p> SNP counts: </p>" >> "${baseDir}"/email_log.html
cat "${baseDir}"/section4.txt \
    | awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' \
    >> "${baseDir}"/email_log.html
echo -e "\n****************************************************\n" >> "${baseDir}"/email_log.html

echo "<p> AC1 called SNPs: </p>" >> "${baseDir}"/email_log.html
cat "${baseDir}"/all_vcfs/all_vcfs_AC1Report.txt \
    | awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' \
    >> "${baseDir}"/email_log.html

echo "</Body>" >> "${baseDir}"/email_log.html
echo "</html>" >> "${baseDir}"/email_log.html

#Cleanup
[ -f "${baseDir}"/section1.txt ] && rm "${baseDir}"/section1.txt
[ -f "${baseDir}"/csection1.txt ] && rm "${baseDir}"/csection1.txt
rm "${baseDir}"/section2.txt
rm "${baseDir}"/section3.txt
rm "${baseDir}"/section4.txt
rm "${baseDir}"/section5.txt
rm "${baseDir}"/sectiontime.txt
rm -r "$filterdir"


#####################
#                   #
#    Send report    #
#                   #
#####################


if [ "$mflag" ]; then
    email_list="marc-olivier.duceppe@inspection.gc.ca"
    mail -s "Resuts from script 2" -A "${baseDir}"/email_log.html -t "$email_list" <<< "${baseDir}"/email_log.html
else
    mail -s "Resuts from script 2" -A "${baseDir}"/email_log.html -t "$email_list" <<< "${baseDir}"/email_log.html
fi

#Cleanup
# rm "${baseDir}"/email_log.html

#Closing comment
echo -e "\n****************************** END ******************************\n"
