#!/bin/bash

: <<'END'
This script is the second script in a two script workflow. 
Script 2 genotypes Mycobacterium tuberculosis complex and Brucella species from SNP data contained in VCFs.
It operates on VCFs generated with the same reference output from script 1.
VCFs are collected into a single working directory.
Comparisons are output as SNP tables and alignment files to view as trees in your program of choice.

Script 2 will run and output tables and alignments from just the data generated from script 1, 
however the data will be more informative if additional files are provide.  Those files are:
1) A file that contains positions to cluster individual isolates/VCFs into groups, subgroups and clades.
2) Files that contain positions to remove from the analysis.

Paradigm
1) Once a SNP occurs and establishes in a population it does not revert back
2) Observations of homoplasy are rare
3) Group, subgroup and clade clusters only show parsimony informative SNPs for the isolates within that cluster
4) SNPs observed in a single isolate are less informative than SNPs seen in multiple isolates 
and therefore established in a population

END


######################
#                    #
#    User Defined    #
#                    #
######################


#Where analysis will take place
baseDir=""${HOME}"/analyses/mbovis_script2"

#Where the VCF files are
vcfPath="/home/bioinfo/Desktop/vcf"

#Coverage below this value will be changed to -
Ncov="1" 


#####################
#                   #
#   Data Stucture   #
#                   #
#####################


#script dependencies
dependents=""${HOME}"/prog/snp_analysis/script_dependents"

#time stamp of analysis
uniqdate="$(date "+%Y-%m-%dat%Hh%Mm%Ss")"
# uniqdate="2016-08-24at17h08m21s"
filterdir=""${baseDir}"/"${uniqdate}"-FilterFiles"

#Files containing positions to filter
FilterDirectory="$filterdir"


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


echo -e "\n****************************** START ******************************\n"

echo "Start Time: "$(date)"" > "${baseDir}"/sectiontime.txt
starttime=$(date +%s)
argUsed="$1"
dircalled="$vcfPath"
echo "start time: "$uniqdate""


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



function filterFileCreations ()
{
    # Use to make filter files from the text pasted from the Excel worksheet.
    # working directory does not need to be set.
    #   Set variables:

    # Path to txt file containing paste from Excel worksheet.
    filterFile=""${filterdir}"/filterFile.txt"

    # Number of columns in Excel worksheet
    columns=$(cat "$filterFile" | head | awk 'BEGIN{ FS="\t"; OFS="\t" }  END {print NF}')

    # Location filter files are output to.
    output="$filterdir"

    let columns=columns+1

    echo "Filter sets: "$columns""
    echo "Extracting from Excel to text files..."

    count=1
    while [ "$count" -lt "$columns" ]; do
        #echo "$count"
        filename=$(cat "$filterFile" | awk -v x="$count" 'BEGIN{FS=OFS="\t"}{print $x}' | head -n 1)
        #echo "Filename: "$filename""
        cat "$filterFile" | awk -v x="$count" 'BEGIN{FS=OFS="\t"} FNR>1 {print $x}' | grep -v "^$" > "${output}"/"${filename}".list
        let count=count+1
    done

    #cleanup
    rm "$filterFile"

    #expand the ranges
    for i in $(find "$output" -type f | grep -F ".list"); do
        # i="/home/bioinfo/analyses/mbvis_script2/2016-08-18at16h26m36s-FilterFiles/Subgroup-23D.list.txt"
        base=$(basename "$i")
        readyfile=$(echo "$base" | sed 's/\..*//')

        mylist=($(cat "$i")) #array
        # echo "${mylist[@]}"

        for l in "${mylist[@]}"; do
            pos1=$(echo  "$l" | sed 's/-/ /g' | awk '{print $1}')
            pos2=$(echo  "$l" | sed 's/-/ /g' | awk '{print $2}')
            #echo $pos2 #
            if [ -z "$pos2" ]; then
            let pos2=pos1+1
                while [ "$pos1" -lt "$pos2" ]; do
                    #echo $pos1 #
                    echo "$pos1" >> "${output}"/"${readyfile}".txt
                    let pos1=pos1+1
                done
            else
                let pos2=pos2+1
                while [ "$pos1" -lt "$pos2" ]; do
                    #echo $pos1 #
                    echo "$pos1" >> "${output}"/"${readyfile}".txt
                    let pos1=pos1+1
                done
            fi
        done
    done

    rm "${output}"/*.list
}

function parseXLS ()
{
#install python module without su rights
# mkdir -p "${HOME}"/local/lib/python2.7/site-packages
# easy_install --install-dir=/home/CFIA-ACIA/duceppem/local/lib/python2.7/site-packages/ xlrd

    cat >"${baseDir}"/inputXLS.py <<'EOL'
#!/usr/bin/env python

import os
import xlrd
from sys import argv

script, input = argv

wb = xlrd.open_workbook(input)
wb.sheet_names()
#sh = wb.sheet_by_index(1)
sh = wb.sheet_by_name(u'New groupings')
for rownum in range(sh.nrows):
    print sh.row_values(rownum)

EOL

    chmod 755 "${baseDir}"/inputXLS.py
    python "${baseDir}"/inputXLS.py "$excelinfile"
    rm "${baseDir}"/inputXLS.py

}



# Environment controls:

if [ "$1" == "ab1" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt" #to change the name of vcf files
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_abortus/Abortus1_Defining_SNPs.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_abortus/FilterFiles" #Files containing positions to filter
    
    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using Brucella abortus bv 1, 2 or 4 variables" | tee "${baseDir}"/section5.txt

    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"

elif [ "$1" == "mel" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_melitensis/Mel_Defining_SNPs.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_melitensis/FilterFiles" #Files containing positions to filter
    
    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. melitensis variables" | tee "${baseDir}"/section5.txt

    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"

elif [ "$1" == "suis1" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_suis_bv1/Suis1_Defining_SNPs.txt"

    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_suis_bv1/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. suis bv1 variables" | tee "${baseDir}"/section5.txt

    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"

elif [ "$1" == "suis2" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_suis_bv2/suis2_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_suis_bv2/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. suis bv2 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"

elif [ "$1" == "suis3" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_suis_bv3/Suis3_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_suis_bv3/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. suis bv3 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"

elif [ "$1" == "suis4" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_suis_bv4/Suis4_Defining_SNPs.txt"

    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_suis_bv/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. suis bv4 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"

elif [ "$1" == "canis" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_canis/Canis_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_canis/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "vcftofasta.sh ran as B. canis"
    echo "Script vcftofasta.sh ran using B. canis variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"


elif [ "$1" == "ceti1" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_ceti-grp1/Ceti1_Defining_SNPs.txt"

    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_ceti-grp1/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B ceti group 1 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"


elif [ "$1" == "ceti2" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_ceti-grp2/Ceti2_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_ceti-grp2/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B ceti group 2 variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"


elif [ "$1" == "ovis" ]; then

    genotypingcodes=""${dependents}"/genotyping_codes/bruc_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Brucella_ovis/Ovis_Defining_SNPs.txt"

    FilterAllVCFs="no" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades
    FilterDirectory=""${dependents}"/Brucella_ovis/FilterFiles" #Files containing positions to filter

    QUAL=300 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=350 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using B. ovis variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"

elif [ "$1" == "bovis" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/Mycobacterium_bovis/DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using M. bovis variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    excelinfile=""${dependents}"/Mycobacterium_bovis/Filtered_Regions.xlsx"

    if [ "$eflag" ]; then
        echo "Only the "elite" bovis isolates are being ran"
    else
        echo "All bovis are being ran"
        echo "Like to run selected isolates? Use... vcftofasta.sh -e bovis"
    fi

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script2/Filtered_Regions.xlsx
    # Excel tab label "New groupings"

    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    wait

    filterFileCreations
    wait

elif [ "$1" == "tb1" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB1/tb1DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB1/tb1Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    wait
    filterFileCreations
    wait

elif [ "$1" == "tb2" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB2/tb2DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB2/tb2Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    wait
    filterFileCreations
    wait

elif [ "$1" == "tb3" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB3/tb3DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="no" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB3/tb3Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    filterFileCreations

elif [ "$1" == "tb4a" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB4a/tb4aDefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB4a/tb4aFiltered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    wait
    filterFileCreations
    wait

elif [ "$1" == "tb4b" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB4b/tb4bDefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB4b/tb4bFiltered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    wait
    filterFileCreations
    wait

elif [ "$1" == "tb5" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB5/tb5DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB5/tb5Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    wait
    filterFileCreations
    wait

elif [ "$1" == "tb6" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/TB6/tb6DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using ${1} variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"
    excelinfile=""${dependents}"/TB6/tb6Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    wait
    filterFileCreations
    wait

elif [ "$1" == "para" ]; then
    genotypingcodes=""${dependents}"/genotyping_codes/tb_tags.txt"
    # This file tells the script how to cluster VCFs
    DefiningSNPs=""${dependents}"/paraTB/DefiningSNPsGroupDesignations.txt"
    FilterAllVCFs="yes" #(yes or no), Do you want to filter all VCFs?
    FilterGroups="yes" #(yes or no), Do you want to filter VCFs withing their groups, subgroups, and clades

    QUAL=150 # Minimum quality for calling a SNP
    lowEnd=1
    highEnd=200 # QUAL range to change ALT to N

    echo "Script vcftofasta.sh ran using para variables" | tee "${baseDir}"/section5.txt
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaia@inspection.gc.ca"

    # For tb inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    # Excel file that is being used is at: /bioinfo11/TStuber/Results/mycobacterium/vcfs/Filtered_Regions.xlsx
    # Excel tab label "New groupings"

    excelinfile=""${dependents}"/paraTB/Filtered_Regions.xlsx"
    parseXLS | sed 's/ u//g' | tr "," "\t" \
        | sed -e 's/\[//g' -e 's/\]//g' -e 's/ //g' -e 's/^u//g' -e 's/\.0//g' \
        | tr -d "'"  \
        > "${filterdir}"/filterFile.txt
    wait
    filterFileCreations
    wait

else

    echo ""
    echo "Incorrect argument!  Must use one of the following arguments: ab1, mel, suis1, suis2, suis3, suis4, canis, ceti1, ceti2, ovis, bovis, tb1, tb2, tb3, tb4a, tb4b, tb5, tb6, para"
    echo ""
    echo "Set optional flags"
    echo "flag -c with look for positions to filter.  By default, with no -c, this will not be done."
    echo "flag -m will email just "M"e"
    echo "flag -e will run the bovis "E"lite representative samples"
    echo "flag -a get "a"ll_vcf alignment table"
    echo ""
    echo "Example: [prompt]$ vcftofasta.sh -mea bovis"
    echo ""
    rm "${baseDir}"/sectiontime.txt
    exit 1

fi


# If there are 2 vcf files with the same name one of the files might unknowingly
# get cut out of the analysis and keep the undesired vcf instead.  This will
# alert if 2 vcf with the same TB number are present.
# The regular expression used in sed should be changed based on vcf naming convention

function testDuplicates ()
{
    echo "Checking for empty or duplicated VCF files."

    directorytest="${baseDir##*/}" #name of directory where script was launched (or name of "$baseDir")
    if [ "$directorytest" == "VCF_Source_All" ]; then # if where all the reference vcf files are stored
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

# Looks for defining positions in VCF files.
# If an AC=1 is found at a defined position it is flagged as a posible mixed infection.
# These defining positions must be SNPs found cluster's main branch

function AConeCallPosition ()
{
    positionList=$(cat "$DefiningSNPs" | awk ' { print $2 }' | awk ' NF > 0 ')

    echo "AConeCallPosition is running, started -->  "$(date)""
    #echo "*********************************************************************" >> "${baseDir}"/section2.txt
    #echo "Possible Mixed Isolates" > "${baseDir}"/section2.txt
    #echo "Defining SNPs that are called as AC=1" >> "${baseDir}"/section2.text
    [ -e "${baseDir}"/section2.txt ] && rm "${baseDir}"/section2.txt ]
    for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
        for pos in "$positionList"; do
            awk -v x="$pos" 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ "^"x"$" ) print FILENAME, "Pos:", $2, "QUAL:", $6, $8 }' "$i"
        done \
            | grep "AC=1;A" \
            | awk 'BEGIN {FS=";"} {print $1, $2}' \
            >> "${baseDir}"/section2.txt
    done
}

#################################################################################

# This function prepares the filter files.
# awk needs to see a number in the file, so if the file is blank 2 matching numbers are added.  2 numbers because duplicates are kept therefore allowing a number to be pasting into awk when comparing files.

function filterFilespreparation ()
{
    # For tb, inputXLS.py creates text files with positions to be filetered, and places them in FilterDirectory
    echo "Waiting for filter file creation to complete"

    echo "Preparing Filter Files"
    for i in $(find -L "$FilterDirectory" -type f | grep -F ".txt"); do
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
    echo "Changeing low calls, started --> "$(date)""
    for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
        base=$(basename "$i")
        baseNoExt="${base%.*}"

        cat "$i" | awk -v x="$lowEnd" -v y="$highEnd" 'BEGIN {OFS="\t"} { if ($6 >= x && $6 <= y) print $1, $2, $3, $4, "N", $6, $7, $8; else print $0 }' \
            > "${baseDir}"/"${baseNoExt}".txt

        rm "$i"
        mv "${baseDir}"/"${baseNoExt}".txt "${baseDir}"/"${baseNoExt}".vcf
    done
}

#################################################################################

function findpositionstofilter ()
{
    echo "Finding positions to filter --> "$(date)""

    d="$1"

    # positions have already been filtered via cutting specific positions.
    cp "${d}"/filtered_total_pos.txt "${d}"/total_pos.txt
    cat "${d}"/total_pos.txt | awk '{print $1}' > "${d}"/prepositionlist.txt

    for n in $(cat "${d}"/prepositionlist.txt); do
        front=$(echo "$n" | sed 's/\(.*\)-\([0-9]*\)/\1/')
        back=$(echo "$n" | sed 's/\(.*\)-\([0-9]*\)/\2/')
        # echo "front: $front"
        # echo "back: $back"

        positioncount=$(cat $(find "${d}" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" ' $1 == f && $2 == b {count++} END {print count}')

        #echo "position count: $positioncount"
        if [ "$positioncount" -gt 2 ]; then
            #printf "%s\t%s\n" "$front" "$back"
            echo "$n" >> "${d}"/positionlist.txt
        else
            echo "$n" >> "${d}"-DONOT_filtertheseposition.txt # "$d" from fasta_table ()
        fi
    done

    echo "Filtering --> "$(date)""

    for p in $(cat "${d}"/positionlist.txt); do
        front=$(echo "$p" | sed 's/\(.*\)-\([0-9]*\)/\1/')
        back=$(echo "$p" | sed 's/\(.*\)-\([0-9]*\)/\2/')
        #echo "front: $front"
        #echo "back: $back"

        maxqual=$(cat $(find "$d" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" 'BEGIN{max=0} $1 == f && $2 == b {if ($6>max) max=$6} END {print max}' \
            | sed 's/\..*//')

        avequal=$(cat $(find "${d}" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" '$6 != "." && $1 == f && $2 == b {print $6}' \
            | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' \
            | sed 's/\..*//')

        maxmap=$(cat $(find "$d" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" ' $1 == f && $2 == b {print $8}' \
            | sed 's/.*MQ=\(.....\).*/\1/' | awk 'BEGIN{max=0}{if ($1>max) max=$1} END {print max}' \
            | sed 's/\..*//')

        avemap=$(cat $(find "$d" -maxdepth 1 -type f | grep -F ".vcf") \
            | awk -v f="$front" -v b="$back" '$6 != "." && $1 == f && $2 == b {print $8}' \
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
    # fasta_table "${baseDir}"/all_groups

    # Loop through the directories
    directories=($(find "$1" -maxdepth 1 -type d | awk 'NR > 1'))
    # echo "${directories[@]}" | tr " " "\n"

    for d in "${directories[@]}"; do
        # d=""${baseDir}"/all_groups/Group-9"
        dName=$(basename "$d")

        #backup vcf files
        [ -d "${d}"/starting_files ] || mkdir "${d}"/starting_files
        cp "${d}"/*.vcf "${d}"/starting_files

        if [ "$FilterGroups" == "yes" ]; then
            if [ "$chromCount" -eq 1 ]; then
                #Mark vcf allowing areas of the genome to be removed from the SNP analysis
                for i in $(find "$d" -maxdepth 1 -type f | grep -F ".vcf"); do
                    # i=""${d}"/MBWGS083.SNPsZeroCoverage.vcf"
                    m=$(basename "$i")
                    n=$(echo "$m" | sed "$dropEXT")

                    #SNP without missing genotypes (./.)
                    cat "$i" \
                        | awk '$1 !~ /#/ && $10 !~ /\.\/\./ {print $2}' \
                        > "${i}".file

                    cat "${FilterDirectory}"/"${dName}".txt "${i}".file \
                        >> "${i}".catFile

                    cat "${i}".catFile \
                        | sort | uniq -d \
                        > "${i}".txt

                    #create the regex expressions to filter out the SNPs based on their positions
                    pos=$(cat "${i}".txt \
                        | tr "\n" "W" \
                        | sed -e 's/W/\$\|\^/g' -e 's/\$\|\^$//' -e 's/$/\$/' -e 's/^/\^/' -e 's/|$$//')

                    #If "$pos is empty, assign a default value larger than the chromosome"
                    [ -z "$pos" ] && pos="^1000000000$"

                    #remove all "$pos" from "$i"
                    cat "$i" \
                        | awk -v x="$pos" 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' \
                        > "${d}"/"${n}".filtered.vcf

                    rm "${i}".file
                    rm "${i}".catFile
                    rm "${i}".txt
                    cat "${d}"/"${n}".filtered.vcf \
                        | grep -v "Not_Included" \
                        > "$i"
                done

            else
                #Mark vcf allowing areas of the genome to be removed from the SNP analysis
                for i in $(find "$d" -maxdepth 1 -type f | grep -F ".vcf"); do
                    m=$(basename "$i")
                    n=$(echo "$m" | sed "$dropEXT") # n is name with all right of "_" and "." removed.

                    cat "$i" | grep '^#' > "${i}".header
                    cat "$i" | grep -v '^#' > "${i}".body

                    #Mark vcf allowing areas of the genome to be removed from the SNP analysis
                    # Iterate through chrom number range
                    COUNTER=0
                    for c in $(cat "${baseDir}"/chroms.txt); do
                        let COUNTER=COUNTER+1

                        cat "${i}".body \
                            | awk -v c="$c" 'BEGIN{OFS="\t"} $1 !~ /#/ && $10 !~ /\.\/\./ && $1 == c {print $2}' \
                            > "${i}".filepositions

                        cat "${FilterDirectory}"/"${d}".txt \
                            | awk -v c="$c" ' $1 == c {print $2}' \
                            > "${i}".positionstofilter

                        cat "${i}".positionstofilter "${i}".filepositions \
                            | sort -k1,1 \
                            | uniq -d \
                            > "${i}".foundpositions

                        pos=$(cat "${i}".foundpositions \
                            | tr "\n" "W" \
                            | sed -e 's/W/\$\|\^/g' -e 's/\$\|\^$//' -e 's/$/\$/' -e 's/^/\^/' -e 's/|$$//')

                        if [ -n "$pos" ]; then
                            echo "pos: "$pos"" > /dev/null 2>&1
                        else
                        #    echo "string is zero; no findings for pos; giving pos=1"
                            pos="^1$"
                        fi

                        cat "${i}".body \
                            | awk -v var1="$c" -v var2="$pos" 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' \
                            | grep "$c" \
                            > "${d}"/"${n}".filterchrom"${COUNTER}".vcf
                    done

                    cat "${i}".header "${d}"/"${n}".filterchrom*.vcf > "${d}"/"${n}".filtered.vcf
                    cat "${d}"/"${n}".filtered.vcf \
                        | grep -v "Not_Included"  \
                        > "$i"

                    rm "${i}".header
                    rm "${i}".body
                    rm "${i}".filepositions
                    rm "${i}".positionstofilter
                    rm "${d}"/"${n}".filterchrom*.vcf
                    rm "${i}".foundpositions

                done
            fi
            rm "${d}"/*.filtered.vcf
        fi

        # Make concatemer with the position and REF call.
        # Factor in possible multiple chromosomes
        # Get rid of duplicates in concatemer and "${baseDir}"/list.txt all the positions and REF calls
        for i in $(find "$d" -maxdepth 1 -type f | grep -F ".vcf"); do
            cat "$i" | awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $4}' >> "${d}"/concatemer.txt
        done

        # Get rid of duplicates in concatemer and "${baseDir}"/list.txt all the positions and REF calls
        cat "${d}"/concatemer.txt \
            | sort -k1,1 \
            | uniq \
            > "${d}"/filtered_total_alt.txt

        cat "${d}"/filtered_total_alt.txt \
            | awk '{print $1}'  \
            > "${d}"/filtered_total_pos.txt

        # Count the number of SNPs
        totalSNPs=$(cat "${d}"/filtered_total_pos.txt | wc -l)
        echo "Total SNPs: "$totalSNPs""

        ######################## FILTER FILE CREATOR ###########################

        if [ "$cflag" ]; then
            findpositionstofilter "$d"
            wait
        fi

        #########################################################################

        # Find AC1 positions also found in "${d}"/total_pos.txt
        cat "${d}"/filtered_total_pos.txt \
            | awk '{print $1}' \
            > "${d}"/total.list.txt

        for i in $(find "$d" -maxdepth 1 -type f | grep -F ".vcf"); do
            # i=""${d}"/MBWGS083.SNPsZeroCoverage.vcf"
            m=$(basename "$i")
            n=$(echo "$m" | sed 's/\..*//')

            # search for AC1 positions
            cat "$i" \
                | awk ' $0 !~ /^#/ && $8 ~ /^AC=1/ && $6 > 0 {print $1 "-" $2}' \
                > "${d}"/"${n}".list.txt

            # AC1 positions that are being found in this group
            positionsfound=$(cat "${d}"/"${n}".list.txt "${d}"/total.list.txt | sort -n | uniq -d)
            countfind=$(echo "$positionsfound" | wc -w)

            #echo "positonsfound: $positionsfound  countfind: $countfind"
            rm "${d}"/"${n}".list.txt

            [ -z "$positionsfound" ] && positionsfound="No positions found"
            
            if [ "$countfind" -gt 2  ]; then
                searchname=$(echo "$n" | sed 's/_.*//')

                # if [ "$argUsed" == "para" ]; then
                #     unmappedContigs=$(grep -A 1 "Unmapped contig count" /bioinfo11/TStuber/Results/mycobacterium/mac/para_cattle-bison/data/"${searchname}"/bwamem-gatk/qualityvalues/*stats.txt)
                # elif [  "$argUsed" == "bovis" ]; then
                #     unmappedContigs=$(grep -A 1 "Unmapped contig count" /bioinfo11/TStuber/Results/mycobacterium/tbc/tbbov/script1/"${searchname}"/bwamem-gatk/qualityvalues/*stats.txt)
                # else
                # contigMessage="possibly set a new contig path at script line: $LINENO"
                # # fi

                # if [ -z "$unmappedContigs" ]; then
                #     unmappedContigs="Contig counts not available"
                # fi

                echo -e ""$dName"\n" > "${d}"/"${dName}"-AC1postions.txt

                echo  ""$dName" Sample: "$n"  AC1 findings: "$countfind"  "$unmappedContigs" "$contigMessage"" > "${d}"/delete.txt

                cat "${d}"/delete.txt \
                    | tr -d "\n" \
                    >> "${d}"/"${dName}"-AC1postions.txt

                echo "" >> "${d}"/"${dName}"-AC1postions.txt

                cat "${d}"/delete.txt \
                    | tr -d "\n" \
                    >> "${baseDir}"/emailAC1counts.txt

                echo "" >> "${baseDir}"/emailAC1counts.txt
                
                for p in $(echo "$positionsfound"); do
                    position=$(echo "$p" | sed 's/chrom[0-9]*-//')

                    cat "$i" \
                        | awk -v p="$position" '$2 == p {print $0}' \
                        >> "${d}"/"${dName}"-AC1postions.txt
                done

                rm "${d}"/delete.txt
            fi
        done
        
        rm "${d}"/total.list.txt

        echo "***Creating normalized vcf using AC2, QUAL > "$QUAL""
        # Grab the name of the vcf file

        #########################################################################

        # Count the number of SNPs
        filteredSNPs=$(cat "${d}"/filtered_total_pos.txt | wc -l)
        echo ""$dName" total SNPs after filtering: "$filteredSNPs"" | tee -a "${baseDir}"/section4.txt

        for i in $(find "$d" -maxdepth 1 -type f | grep -F ".vcf"); do
            # i=""${d}"/MBWGS083.SNPsZeroCoverage.vcf"
            n="${i%.vcf}"

            cat "$i" \
                | awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $5}' \
                > "${n}".allsnps_alt

            #get SNPs of interest
            cat "${n}".allsnps_alt \
                | fgrep -f "${d}"/filtered_total_pos.txt \
                > "${n}".targetsnps_alt

            #if SNP not found in sample default call to reference, normalize.
            cat "${n}".targetsnps_alt "${d}"/filtered_total_alt.txt  \
                | awk '{ if (a[$1]++ == 0) print $0; }' \
                | sort -nk1,1 \
                > "${n}".filteredsnps_alt

            # If position has zero map quality change alt call to -
            # get positions being used
            cat "${n}".filteredsnps_alt \
                | awk '{print $1}' \
                > "${n}".filteredsnps_pos

            # Get zero coverage positions.
            cat "$i" \
                | awk ' $0 !~ /^#/ && $10 ~ /\.\/\./ {print $1 "-" $2}' \
                > "${n}".zeropositions

            # if duplicate then zero mapped position found for sample
            cat "${n}".filteredsnps_pos "${n}".zeropositions \
                | sort \
                | uniq -d \
                | awk '{print $1, "-"}' \
                > "${n}".zerotomerge_alt #the - makes it and alt file

            #if zero positions found merge them to the SNPs found
            if [ -s "${n}".zerotomerge_alt ]; then
                # merge zero updates to SNP file
                cat "${n}".zerotomerge_alt "${n}".filteredsnps_alt \
                    | awk '{ if (a[$1]++ == 0) print $0; }' \
                    | sort -nk1,1 \
                    > "${n}".zerofilteredsnps_alt

                #echo "***Found zero postions: $n"
                rm "${n}".filteredsnps_alt
            else
                #echo "no zero positions found for $n"
                mv "${n}".filteredsnps_alt "${n}".zerofilteredsnps_alt
            fi

            rm "${n}".allsnps_alt
            rm "${n}".filteredsnps_pos
            rm "${n}".targetsnps_alt
            rm "${n}".zeropositions
            rm "${n}".zerotomerge_alt
        done

        echo "Finding parsimony informative positions --> "$(date)""

        ### Not useful if only one sample in group

        # Capture only positions that have more than one SNP type called at a position

        cat "${d}"/*zerofilteredsnps_alt \
            | sort -nk1,1 \
            | uniq \
            | awk '{print $1}' \
            | uniq -d \
            > "${d}"/parsimony_informative.txt #empty if only one sample in the group

        # This removes calls that are the same for all isolates being analyzed

        # If many SNPs fgrep may not do much and be slow
        cat "${d}"/filtered_total_alt.txt \
            | fgrep -f "${d}"/parsimony_informative.txt \
            | sort -k1,1n \
            > "${d}"/parsimony_filtered_total_alt.txt

        cat  "${d}"/parsimony_filtered_total_alt.txt \
            | awk '{print $1}' \
            > "${d}"/parsimony_filtered_total_pos.txt

        # Create table and fasta
        cat "${d}"/parsimony_filtered_total_alt.txt \
            | awk '{print $1}' \
            | awk 'BEGIN{print "reference_pos"}1' \
            | tr '\n' '\t' \
            | sed 's/$//' \
            | awk '{print $0}' \
            >> "${d}"/"${dName}".table.txt

        cat "${d}"/parsimony_filtered_total_alt.txt \
            | awk '{print $2}' \
            | awk 'BEGIN{print "reference_call"}1' \
            | tr '\n' '\t' \
            | sed 's/$//' \
            | awk '{print $0}' \
            >> "${d}"/"${dName}".table.txt


        #Make the fasta files
        for i in "${d}"/*zerofilteredsnps_alt; do
            m=$(basename "$i")
            n=$(echo "$m" | sed 's/\..*//')

            cat "$i" \
                | fgrep -f "${d}"/parsimony_filtered_total_pos.txt  \
                | sort -k1,1n \
                > "${d}"/"${n}".pretod

            ##############################################################
            # Change AC1s to IUPAC

            # get positions being used
            cat "${d}"/"${n}".pretod \
                | awk '{print $1}' \
                > "${d}"/"${n}".usedpostions

            # get AC1 positions and iupac calls  that were changed to iupac
            cat ${i%zerofilteredsnps_alt}vcf \
                | awk -v Q="$QUAL" ' $0 !~ /#/ && $6 > Q && $8 ~ /^AC=1;/ {print $1 "-" $2, $5}' \
                > "${d}"/"${n}".ac

            # get just positions of those AC1 grabbed above
            cat "${d}"/"${n}".ac \
                | awk '{print $1}' \
                > "${d}"/"${n}".acpositions

            # AC duplicate positions will need to be kept
            cat "${d}"/"${n}".usedpostions "${d}"/"${n}".acpositions \
                | sort | uniq -d \
                > "${d}"/"${n}".actokeep

            # get AC1 position with iupac, these are only positions already in the pretod

            if [ -s "${d}"/"${n}".actokeep ]; then
                cat "${d}"/"${n}".ac \
                    | fgrep -f "${d}"/"${n}".actokeep \
                    > "${d}"/"${n}".actomerge

                # merge iupac updates to filledcut
                cat "${d}"/"${n}".actomerge "${d}"/"${n}".pretod \
                    | awk '{ if (a[$1]++ == 0) print $0; }' \
                    | sort -nk1,1 \
                    > "${d}"/"${n}".tod

                rm "${d}"/"${n}".pretod
                rm "${d}"/"${n}".actomerge
            else
                #echo "else done"
                mv "${d}"/"${n}".pretod "${d}"/"${n}".tod
            fi
            
            #Cleanup
            rm "${d}"/"${n}".usedpostions
            rm "${d}"/"${n}".ac
            rm "${d}"/"${n}".acpositions
            rm "${d}"/"${n}".actokeep
            ##############################################################

            cat "${d}"/"${n}".tod \
                | awk '{print $2}' \
                | tr -d [:space:] \
                | sed "s/^/>$n;/" \
                | tr ";" "\n" \
                | sed 's/[A-Z],[A-Z]/N/g' \
                > "${d}"/"${n}".fas

            echo "" >> "${d}"/"${n}".fas

            # Add each isolate to the table
            cat "${d}"/"${n}".tod \
                | awk '{print $2}' \
                | awk -v number="$n" 'BEGIN{print number}1' \
                | tr '\n' '\t' \
                | sed 's/$//' \
                | awk '{print $0}' \
                >> "${d}"/"${dName}".table.txt 
        done

        #Create root sequence
        cat "${d}"/parsimony_filtered_total_alt.txt \
            | awk '{print $2}' \
            > "${d}"/root.txt

        cat "${d}"/root.txt \
            | tr -cd "[:print:]" \
            | sed "s/^/>root;/" \
            | tr ";" "\n" \
            | sed 's/[A-Z],[A-Z]/N/g' \
            > "${d}"/root.fas

        echo "" >> "${d}"/root.fas

        totalSNPs=$(cat "${d}"/parsimony_filtered_total_pos.txt | grep -c ".*")
        echo "Total informative SNPs: "$totalSNPs""

        # Make a file containing all fasta files. Used awk instead of cat to insure newline between files
        cat "${d}"/*.fas \
            | awk '{print $0}' \
            > "${d}"/"${dName}"_alignment.fasta

        #Clean-up
        rm "${d}"/concatemer.txt
        rm "${d}"/*.tod
        mkdir "${d}"/fasta
        mv "${d}"/*.fas "${d}"/fasta
        #rm "${d}"/root.txt
        rm "${d}"/*vcf
        rm "${d}"/filtered_total_alt.txt
        rm "${d}"/filtered_total_pos.txt
        rm "${d}"/parsimony_filtered_total_alt.txt
        rm "${d}"/parsimony_filtered_total_pos.txt
        rm "${d}"/parsimony_informative.txt
        rm "${d}"/*zerofilteredsnps_alt
    done
}

#****************************************************************

function addQuality ()
{

    # Add map qualities to sorted table

    # Get just the position.  The chromosome must be removed
    cat "$1" \
        | awk ' NR == 1 {print $0}' \
        | tr "\t" "\n" \
        | sed "1d" \
        | awk '{print NR, $0}' \
        > "${d}"/"${dName}"-positions.txt

    echo "map-quality map-quality" > "${d}"/quality.txt
    echo "Sorted table map quality gathering for "$dName" --> "$(date)""

    #This loop takes long to run
    while IFS= read -r line || [[ -n "$line" ]]; do
        rownumber=$(echo "$line" | awk '{print $1}')
        front=$(echo "$line" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\1/')
        back=$(echo "$line" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\2/')

        [ -z "$front" ] && continue #skip the last line where there is only a line number and no SNP info

        #echo "rownumber: $rownumber"
        #echo "front: $front"
        #echo "back: $back"
        avemap=$(cat "${parent}"/starting_files/*vcf \
            | awk -v f="$front" -v b="$back" '$6 != "." && $1 == f && $2 == b {print $8}' \
            | sed 's/.*MQ=\(.....\).*/\1/' \
            | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' \
            | sed 's/\..*//')

        echo ""$rownumber" "$avemap"" >> "${d}"/quality.txt

    done < "${d}"/"${dName}"-positions.txt
    
    wait

    cat "${d}"/quality.txt \
        | sort -nk1,1 \
        | awk '{print $2}' \
        | tr "\n" "\t" \
        > "${d}"/qualitytransposed.txt

    cat "$1" "${d}"/qualitytransposed.txt \
        | grep -v '^$' \
        > "${d}"/"${dName}"-mapquality-sortedtable.txt

    #crush original file
    mv "${d}"/"${dName}"-mapquality-sortedtable.txt "$1"

}



function alignTable ()
{

    d="$1"
    # d=""${baseDir}"/all_groups/Group-18/fasta"
    # d=""${baseDir}"/all_vcfs/fasta"

    #Group/Subgroup/Clade name
    parent=$(dirname "$d")
    dName="$(basename "$parent")"
    
    # Beginning in fasta folder
    echo "RAxML started on "$dName" --> "$(date)""


    cat "${d}"/*.fas \
        | awk '{print $0}'  \
        | sed '/root/{N;d;}' \
        >> "${d}"/fastaGroup.txt

    cat "${d}"/*.fas \
        | awk '{print $0}' \
        >> "${d}"/RAxMLfastaGroup.txt

    #Make the tree
    raxmlHPC-AVX \
        -T "$cpu" \
        -s "${d}"/RAxMLfastaGroup.txt \
        -w "$d" \
        -n "$dName" \
        -m GTRCAT \
        -p 12345 \
        &>/dev/null
    wait

    #reroot tree
    nw_reroot \
        "${d}"/RAxML_bestTree."$dName" \
        root \
        | tee >(tee "${d}"/tableinput."$dName" "${d}"/rooted_RAxML_bestTree."$dName" 2&> /dev/null) \
        | tee >(nw_display \
                -s \
                -w 1000 \
                -v 20 \
                -b 'opacity:0' \
                -i 'font-size:8' \
                -l 'font-family:serif;font-style:italic' \
                -d 'stroke-width:2;stroke:blue' \
                /dev/stdin \
                    | tee >(tee "${d}"/"${dName}"-tree.svg 2&> /dev/null) \
                    | tee >(rsvg-convert \
                            -f pdf \
                            "${d}"/"${dName}"-tree.svg \
                            > "${d}"/"${dName}"-tree.pdf))
    wait

    mv "${d}"/rooted_RAxML_bestTree."$dName" "${d}"/RAxML_bestTree."$dName"

    #cleanup
    rm "${d}"/RAxML_parsimonyTree*

    cat "${d}"/tableinput."$dName" \
        | tr ":" "\n" \
        | tr "," "\n" \
        | sed -e 's/(//g' -e 's/)//g' \
        | grep -v "\.[0-9]*" \
        | grep -v "root" \
        > "${d}"/cleanedAlignment.txt

    awk 'NR==FNR{o[FNR]=$1; next} {t[$1]=$0} END{for(x=1; x<=FNR; x++){y=o[x]; print t[y]}}' "${d}"/cleanedAlignment.txt "${parent}"/"${dName}".table.txt \
        > "${d}"/joined.txt
    
    cat "${parent}"/"${dName}".table.txt \
        | grep "reference" \
        > "${d}"/references.txt

    cat "${d}"/references.txt "${d}"/joined.txt >> "${d}"/joined2.txt

    mv "${d}"/joined2.txt "${d}"/"${dName}".sortedTable.txt

    rm "${d}"/joined.txt
    rm "${d}"/references.txt

    #get the number of columns
    columnCount=$(cat "${d}"/"${dName}".sortedTable.txt | awk '$0 !~ /^$/ {print $0}' | awk '{print NF-1; exit}')

    #number=`jot - 1 $columnCount`
    number=($(seq "$columnCount"))
    #echo "Numbers in "${baseDir}"/list.txt: $number"

    #get the reference row
    cat "${d}"/"${dName}".sortedTable.txt \
        | sed -n 2p \
        | awk 'BEGIN{FS=OFS="\t"} {$1="";sub("\t","")}1' \
        | tee "${d}"/referenceRow.txt "${d}"/out2.txt > /dev/null

    #remove first column (Sample name)
    cat "${d}"/"${dName}".sortedTable.txt \
        | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' \
        > "${d}"/table.txt

    echo "countDif" > "${d}"/countOutput.txt
    echo "countFirst" > "${d}"/firstOutput.txt

    #iterate numbers up to number of columns
    for n in "${number[@]}"; do
        #use number to grab character in reference (first) row i.e. C
        letter=$(cat "${d}"/referenceRow.txt | awk -v x="$n" '{print $x}')

        #SNP sequence for all samples
        cat "${d}"/table.txt \
            | awk -v var1="$n" '{print $var1}' \
            > "${d}"/column.txt

        #number of samples with the reference SNP
        cat "${d}"/column.txt \
            | grep "$letter" \
            | wc -l \
            >> "${d}"/countOutput.txt

        #remove 2 fisrt lines of SNP sequence. Why???
        cat "${d}"/column.txt \
            | sed '1,2d' \
            > "${d}"/column2.txt

        #Get the position of the fisrt occurence of the alternative allele
        cat "${d}"/column2.txt \
            | awk -v var2="$letter" ' $0 !~ var2 {print NR; exit}' \
            >> "${d}"/firstOutput.txt

        # cat "${d}"/column.txt \
        #     | awk -v var2="$letter" ' $0 !~ var2 {print NR; exit}' \
        #     >> "${d}"/firstOutput.txt
    done

    #copy content of sorted table
    cat "${d}"/"${dName}".sortedTable.txt > "${d}"/table2.txt

    #Prepare count line
    cat "${d}"/countOutput.txt \
        | sed 's/ //g' \
        | tr "\n" "\t" \
        > "${d}"/readyline.txt

    #Create readytable.txt with count line.
    cat "${d}"/table2.txt "${d}"/readyline.txt \
        > "${d}"/orginizedTable2.txt

    #Remove empty lines
    cat "${d}"/orginizedTable2.txt \
        | grep -vE "^$" \
        > "${d}"/orginizedTable3.txt

    #Prepare firstOut line
    cat "${d}"/firstOutput.txt \
        | sed 's/ //g' \
        | tr "\n" "\t" \
        > "${d}"/readyFirstOut.txt

    cat "${d}"/orginizedTable3.txt "${d}"/readyFirstOut.txt \
        > "${d}"/orginizedTable4.txt

    #Cleanup
    rm "${d}"/referenceRow.txt
    rm "${d}"/out2.txt
    rm "${d}"/table.txt
    rm "${d}"/countOutput.txt
    rm "${d}"/table2.txt
    rm "${d}"/readyline.txt
    rm "${d}"/orginizedTable2.txt

    #move last line fist
    cat "${d}"/orginizedTable4.txt \
        | awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' \
        > "${d}"/orginizedTable5.txt

    #move last line fist
    cat "${d}"/orginizedTable5.txt \
        | awk '{a[NR]=$0} END {print a[NR]; for (i=1;i<NR;i++) print a[i]}' \
         > "${d}"/orginizedTable6.txt

    #Transpose
    cat "${d}"/orginizedTable6.txt \
        | awk '
            {
              for(c = 1; c <= NF; c++) {
                a[c, NR] = $c
              }
              if(max_nf < NF) {
                max_nf = NF
              }
            }
            END {
              for(r = 1; r <= max_nf; r++) {
                for(c = 1; c <= NR; c++) {
                  printf("%s ", a[r, c])
                }
                print ""
              }
            }
            ' > "${d}"/orginizedTable7.txt

    wait

    #Orgainize file based on 1st 2 columns
    cat "${d}"/orginizedTable7.txt \
        | sort -n -k1 \
        | sort -n -k2 \
        > "${d}"/orginizedTable8.txt

    #Convert spaces to tabs
    cat "${d}"/orginizedTable8.txt \
        | awk -v OFS="\t" '$1=$1' \
        > "${d}"/orginizedTable9.txt

    cat "${d}"/orginizedTable9.txt \
        | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' \
        | awk 'BEGIN{FS=OFS="\t"}{$1="";sub("\t","")}1' \
        > "${d}"/orginizedTable10.txt

    #Transpose back
    cat "${d}"/orginizedTable10.txt \
        | awk '
            {
              for(c = 1; c <= NF; c++) {
                a[c, NR] = $c
              }
              if(max_nf < NF) {
                max_nf = NF
              }
            }
            END {
              for(r = 1; r <= max_nf; r++) {
                for(c = 1; c <= NR; c++) {
                  printf("%s ", a[r, c])
                }
                print ""
              }
            }
            ' >  "${d}"/orginizedTable11.txt

    wait

    #Convert spaces to tabs
    cat "${d}"/orginizedTable11.txt \
        | awk -v OFS="\t" '$1=$1' \
        > "${d}"/"${dName}".organizedTable.txt


    rm "${d}"/orginizedTable3.txt
    rm "${d}"/orginizedTable4.txt
    rm "${d}"/orginizedTable5.txt
    rm "${d}"/orginizedTable6.txt
    rm "${d}"/orginizedTable7.txt
    rm "${d}"/orginizedTable8.txt
    rm "${d}"/orginizedTable9.txt
    rm "${d}"/orginizedTable10.txt
    rm "${d}"/orginizedTable11.txt
    rm "${d}"/column.txt
    rm "${d}"/column2.txt
    rm "${d}"/readyFirstOut.txt
    rm "${d}"/firstOutput.txt


    echo "Adding map qualities..."


    ######################
    #                    #
    #    Sorted table    #
    #                    #
    ######################


    
    # Add map qualities to sorted table
    addQuality "${d}"/"${dName}".sortedTable.txt
    addQuality "${d}"/"${dName}".organizedTable.txt

    # # Get just the position.  The chromosome must be removed
    # cat "${d}"/"${dName}".sortedTable.txt \
    #     | awk ' NR == 1 {print $0}' \
    #     | tr "\t" "\n" \
    #     | sed "1d" \
    #     | awk '{print NR, $0}' \
    #     > "${d}"/"${dName}"-positions.txt

    # echo "map-quality map-quality" > "${d}"/quality.txt
    # echo "Sorted table map quality gathering for "$dName" --> "$(date)""

    # #This loop takes long to run
    # while IFS= read -r line || [[ -n "$line" ]]; do
    #     rownumber=$(echo "$line" | awk '{print $1}')
    #     front=$(echo "$line" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\1/')
    #     back=$(echo "$line" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\2/')

    #     [ -z "$front" ] && continue #skip the last line where there is only a line number and no SNP info

    #     #echo "rownumber: $rownumber"
    #     #echo "front: $front"
    #     #echo "back: $back"
    #     avemap=$(cat "${parent}"/starting_files/*vcf \
    #         | awk -v f="$front" -v b="$back" '$6 != "." && $1 == f && $2 == b {print $8}' \
    #         | sed 's/.*MQ=\(.....\).*/\1/' \
    #         | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' \
    #         | sed 's/\..*//')

    #     echo ""$rownumber" "$avemap"" >> "${d}"/quality.txt

    # done < "${d}"/"${dName}"-positions.txt
    
    # wait

    # cat "${d}"/quality.txt \
    #     | sort -nk1,1 \
    #     | awk '{print $2}' \
    #     | tr "\n" "\t" \
    #     > "${d}"/qualitytransposed.txt

    # cat "${d}"/"${dName}".sortedTable.txt "${d}"/qualitytransposed.txt \
    #     | grep -v '^$' \
    #     > "${d}"/"${dName}"-mapquality-sortedtable.txt

    # #crush original file
    # mv "${d}"/"${dName}"-mapquality-sortedtable.txt "${d}"/"${dName}".sortedTable.txt


    #########################
    #                       #
    #    Organized table    #
    #                       #
    #########################


    # # Add map qualities to organized table
    # cat "${d}"/"${dName}".organizedTable.txt \
    #     | awk ' NR == 1 {print $0}' \
    #     | tr "\t" "\n" \
    #     | sed "1d" \
    #     | awk '{print NR, $0}' \
    #     > "${d}"/"${dName}"-positions.txt

    # echo "map-quality map-quality" > "${d}"/quality.txt
    # echo "Organized table map quality gathering for "$dName" --> "$(date)""

    # #This loop takes long to run
    # while IFS= read -r line || [[ -n "$line" ]]; do
    #     rownumber=$(echo "$line" | awk '{print $1}')
    #     front=$(echo "$line" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\1/')
    #     back=$(echo "$line" | awk '{print $2}' | sed 's/\(.*\)-\([0-9]*\)/\2/')

    #     [ -z "$front" ] && continue #skip the last line where there is only a line number and no SNP info

    #     #echo "rownumber: $rownumber"
    #     #echo "front: $front"
    #     #echo "back: $back"
    #     avemap=$(cat "${parent}"/starting_files/*vcf \
    #         | awk -v f="$front" -v b="$back" '$6 != "." && $1 == f && $2 == b {print $8}' \
    #         | sed 's/.*MQ=\(.....\).*/\1/' \
    #         | awk '{ sum += $1; n++ } END { if (n > 0) print sum / n; }' \
    #         | sed 's/\..*//')

    #     echo ""$rownumber" "$avemap"" >> "${d}"/quality.txt
    # done < "${d}"/"${dName}"-positions.txt

    # wait

    # cat "${d}"/quality.txt \
    #     | sort -nk1,1  \
    #     | awk '{print $2}' \
    #     | tr "\n" "\t" \
    #     > "${d}"/qualitytransposed.txt

    # cat "${d}"/"${dName}".organizedTable.txt "${d}"/qualitytransposed.txt \
    #     | grep -v '^$' \
    #     > "${d}"/"${dName}"-mapquality-orgainizedtable.txt

    # mv "${d}"/"${dName}"-mapquality-orgainizedtable.txt "${d}"/"${dName}".organizedTable.txt

    #Cleanup
    rm "${d}"/quality.txt
    rm "${d}"/qualitytransposed.txt
    rm "${d}"/"${dName}"-positions.txt
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
echo "Counting the number of chromosomes in first 100 samples, started -->  "$(date)""
chromCount=$(cat $(find -L "$baseDir" -type f | grep -F ".vcf" | head -100) \
                | awk ' $0 !~ /^#/ {print $1}' \
                | sort | uniq -d \
                | awk 'END {print NR}')

echo "The number of chromosomes/segments seen in VCF files: "$chromCount""
cat $(find -L "$baseDir" -type f | grep -F ".vcf" | head -100) \
    | awk ' $0 !~ /^#/ {print $1}' \
    | sort | uniq -d \
    > "${baseDir}"/chroms.txt

echo "These are the chromosomes/segments found:"
cat "${baseDir}"/chroms.txt


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


########################################################################

printf "%s\t%s\t%s\t%s\n" "TB Number" "Group" "Subgroup" "Clade" \
    > "${baseDir}"/FileMakerGroupImport.txt

AConeCallPosition
wait

# Change low QUAL SNPs to N, see set variables
#gets rid of symbolic links and saves "real" files in "$baseDir"
changeLowCalls
wait

######################## Change AC1s to IUPAC ########################

echo "Changing AC=1 to IUPAC, started -->  "$(date)""

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

######################## Mark Files and Remove Marked Regions ########################
echo "chromCount:  "$chromCount""

if [ "$FilterAllVCFs" == yes ]; then
    echo "$(date) --> Marking all VCFs and removing filtering region"
    # Label filter field for positions to be filtered in all VCFs
    if [ "$chromCount" -eq 1 ]; then
        for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
            m=$(basename "$i")
            n=$(echo "$m" | sed "$dropEXT")

            # Get usable positions in the VCF
            cat "$i" | awk '$1 !~ /#/ && $10 !~ /\.\/\./ {print $2}' > "${i}".file

            # Combine with positions that will be filtered
            cat "${FilterDirectory}/FilterToAll.txt" "${i}".file >> "${i}".catFile
            
            # Output any duplicate positions, aka decreasing positions to be marked and used by awk
            cat "${i}".catFile | sort | uniq -d > "${i}".txt

            # preparing postions
            pos=$(cat "${i}".txt \
                | tr "\n" "W" \
                | sed -e 's/W/\$\|\^/g' -e 's/\$\|\^$//' -e 's/$/\$/' -e 's/^/\^/' -e 's/|$$//')
    
            # If no positions found to be filtered a filler is needed or all positions will get marked as "Not_Included"
            if [ -z "$pos" ]; then
                pos="100000000"
            fi

            # Making a vcf with positions marked that should not be included based on filter file
            cat "$i" \
                | awk -v x="$pos" 'BEGIN {FS="\t"; OFS="\t"} { if($2 ~ x ) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' \
                > "${baseDir}"/"${n}".filtered.vcf
            
            # Clean up
            rm "${i}".file
            rm "${i}".catFile
            rm "${i}".txt

            # Removed positions and clobber origingal vcf
            cat "${baseDir}"/"${n}".filtered.vcf \
                | grep -v "Not_Included"  \
                > "$i"
        done
        wait

    elif [ "$chromCount" -gt 1 ]; then
        #echo "multiple chromosomes"
        for i in *.vcf; do m=$(basename "$i"); n=$(echo "$m" | sed "$dropEXT") # n is name with all right of "_" and "." removed.
            cat "$i" | grep '^#' > "${i}".header
            cat "$i" | grep -v '^#' > "${i}".body

            #Mark vcf allowing areas of the genome to be removed from the SNP analysis
            # Iterate through chrom number range
            COUNTER=0
            for c in $(cat "${baseDir}"/chroms.txt); do
                let COUNTER=COUNTER+1
                #echo The counter is $COUNTER
                #echo "********* In all_vcfs $n working on chromos $c **********"
                cat "${i}".body \
                    | awk -v c="$c" 'BEGIN{OFS="\t"} $1 !~ /#/ && $10 !~ /\.\/\./ && $1 == c {print $2}' \
                    > "${i}".filepositions

                cat "${FilterDirectory}"/FilterToAll.txt \
                    | awk -v c="$c" ' $1 == c {print $2}' \
                    > "${i}".positionstofilter

                cat "${i}".positionstofilter "${i}".filepositions \
                    | sort -k1,1 \
                    | uniq -d \
                    > "${i}".foundpositions

                pos=$(cat "${i}".foundpositions \
                    | tr "\n" "W" \
                    | sed -e 's/W/\$\|\^/g' -e 's/\$\|\^$//' -e 's/$/\$/' -e 's/^/\^/' -e 's/|$$//')

                if [ -n "$pos" ]; then
                    echo "pos: "$pos"" > /dev/null 2>&1
                else
                    #echo "string is zero; no findings for pos; giving pos=1"
                    pos="^1$"
                    #echo $pos
                fi

                cat "${i}".body \
                    | awk -v var1=$c -v var2=$pos 'BEGIN {FS="\t"; OFS="\t"} { if($1 ~ var1 && $2 ~ var2) print $1, $2, $3, $4, $5, $6, "Not_Included", $8, $9, $10; else print $0}' \
                    | grep "$c" \
                    > "${n}".filterchrom"${COUNTER}".vcf
            done

            #Create ouptut VCF file
            cat "${i}".header "${n}".filterchrom*.vcf \
                > "${n}".filtered.vcf

            #Replace input file with output
            cat "${n}".filtered.vcf \
                | grep -v "Not_Included" \
                > "$i"

            #Cleanup
            rm "${i}".header
            rm "${i}".body
            rm "${i}".filepositions
            rm "${i}".positionstofilter
            rm "${n}".filterchrom*.vcf
            rm "${i}".foundpositions
        done
    else
        echo "Check chromosome count numbers at line "$LINENO".  Exiting script."
        exit 1
    fi

    wait

    rm "${baseDir}"/*.filtered.vcf 

else
    echo "***All VCF filtering was NOT done." | tee -a "${baseDir}"/section5.txt
fi

wait

#################### Categorize VCFs into Groups, Subgroups and Clades #####################

# Print header
echo -e "NAME\tGROUP\tSUBGROUP\tCLADE" > "${baseDir}"/section3.txt

for i in $(find -L "$baseDir" -type f | grep -F ".vcf"); do
    # If there is one chromosome present just get the position.
    # If multiple chromosomes are present than the chromsome identification needs to be identified.
    # The filter file needs to sync with this chromosome identification.
    # If multiple chromosomes the filter file will be kept as a text file.
    # If a single chromosome an linked Excel file can be used.

    if [ "$chromCount" -eq 1 ]; then
        # Get quality positions in VCF
        cat "$i" | awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/{print $2}' > "${i%.vcf}.quality"
    else
        # Get quality positions in VCF and include chromosome identification
        cat "$i" | awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2}' > "${i%.vcf}.quality"
    fi

    echo ""${i%.vcf}.quality":"

    ##----Group

    # If a group number matches a quality position in the VCF (formatedpos) then print the position
    cat "$DefiningSNPs" | grep "Group" > "${baseDir}"/groupsnps.txt

    awk 'NR==FNR{a[$0];next}$2 in a' "${i%.vcf}.quality" "${baseDir}"/groupsnps.txt | awk '{print $1}' > "${i%.vcf}.group-foundpositions"

    echo "This is the Group Numbers: $(cat "${i%.vcf}.group-foundpositions")"

    # Typically a single group position is found, and the VCF will be placed into just one group.  It is posible that an isolate will need to go in more than one group because of were it falls on the tree.  In this case there may be 2 group, or more, group positions found.  The number of group positions found is captured in sizeGroup.
    sizeGroup=$(cat "${i%.vcf}.group-foundpositions" | wc -l)

    # Loop through the number of groups positions found
    loops=$(cat "${i%.vcf}.group-foundpositions")

    if [ "$sizeGroup" -lt 1 ]; then # There was not a position found that places VCF into group
        echo ""$i" Group not found" >> "${baseDir}"/section3.txt
        echo ""$i" was not assigned a Group"
    elif [ "$sizeGroup" -gt 1 ]; then
        echo ""$i" has multiple groups" >> "${baseDir}"/section3.txt
        echo ""$i" has multiple groups"
        for l in "$loops"; do
            echo "Making group "$i""
            mkdir -p "${baseDir}"/Group-"${l}" #Make groupNumber folder if one does not exist.
            cp "$i" "${baseDir}"/Group-"${l}"/ #Then copy to each folder
        done
    else
        echo "Just one group"
        mkdir -p "${baseDir}"/Group-"${loops}" #Make groupNumber folder if one does not exist.
        cp "$i" "${baseDir}"/Group-"${loops}"/ #Then copy to each folder
    fi

    ##----Subgroup

    # If a group number matches a quality position in the VCF (formatedpos) then print the position
    cat "$DefiningSNPs" | grep "Subgroup" > "${baseDir}"/subgroupsnps.txt

    awk 'NR==FNR{a[$0];next}$2 in a' "${i%.vcf}.quality" "${baseDir}"/subgroupsnps.txt | awk '{print $1}' > "${i%.vcf}.subgroup-foundpositions"

    echo "This is the Subgroup Numbers: $(cat "${i%.vcf}.subgroup-foundpositions")"

    # Typically a single group position is found, and the VCF will be placed into just one group.  It is posible that an isolate will need to go in more than one group because of were it falls on the tree.  In this case there may be 2 group, or more, group positions found.  The number of group positions found is captured in sizeGroup.
    sizeGroup=$(cat "${i%.vcf}.subgroup-foundpositions" | wc -l)

    # Loop through the number of groups positions found
    loops=$(cat "${i%.vcf}.subgroup-foundpositions")

    if [ "$sizeGroup" -lt 1 ]; then # There was not a position found that places VCF into group
        echo ""$i" was not assigned a Subgroup"
    elif [ "$sizeGroup" -gt 1 ]; then
        echo ""$i" has multiple subgroups" >> "${baseDir}"/section3.txt
        echo ""$i" has multiple subgroups"
        for l in "$loops"; do
            echo "Making subgroup "$i""
            mkdir -p "${baseDir}"/Subgroup-"${l}" #Make groupNumber folder if one does not exist.
            cp "$i" "${baseDir}"/Subgroup-"${l}"/ #Then copy to each folder
        done
    else
        echo "Just one Subgroup"
        mkdir -p "${baseDir}"/Subgroup-"${loops}" #Make groupNumber folder if one does not exist.
        cp "$i" "${baseDir}"/Subgroup-"${loops}"/ #Then copy to each folder
    fi

    ##----Clade

    # If a group number matches a quality position in the VCF (formatedpos) then print the position
    cat "$DefiningSNPs" | grep "Clade" > "${baseDir}"/cladesnps.txt

    awk 'NR==FNR{a[$0];next}$2 in a' "${i%.vcf}.quality" "${baseDir}"/cladesnps.txt | awk '{print $1}' > "${i%.vcf}.clade-foundpositions"

    echo "This is the Clade Numbers: $(cat "${i%.vcf}.clade-foundpositions")"

    # Typically a single group position is found, and the VCF will be placed into just one group.
    #It is posible that an isolate will need to go in more than one group because of were it falls on the tree.
    #In this case there may be 2 group, or more, group positions found.  The number of group positions found is captured in sizeGroup.
    sizeGroup=$(cat "${i%.vcf}.clade-foundpositions" | wc -l)

    # Loop through the number of groups positions found
    loops=$(cat "${i%.vcf}.clade-foundpositions")

    if [ "$sizeGroup" -lt 1 ]; then # There was not a position found that places VCF into group
        echo ""$i" was not assigned a Clade"
    elif [ "$sizeGroup" -gt 1 ]; then
        echo ""$i" has multiple clades" >> "${baseDir}"/section3.txt
        echo ""$i" has multiple clades"
        for l in "$loops"; do
            echo "making clade "$i""
            mkdir -p "${baseDir}"/Clade-"${l}" #Make groupNumber folder if one does not exist.
            cp "$i" "${baseDir}"/Clade-"${l}"/ #Then copy to each folder
        done
    else
        echo "Just one clade"
        mkdir -p "${baseDir}"/Clade-"${loops}" #Make groupNumber folder if one does not exist.
        cp "$i" "${baseDir}"/Clade-"${loops}"/ #Then copy to each folder
    fi

    echo "${i%.vcf} $(cat "${i%.vcf}.group-foundpositions" "${i%.vcf}.subgroup-foundpositions" "${i%.vcf}.clade-foundpositions")" \
        | tr "\n" "\t" >> "${baseDir}"/section3.txt
    echo "" >> "${baseDir}"/section3.txt

    #cleanup
    rm "${i%.vcf}.quality"
    rm "${baseDir}"/*snps.txt
    rm "${baseDir}"/*foundpositions

    ######

    [ -d "${baseDir}"/all_vcfs ] || mkdir -p "${baseDir}"/all_vcfs #Make all_vcfs folder if one does not exist.
    mv "$i" "${baseDir}"/all_vcfs/
done


################### Organize folders #####################


[ -d "${baseDir}"/all_groups ] || mkdir -p "${baseDir}"/all_groups
mv "${baseDir}"/Group-* "${baseDir}"/all_groups

[ -d "${baseDir}"/all_subgroups ] || mkdir -p "${baseDir}"/all_subgroups
mv "${baseDir}"/Subgroup*/ "${baseDir}"/all_subgroups/

[ -d "${baseDir}"/all_clades ] || mkdir -p "${baseDir}"/all_clades
mv "${baseDir}"/Clade*/ "${baseDir}"/all_clades/


##################### Start: All vcf folder #####################


[ -d "${baseDir}"/all_vcfs/starting_files ] || mkdir "${baseDir}"/all_vcfs/starting_files
cp "${baseDir}"/all_vcfs/*vcf "${baseDir}"/all_vcfs/starting_files

# Make concatemer with the position and REF call.
# Factor in possible multiple chromosomes
# Get rid of duplicates in concatemer and "${baseDir}"/list.txt all the positions and REF calls
echo "Gathering SNP positions --> "$(date)""


for i in $(find "${baseDir}"/all_vcfs -maxdepth 1 -type f | grep -F ".vcf"); do
    cat "$i" | awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $4}' >> "${baseDir}"/all_vcfs/concatemer.txt
done

# Get rid of duplicates in concatemer and "${baseDir}"/list.txt all the positions and REF calls
cat "${baseDir}"/all_vcfs/concatemer.txt | sort -k1,1 | uniq > "${baseDir}"/all_vcfs/filtered_total_alt.txt
cat "${baseDir}"/all_vcfs/filtered_total_alt.txt | awk '{print $1}' > "${baseDir}"/all_vcfs/filtered_total_pos.txt

# Count the number of SNPs
totalSNPs=$(cat "${baseDir}"/all_vcfs/filtered_total_pos.txt | wc -l)
echo "Total filtered SNPs: "$totalSNPs""


for i in $(find "${baseDir}"/all_vcfs -maxdepth 1 -type f | grep -F ".vcf"); do
    n="${i%.vcf}"

    cat "$i" | awk -v Q="$QUAL" ' $0 !~ /^#/ && $6 > Q && $8 ~ /^AC=2;/ {print $1 "-" $2, $5}' > "${n}".allsnps_alt

    #get SNPs of interest
    cat "${n}".allsnps_alt | fgrep -f "${baseDir}"/all_vcfs/filtered_total_pos.txt > "${n}".targetsnps_alt

    #if SNP not found in sample default call to reference, normalize.
    cat "${n}".targetsnps_alt "${baseDir}"/all_vcfs/filtered_total_alt.txt \
        | awk '{ if (a[$1]++ == 0) print $0; }' \
        | sort -nk1,1 \
        > "${n}".filteredsnps_alt

    # If position has zero map quality change alt call to -
    # get positions being used
    cat "${n}".filteredsnps_alt \
        | awk '{print $1}' \
        > "${n}".filteredsnps_pos

    # Get zero coverage positions.
    cat "$i" \
        | awk ' $0 !~ /^#/ && $10 ~ /\.\/\./ {print $1 "-" $2}' \
        > "${n}".zeropositions

    # if duplicate then zero mapped position found for sample
    cat "${n}".filteredsnps_pos "${n}".zeropositions \
        | sort | uniq -d \
        | awk '{print $1, "-"}' \
        > "${n}".zerotomerge_alt #the - makes it and alt file

    #if zero positions found merge them to the SNPs found
    if [ -s "${n}".zerotomerge_alt ]; then
        # merge zero updates to SNP file
        cat "${n}".zerotomerge_alt "${n}".filteredsnps_alt \
            | awk '{ if (a[$1]++ == 0) print $0; }' \
            | sort -nk1,1 \
            > "${n}".zerofilteredsnps_alt

        #echo "***Found zero postions: $n"
        rm "${n}".filteredsnps_alt

    else
        #echo "no zero positions found for $n"
        mv "${n}".filteredsnps_alt "${n}".zerofilteredsnps_alt
    fi

    rm "${n}".allsnps_alt
    rm "${n}".filteredsnps_pos
    rm "${n}".targetsnps_alt
    rm "${n}".zeropositions
    rm "${n}".zerotomerge_alt
done

wait



echo "Finding parsimony informative positions --> "$(date)""

# Capture only positions that have more than one SNP type called at a position
cat "${baseDir}"/all_vcfs/*zerofilteredsnps_alt \
    | sort -nk1,1 \
    | uniq \
    | awk '{print $1}' \
    | uniq -d \
    > "${baseDir}"/all_vcfs/parsimony_informative.txt

# This removes calls that are the same for all isolates being analyzed
# If many SNPs fgrep may not do much and be slow
cat "${baseDir}"/all_vcfs/filtered_total_alt.txt \
    | fgrep -f "${baseDir}"/all_vcfs/parsimony_informative.txt \
    | sort -k1,1n \
    > "${baseDir}"/all_vcfs/parsimony_filtered_total_alt.txt

cat "${baseDir}"/all_vcfs/parsimony_filtered_total_alt.txt \
    | awk '{print $1}' \
    > "${baseDir}"/all_vcfs/parsimony_filtered_total_pos.txt


######################## FILTER FILE CREATOR ###########################


if [ "$cflag" ]; then
    findpositionstofilter "${baseDir}"/all_vcfs
    wait
fi


#########################################################################


# Create table and fasta
cat "${baseDir}"/all_vcfs/parsimony_filtered_total_alt.txt \
    | awk '{print $1}' \
    | awk 'BEGIN{print "reference_pos"}1' \
    | tr '\n' '\t' \
    | sed 's/$//' \
    | awk '{print $0}' \
    > "${baseDir}"/all_vcfs/all_vcfs.table.txt

cat "${baseDir}"/all_vcfs/parsimony_filtered_total_alt.txt \
    | awk '{print $2}' \
    | awk 'BEGIN{print "reference_call"}1' \
    | tr '\n' '\t' \
    | sed 's/$//' \
    | awk '{print $0}' \
    >> "${baseDir}"/all_vcfs/all_vcfs.table.txt


for i in "${baseDir}"/all_vcfs/*zerofilteredsnps_alt; do
    m=$(basename "$i")
    n=$(echo "$m" | sed 's/\..*//')
    # n="${i%.vcf}"

    cat "$i" | fgrep -f "${baseDir}"/all_vcfs/parsimony_filtered_total_pos.txt \
        | sort -k1,1n \
        > "${baseDir}"/all_vcfs/"${n}".pretod

    ##############################################################
    # Change AC1s to IUPAC

    # get positions being used
    cat "${baseDir}"/all_vcfs/"${n}".pretod \
        | awk '{print $1}' \
        > "${baseDir}"/all_vcfs/"${n}".usedpostions

    # get AC1 positions and iupac calls  that were changed to iupac
    cat "${i%zerofilteredsnps_alt}"vcf \
        | awk -v Q="$QUAL" ' $0 !~ /#/ && $6 > Q && $8 ~ /^AC=1;/ {print $1 "-" $2, $5}' \
        > "${baseDir}"/all_vcfs/"${n}".ac

    # get just positions of those AC1 grabbed above
    cat "${baseDir}"/all_vcfs/"${n}".ac \
        | awk '{print $1}' \
        > "${baseDir}"/all_vcfs/"${n}".acpositions

    # AC duplicate positions will need to be kept
    cat "${baseDir}"/all_vcfs/"${n}".usedpostions "${baseDir}"/all_vcfs/"${n}".acpositions \
        | sort | uniq -d \
        > "${baseDir}"/all_vcfs/"${n}".actokeep

    # get AC1 position with iupac, these are only positions already in the pretod

    if [ -s "${baseDir}"/all_vcfs/"${n}".actokeep ]; then
        cat "${baseDir}"/all_vcfs/"${n}".ac \
            | fgrep -f "${baseDir}"/all_vcfs/"${n}".actokeep \
            > "${baseDir}"/all_vcfs/"${n}".actomerge

        # merge iupac updates to filledcut
        cat "${baseDir}"/all_vcfs/"${n}".actomerge "${baseDir}"/all_vcfs/"${n}".pretod \
            | awk '{ if (a[$1]++ == 0) print $0; }' \
            | sort -nk1,1 \
            > "${baseDir}"/all_vcfs/"${n}".tod

        rm "${baseDir}"/all_vcfs/"${n}".pretod
        rm "${baseDir}"/all_vcfs/"${n}".actomerge
    else
        #echo "else done"
        mv "${baseDir}"/all_vcfs/"${n}".pretod "${baseDir}"/all_vcfs/"${n}".tod
    fi

    rm "${baseDir}"/all_vcfs/"${n}".usedpostions
    rm "${baseDir}"/all_vcfs/"${n}".ac
    rm "${baseDir}"/all_vcfs/"${n}".acpositions
    rm "${baseDir}"/all_vcfs/"${n}".actokeep

    ##############################################################

    cat "${baseDir}"/all_vcfs/"${n}".tod \
        | awk '{print $2}' \
        | tr -d [:space:] \
        | sed "s/^/>$n;/" \
        | tr ";" "\n" \
        | sed 's/[A-Z],[A-Z]/N/g' \
        > "${baseDir}"/all_vcfs/"${n}".fas

    echo "" >> "${baseDir}"/all_vcfs/"${n}".fas

    # Add each isolate to the table
    cat "${baseDir}"/all_vcfs/"${n}".tod \
        | awk '{print $2}' \
        | awk -v number="$n" 'BEGIN{print number}1' \
        | tr '\n' '\t' \
        | sed 's/$//' \
        | awk '{print $0}' \
        >> "${baseDir}"/all_vcfs/all_vcfs.table.txt
done


#Create root sequence
cat "${baseDir}"/all_vcfs/parsimony_filtered_total_alt.txt \
    | awk '{print $2}' \
    > "${baseDir}"/all_vcfs/root.txt

cat "${baseDir}"/all_vcfs/root.txt \
    | tr -cd "[:print:]" \
    | sed "s/^/>root;/" \
    | tr ";" "\n" \
    | sed 's/[A-Z],[A-Z]/N/g' \
    > "${baseDir}"/all_vcfs/root.fas

echo "" >> "${baseDir}"/all_vcfs/root.fas

totalSNPs=$(cat "${baseDir}"/all_vcfs/parsimony_filtered_total_pos.txt | grep -c ".*")
echo "Total informative SNPs: "$totalSNPs""


#Clean-up
rm "${baseDir}"/all_vcfs/concatemer.txt
rm "${baseDir}"/all_vcfs/*.tod

[ -d "${baseDir}"/all_vcfs/fasta ] || mkdir "${baseDir}"/all_vcfs/fasta
mv "${baseDir}"/all_vcfs/*.fas "${baseDir}"/all_vcfs/fasta

rm "${baseDir}"/all_vcfs/root.txt
rm "${baseDir}"/all_vcfs/*vcf
rm "${baseDir}"/all_vcfs/filtered_total_alt.txt
rm "${baseDir}"/all_vcfs/filtered_total_pos.txt
rm "${baseDir}"/all_vcfs/parsimony_filtered_total_alt.txt
rm "${baseDir}"/all_vcfs/parsimony_filtered_total_pos.txt
rm "${baseDir}"/all_vcfs/parsimony_informative.txt
rm "${baseDir}"/all_vcfs/*zerofilteredsnps_alt


if [ "$eflag" -o "$aflag" ]; then
    # d="all_vcfs"
    # cd ./fasta
    alignTable "${baseDir}"/all_vcfs/fasta
else
    echo "Tree not ran for all_vcfs"
fi

##################### End: All vcf folder #####################


######################
#                    #
#    Fasta Tables    #
#                    #
######################


#all_groups
if [ -n "$(ls -A "${baseDir}"/all_groups)" ]; then 
    fasta_table "${baseDir}"/all_groups
    wait
else
    echo "No Groups"
fi

#all_subgroups
if [ -n "$(ls -A "${baseDir}"/all_subgroups)" ]; then 
    fasta_table "${baseDir}"/all_subgroups
    wait
else
    echo "No SubGroup"
fi

#all_clades
if [ -n "$(ls -A "${baseDir}"/all_clades)" ]; then 
    fasta_table "${baseDir}"/all_clades
    wait
else
    echo "No Clades"
fi


#For QA
# cp "$DefiningSNPs" "$baseDir"
# cp "${dependents}"/Table_Template.xlsx "$baseDir"
# cp "$0" "$baseDir" #$0 is the name of the script itself


#########################
#                       #
#    Table alignment    #
#                       #
#########################


directories=($(find "$baseDir" -maxdepth 1 -type d \
    | awk 'NR > 1' \
    | grep -v "$uniqdate" \
    | grep -vF "all_vcfs"))

# echo "${directories[@]}" | tr " " "\n"

#all_groups/all_subgroups/all_clades
for d in "${directories[@]}"; do
    # echo ""$d""
    if [ -n "$(ls -A "$d")" ]; then
        subdirectories=($(find "$d" -maxdepth 1 -type d | awk 'NR > 1'))

        for sd in "${subdirectories[@]}"; do
            # echo -e "\t"$sd""
            alignTable "${sd}"/fasta
            wait
        done

    else
        echo "Nothing to do"
    fi
done


#################################################################


cat "${baseDir}"/section1.txt \
    column \
    > "${baseDir}"/csection1.txt

cat "${baseDir}"/section4.txt \
    | sort -nr \
    > "${baseDir}"/ssection4.txt

#Elapsed time
echo "End Time: "$(date)"" >> "${baseDir}"/sectiontime.txt
endtime=$(date +%s)
runtime=$((endtime-starttime))
printf 'Runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) >> "${baseDir}"/sectiontime.txt


cat "${baseDir}"/sectiontime.txt >  "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt
cat "${baseDir}"/section5.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt
echo "These files did not get renamed:" >> "${baseDir}"/log.txt
cat c"${baseDir}"/section1.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt
echo "Possible Mixed Isolates" >> "${baseDir}"/log.txt
echo "Defining SNPs called AC=1" >> "${baseDir}"/log.txt
cat "${baseDir}"/section2.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************\n" >> "${baseDir}"/log.txt
cat "${baseDir}"/section3.txt >> "${baseDir}"/log.txt
echo "" >> "${baseDir}"/log.txt
echo "****************************************************" >> "${baseDir}"/log.txt
echo "SNP counts::" >> "${baseDir}"/log.txt
cat s"${baseDir}"/section4.txt >> "${baseDir}"/log.txt
echo -e "\n****************************************************" >> "${baseDir}"/log.txt
echo "AC1 called SNPs"
cat "${baseDir}"/emailAC1counts.txt | sort -nk1,1 >> "${baseDir}"/log.txt

echo "<html>" > "${baseDir}"/email_log.html
echo "<Body>" >> "${baseDir}"/email_log.html
cat "${baseDir}"/sectiontime.txt | awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' >  "${baseDir}"/email_log.html
echo -e "****************************************************\n" >> "${baseDir}"/email_log.html
cat "${baseDir}"/section5.txt | awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' >> "${baseDir}"/email_log.html
echo "" >> "${baseDir}"/email_log.html
echo -e "****************************************************\n" >> "${baseDir}"/email_log.html
echo "<p> These files did not get renamed: </p>" >> "${baseDir}"/email_log.html
cat c"${baseDir}"/section1.txt | awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" "$i""</td>";print "</tr>"} END{print "</table>"}' >> "${baseDir}"/email_log.html
echo "" >> "${baseDir}"/email_log.html
echo -e "****************************************************\n" >> "${baseDir}"/email_log.html
echo "<p> Possible Mixed Isolates, Defining SNPs called AC=1 </p>" >> "${baseDir}"/email_log.html
cat "${baseDir}"/section2.txt | awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" "$i""</td>";print "</tr>"} END{print "</table>"}' >> "${baseDir}"/email_log.html
echo "" >> "${baseDir}"/email_log.html
echo -e "****************************************************\n" >> "${baseDir}"/email_log.html
cat "${baseDir}"/section3.txt | awk 'BEGIN{print "<table>"} {print "<tr>";for(i=1;i<=NF;i++)print "<td>" "$i""</td>";print "</tr>"} END{print "</table>"}' >> "${baseDir}"/email_log.html
echo "" >> "${baseDir}"/email_log.html
echo -e "****************************************************\n" >> "${baseDir}"/email_log.html
echo "<p> SNP counts: </p>" >> "${baseDir}"/email_log.html
cat s"${baseDir}"/section4.txt | awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' >> "${baseDir}"/email_log.html
echo "" >> "${baseDir}"/email_log.html
echo "****************************************************" >> "${baseDir}"/email_log.html
echo "<p> AC1 called SNPs: </p>" >> "${baseDir}"/email_log.html
cat "${baseDir}"/emailAC1counts.txt | awk 'BEGIN{print "<Body>"} {print "<p style=\"line-height: 40%;\">" $0 "</p>"} END{print "</Body>"}' >> "${baseDir}"/email_log.html
echo "</Body>" >> "${baseDir}"/email_log.html
echo "</html>" >> "${baseDir}"/email_log.html

rm "${baseDir}"/section1.txt
rm "${baseDir}"/section2.txt
rm "${baseDir}"/section3.txt
rm "${baseDir}"/section4.txt
rm "${baseDir}"/section5.txt
rm "${baseDir}"/sectiontime.txt
rm "${baseDir}"/ssection4.txt
rm "${baseDir}"/csection1.txt
# rm -r "${baseDir}"/all_vcfs/starting_files
# zip -r starting_files.zip starting_files
# rm -r starting_files
#rm -r ${FilterDirectory}

#echo "Copy to ${bioinfoVCF}"
#cp -r $PWD ${bioinfoVCF}
fileName=$(basename "$0")

if [ "$mflag" ]; then
    email_list="marc-olivier.duceppe@inspection.gc.ca"
    echo "vcftofasta.sh completed" > "${baseDir}"/mytempfile.txt
    cat "${baseDir}"/mytempfile.txt | mail -s "vcftofasta.sh completed subject" -a "${baseDir}"/email_log.html -- "$email_list"
else
    echo "$fileName "$@" completed, See attachment" > "${baseDir}"/mytempfile.txt
    cat "${baseDir}"/mytempfile.txt | mail -s ""$fileName" "$@" completed" -a "${baseDir}"/email_log.html -- "$email_list"
fi
rm "${baseDir}"/mytempfile.txt
rm "${baseDir}"/email_log.html
echo ""
echo "****************************** END ******************************"
echo ""
#
#  Created by Stuber, Tod P - APHIS on 5/3/2014.
#2015-04-20#
