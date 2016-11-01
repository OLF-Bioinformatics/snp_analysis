#!/bin/bash


#######################################################################################################

#batch replace text in a file using a tab-delimited dictionary file:

#First column is the word to replace (it can be a sentence?) (watchout for special characters?)
#Second column is the new word

#Note that the find command looks for file recursively (will search in all subfolder if exist).
#Add "-maxdepth 1" before "-type f" in the find command of the for loop to don't look into subfolders.


# Usage: bash replaceText.sh <folderWithFilesToRename> <conversionTable>"

#######################################################################################################


#I/O
# folder=$(readlink -e "$1") 
# conversionTable=$(readlink -e "$2")

# folder="$1"
# conversionTable="$2"


folder="/media/bioinfo/SANDISK128/olga/mbovis_16_23"
conversionTable="/media/3tb_hdd/data/Mycobaterium_bovis/TB_tree_transtate_fuller.txt"


###### CHECKS ######

# lack        0;30     Dark Gray     1;30
# Red          0;31     Light Red     1;31
# Green        0;32     Light Green   1;32
# Brown/Orange 0;33     Yellow        1;33
# Blue         0;34     Light Blue    1;34
# Purple       0;35     Light Purple  1;35
# Cyan         0;36     Light Cyan    1;36
# Light Gray   0;37     White         1;37

RED='\033[1;31m'
BLUE='\033[1;34m'
GREEN='\033[1;32m'
NC='\033[0m' # No Color

#test if argument 1 is a folder (and exists)
if [ ! -d "$folder" ]; then
    printf ""${RED}"Invalid folder\n"${GREEN}"The folder specified does not exist\n"${BLUE}"Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>\n"${NC}""
    exit 1
fi

#test if argument 1 is a folder, if it exists and if it's not empty
if [ $(ls "$folder" | wc -l) -eq 0 ]; then
    printf ""${RED}"The specified folder is empty\n"${BLUE}"Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>\n"${NC}""
    exit 1
fi

#test if argument 2 is a valid tab-sperated files with two fileds
if [ ! -f "$conversionTable" ]; then
    printf ""${RED}"Provided conversion table file doest not exist\n"${BLUE}"Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>\n"${NC}""
    exit 1
fi

#test if argument 2 is a valid tab-sperated files with two fileds
if [ $(head -n 1 "$conversionTable" | awk '{print NF}') -ne 2 ]; then
    printf ""${RED}"Wrong conversion table format\n"${GREEN}"The conversion table should be a tab-separated file with two columns:\ncurrent_name\t desired_name\n"${BLUE}"Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>\n"${NC}""
    exit 1
fi

##########################


#The text replacing loop
#To change names in SNP table
# for i in $(find "$folder" -type f); do #doesn't work when path/file names have spaces
find "$folder" -type f | grep -F "organized_table" | while read i; do
    echo "Renaming "$i""
    awk 'NR==FNR {a[$1]=$2;next} {for ( i in a) gsub(i,a[i])}1' "$conversionTable" "$i" > "${i}".renamed
    mv  "${i}".renamed "$i"
done

