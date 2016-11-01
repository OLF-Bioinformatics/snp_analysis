#!/bin/bash

#######################################################################################################

#Rename all files in a folder from a list wich contains 2 tab-separated fields.
#First column is the actual name of the file
#Second column is the new name of the file

# Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>"

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
if [ "$(head -n 1 "$conversionTable" | awk '{print NF}')" -ne 2 ]; then
    printf ""${RED}"Wrong conversion table format\n"${GREEN}"The conversion table should be a tab-separated file with two columns:\ncurrent_name\t desired_name\n"${BLUE}"Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>\n"${NC}""
    exit 1
fi

##########################



#Do hash with conversion table
declare -A myArray=()

while IFS= read -r line || [[ -n "$line" ]]; do
    line="$(sed -e 's/ /\t/' -e 's/\r//g' -e 's/\n//g' <<< "$line")" #transform the space output field separator from read into tabs, remove carriage return
    key="$(cut -f 1 <<< "$line")"
    value="$(cut -f 2 <<< "$line")"
    myArray["${key}"]="${value}"
done < "$conversionTable"

wait

#readarray myArray < "$conversionTable"
# echo "${!myArray[@]}" #keys
# echo "${myArray[@]}" #values

#The renaming loop
#This can't handle the spaces in the path/names...
##########################################
# for i in "$(find "$folder" -type f)"; do
#     echo "$i"
#     pathPart="$(dirname "$i")"
#     echo "$pathPart"
# done
##########################################

find "$folder" -type f | while read i; do
    # echo "$i"
    pathPart="$(dirname "$i")"
    # echo "$pathPart"
    oldName="$(basename "$i")"
    # echo "$oldName"

    #for each file, check if a part of the name matches on
    for j in "${!myArray[@]}"; do
        # echo "$j"

        if [ "$(echo "$oldName" | grep "$j")" ]; then
            newName="$(echo "$oldName" | sed "s/"$j"/"${myArray["$j"]}"/")"
            fullNewName=""${pathPart}"/"${newName}""

            if [ -e "$rename" ]; then
                echo "Cannot rename "$oldName" to "$newName", file already exists. Skipping"
                continue
                # exit 1
            fi

            echo ""$i" -> "$fullNewName""
            mv "$i" "$fullNewName"
            wait
        fi
    done
    wait
done
