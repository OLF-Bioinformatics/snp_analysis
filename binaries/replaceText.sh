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
folder=$(readlink -e "$1") 
conversionTable=$(readlink -e "$2")


###### CHECKS ######


#test if argument 1 is a folder (and exists) 
[ -d "$folder" ] && [ -n "$folder" ] || \
echo -e "Invalid folder\n\
The folder specified does not exist\n\n\
Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>" || exit 1


#test if argument 1 is a folder, if it exists and if it's not empty
[ $(ls "$folder" | wc -l) -gt 0  ] || \
echo -e "Invalid folder\n\
The folder specified is empty\n\n\
Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>" || exit 1


#test if argument 2 is a valid tab-sperated files with two fileds
[ -e "$conversionTable" ] || \
echo -e "Conversion table file doest not exist\n\
The conversion table should be a tab-separated file with two columns:\n\n\
current_name\t desired_name\n\n\
Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>" || exit 1


#test if argument 2 is a valid tab-sperated files with two fileds
[ $(head -n 1 "$conversionTable" | awk '{print NF}') -eq 2 ] || \
echo -e "Wrong conversion table format\n\
The conversion table should be a tab-separated file with two columns:\n\n\
current_name\t desired_name\n\n\
Usage: bash rename.sh <folderWithFilesToRename> <conversionTable>" || exit 1


##########################


#The renaming loop
for i in $(find "$folder" -type f); do
	awk 'NR==FNR {a[$1]=$2;next} {for ( i in a) gsub(i,a[i])}1' "$conversionTable" "$i" > "${i}".renamed
    mv  "${i}".renamed "$i"
done

