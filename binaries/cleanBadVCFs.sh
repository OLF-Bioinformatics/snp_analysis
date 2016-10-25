#!/bin/bash

#convert to Unix
#remove trailing tabs (whitespaces)
# cat "$i" | sed -e 's/[[:space:]]*$//'
#remove  double double quotes ("")
# cat "$i" | grep -F "\"\"" | sed 's/\"\"/\"/g'
#remove heading and trailing "
# cat "$i" | grep -E "^\"" | sed 's/^\"//g'
# cat "$i" | grep -F '>"' | sed 's/>"/>/g'
for i in *.vcf; do
    cat "$i" | dos2unix | sed -e 's/[[:space:]]*$//' -e 's/\"\"/\"/g' -e 's/^\"//g' -e 's/>"/>/g' > "${i}".tmp
    mv "${i}".tmp "$i"
done

#Remove double quotes in genotype field "1/1:0,46:46:99:1530,137,0"
for i in *.vcf; do
    cat "$i" | grep -E "^#" > "${i}".tmp #header
    cat "$i" | grep -vE "^#" | tr -d "\"" >> "${i}".tmp
    mv "${i}".tmp "$i"
done
