#!/bin/bash

#check if more than 2 dot in file name
bad=0
for i in *.vcf; do
    mid=$(basename "$i" | cut -d "." -f 2)
    if [ "$mid" != "SNPsZeroCoverage" ]; then 
        let bad+=1
        echo $(basename "$i")
    fi
done

if [ "$bad" -gt 0 ]; then
    echo 'Please rename the previous files so dots (".") are only used for the ".SNPsZeroCoverage.vcf" extension'
    exit 1
fi


#convert to Unix
#remove trailing tabs (whitespaces)
# cat "$i" | sed -e 's/[[:space:]]*$//'
#remove  double double quotes ("")
# cat "$i" | grep -F "\"\"" | sed 's/\"\"/\"/g'
#remove heading and trailing "
# cat "$i" | grep -E "^\"" | sed 's/^\"//g'
# cat "$i" | grep -F '>"' | sed 's/>"/>/g'
total=$(ls *.vcf | wc -l)
counter=0
for i in *.vcf; do
    let counter+=1
    echo -e "Working on "$(basename "$i" | cut -d "." -f 1)"... ("${counter}"/"${total}")"

    cat "$i" | dos2unix | tr "\r" "\n" | sed -e "s/[[:space:]]*$//" -e 's/\"\"/\"/g' -e 's/^\"//g' -e 's/>"/>/g' > "${i}".tmp
    mv "${i}".tmp "$i"
done

#Remove double quotes in genotype field "1/1:0,46:46:99:1530,137,0"
counter=0
for i in *.vcf; do
    let counter+=1
    echo -e "Working on "$(basename "$i" | cut -d "." -f 1)"... ("${counter}"/"${total}")"

    cat "$i" | grep -E "^##" > "${i}".tmp #header
    cat "$i" | grep -E "^#CHROM" | cut -f 1-10 >> "${i}".tmp #Only keep the first 10 data fields, because only one sample per VCF file
    cat "$i" | grep -vE "^#" | cut -f 1-10 | tr -d "\"" >> "${i}".tmp #Only keep the first 10 data fields
    mv "${i}".tmp "$i"
done

#Change sample name inside vcf file to match the file name (before the first dot ".")
for i in ./*.vcf; do
    vcfName=$(cat "$i" | grep -E "^#CHROM" | cut -f 10)
    fileName=$(basename "$i" | cut -d "." -f 1)

    if [ "$fileName" != "$vcfName" ]; then
        echo -e ""$fileName"\t"$vcfName""
        sed -i "s/"$vcfName"/"$fileName"/" "$i"
    fi
done

#Sort VCF files and remove duplicates
counter=0
for i in ./*.vcf; do
    let counter+=1
    echo -e "Working on "$(basename "$i" | cut -d "." -f 1)"... ("${counter}"/"${total}")"
    
    (cat "$i" | grep -E "^#"; cat "$i" | grep -vE "^#" | sort -k1,1d -uk2,2n;) > "${i}".tmp
    mv "${i}".tmp "$i"
done
