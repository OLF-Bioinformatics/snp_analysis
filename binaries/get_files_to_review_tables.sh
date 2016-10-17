#!/bin/bash



##########################
#                        #
#   script2 -> script2   #
#                        #
##########################


#Various output files needed to review the SNP quality (manual curration)
#Recursive bad quality SNPs are added to the Filter table used by script2
#Script2 has tp be rerun after editing of the filter file
#As a result, we get a better tree and a higher resolution

#Where to put the files needed to review the SNP table
reportFolder=""${HOME}"/Desktop/Olaga_vcf_mbovisCAN"
[ -d "$reportFolder" ] || mkdir -p "$reportFolder" #create it if doesn't exist

#copy by groups
#When lots of samples are analyzed, the table with all the SNPs becomes too hard to work with
#Plus excel is somewhat unhappy with the very large number of columns.
#So better splitting it up in groups
for i in $(find "${HOME}"/analyses/mbovisCAN_script2v2a/all_groups -maxdepth 1 -type d | grep -F "Group");do
    group=$(basename "$i") #group name
    samples=($(find "${i}"/fasta -type f | grep -F ".fas" | grep -vF "root")) #samples in the group

    # echo "${samples[@]}" | tr " " "\n" #Debug, print array content on screen. Replace spaces by carriage return.
# done
    gFolder="${reportFolder}"/"${group}" #output folder for the group
    [ -d "$gFolder" ] || mkdir -p "$gFolder" #create it if doesn't exist


    #copy bam, bai and vcf from script 1
    for s in "${samples[@]}"; do
        sampleName=$(basename "$s")
        sampleNameNoExt="${sampleName%.*}"
        # echo $sampleNameNoExt

        #realigned+recalibrated bam (All reads)
        cp "${HOME}"/analyses/mbovis_script1/"${sampleNameNoExt}"/realigned/"${sampleNameNoExt}"_realigned_recalibrated.bam "$gFolder" #bam
        cp "${HOME}"/analyses/mbovis_script1/"${sampleNameNoExt}"/realigned/"${sampleNameNoExt}"_realigned_recalibrated.bai "$gFolder" #bam index file (not sure it's needed for IGV)

        #bam from from HaplotypeCaller bamout (Only read used to call variants)
        cp "${HOME}"/analyses/mbovis_script1/"${sampleNameNoExt}"/variant/"${sampleNameNoExt}"_hapotypes.bam "$gFolder" #bam
        cp "${HOME}"/analyses/mbovis_script1/"${sampleNameNoExt}"/variant/"${sampleNameNoExt}"_hapotypes.bai "$gFolder" #bam index file (not sure it's needed for IGV)

        #vcf
        cp "${HOME}"/analyses/mbovis_script1/"${sampleNameNoExt}"/variant/"${sampleNameNoExt}".vcf "$gFolder" #vcf (might need the vcf index file too?)
    done


    #copy table and tree file from script2
    #not the best way to loop through all the files, but works!
    for j in $(find "${i}"/fasta -type f); do
        #Only really want to get those two files per group
        a=$(echo "$j" | grep -F ".organizedTable.txt")
        b=$(echo "$j" | grep -F "RAxML_bestTree.")

        #copy it if it's the files that I'm interested in
        if [ -n "$a" ]; then
           cp "$j" "$gFolder"
        elif [ -n "$b" ]; then
           cp "$j" "$gFolder"
        fi
    done
done
