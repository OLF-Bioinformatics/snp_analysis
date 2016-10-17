#!/bin/bash



##########################
#                        #
#   script1 -> script2   #
#                        #
##########################


#where is the base output folder of script1
vcfs=""${HOME}"/Desktop/vcf_mbovisCAN"
[ -d "$vcfs" ] || mkdir -p "$vcfs" #create it if doesn't exist

#Get the "SNPsZeroCoverage" VCF files from script1 to run script2
for d in $(find "${HOME}"/analyses/mbovis_script1 -maxdepth 1 -type d | grep -F "MBWGS");do
    sampleName=$(basename $d)
    cp "${d}"/variant/"${sampleName}".SNPsZeroCoverage.vcf "$vcfs"
done
