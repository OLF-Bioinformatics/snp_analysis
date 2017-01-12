#!/bin/bash



##########################
#                        #
#   script1 -> script2   #
#                        #
##########################


#where is the base output folder of script1
vcfs=""${HOME}"/Desktop/tb2017"
[ -d "$vcfs" ] || mkdir -p "$vcfs" #create it if doesn't exist



#Get the "SNPsZeroCoverage" VCF files from script1 to run script2
for d in $(find "/home/bioinfo/analyses/mbovis_script1_2017" -maxdepth 1 -type d -name "16-*");do
    sampleName=$(basename $d)
    cp "${d}"/variant/"${sampleName}".SNPsZeroCoverage.vcf "$vcfs"
done
