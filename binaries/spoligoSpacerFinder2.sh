#!/bin/bash


######################
#                    #
#     Resources      #
#                    #
######################


#computer performance
cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"


########################
#                      #
#   Counting matches   #
#                      #
########################


bbduk.sh "$memJava" \
    k=25 \
    in1="$r1" \
    in2="$r2" \
    ref="${dependents}"/spoligos.fasta \
    stats=""${baseDir}"/"${n}"/"${n}"_spoligos.txt"


for i in {01..43}; do
    #get the count vaalue from the bbduk output
    eval sp"${i}"=$(cat "${baseDir}"/"${n}"/"${n}"_spoligos.txt | grep -F "sp"${i}"" | cut -f 2)

    #if didn't find the primer, set value to zero
    [ -z $(eval echo \$sp"${i}") ] && eval sp"${i}"=0

    # count=$(eval echo \$sp"${i}")
    # echo "sp$i=$count"
done


#####################
#                   #
#   Spoligotyping   #
#                   #
#####################

spacers=()
for s in {01..43}; do
    spacers+=("spacer"${s}"")
done
echo "${spacers[@]}" > "${baseDir}"/"${n}"/"${n}".spacers.txt

values=()
for v in {01..43}; do
    values+=($(eval echo \$sp"${v}"))
done
echo "${values[@]}" >> "${baseDir}"/"${n}"/"${n}".spacers.txt

#convert to binary
mybinaries=($(cat "${baseDir}"/"${n}"/"${n}".spacers.txt | awk 'NR==2 {for(i=1;i<=NF;i++) if ($i >= 5) print 1; else print 0}' | tr -cd "[:print:]" | fold -w3))
# echo "${mybinaries[@]}"

#convert binary to octal
octalCode=()
for b in "${mybinaries[@]}"; do
    octalCode+=($((2#$b)))
done
# echo "${octalCode[@]}"


#################
#               #
#   Reporting   #
#               #
#################


WGSpoligo=$(echo "${octalCode[@]}" | tr -d " ")
echo -e "Whole genome based spoligotyping\n" | tee -a "${reports}"/spoligoCheck.txt "${reports}"/spoligoCheck_all.txt
echo ""$n": "$WGSpoligo"" | tee -a "${reports}"/spoligoCheck.txt
# echo "WGSpoligo: "$WGSpoligo"" | tee -a "${reports}"/spoligoCheck.txt

#Make FileMaker file
# dateFile=$(date "+%Y%m%d")
# printf "%s\t%s\n" "$n" "$WGSpoligo" >> "${reports}"/"${dateFile}"_FileMakerSpoligoImport.txt

# Add infor to spoligoCheck_all.txt
# echo "<----- "$n" ----->" >> "${reports}"/spoligoCheck_all.txt
# echo "WGSpoligo: "$WGSpoligo"" >> "${reports}"/spoligoCheck_all.txt


#
#  Created by Stuber, Tod P - APHIS on 03/07/2013.
#
