#/bin/bash


#######################
#                     #
#    Report stuff     #
#                     #
#######################


echo "$(date)" >> "${reports}"/coverageReport.txt
echo "" > "${reports}"/dailyReport.txt

# Reset spoligo and bruc mlst check file
echo "WG Spoligo Check" > "${reports}"/spoligoCheck.txt
echo "Brucella MLST Check" > "${reports}"/mlstCheck.txt

#Reset time stamp
dateFile=$(date "+%Y%m%d")

#Prepare report
echo -e "Isolate\tID\tTotal_Bases\tAveDep\t%>Q15" >> "${reports}"/dailyReport.txt
echo ""  > "${reports}"/dailyStats.txt


########################
#                      #
#     Fastq files      #
#                      #
########################


# #populate fastq folder (symbolic links)
# for i in $(find "$fastqPath" -type f | grep -F ".fastq.gz"); do
#     name=$(basename "$i")
#     ln -s "$i" "${baseDir}"/"$name"
# done


##########################
#                        #
#     Main analyses      #
#                        #
##########################


#Find all paired-end files, create a folder with their sample name (everything before the first "_") and move them in
for i in $(find -L "$baseDir" -maxdepth 1 -type f | grep -F "_R1" | grep -F ".fastq.gz"); do
    r1="$i"
    r2="$(echo "$r1" | sed 's/_R1/_R2/')"

    #check if paired-end files are present
    if [ ! -f "$r1" ]; then
        echo "Read 1 file is missing"
        exit 1
    elif [ ! -f "$r2" ]; then
        echo "Read 2 file is missing"
        exit 1
    fi

    name=$(basename "$r1")
    nameNoExt="${name%%.*}" #no file extension (".fastq.gz")
    export n=$(cut -d "_" -f 1 <<< "$nameNoExt")

    [ -d "${baseDir}"/"$n" ] || mkdir -p "${baseDir}"/"$n"
    mv "${baseDir}"/"$n"*.fastq.gz "${baseDir}"/"$n"

    #redefine r1 and r2
    export r1=$(find -L "${baseDir}"/"$n" -type f | grep -F "_R1" | grep -F ".fastq.gz")
    export r2="$(echo "$r1" | sed 's/_R1/_R2/')"

    #Run oligo_identifier2 (bbduk)
    oligo_identifier2.sh | tee "${baseDir}"/"${n}"/tee_oligo_identifier_out1.txt 

    wait

    #returns in "${baseDir}"/"${n}"/variables.txt:
    # bruFound="$bruFound"
    # tbFound="$tbFound"
    # paraFound="$paraFound"
    # ID="$ID"

    #results from  oligo_identifier2.sh
    source "${baseDir}"/"${n}"/variables.txt

    if [ -n "$ID" ]; then
        if [ "$bruFound" -eq 1 ]; then
            echo "Brucella found in sample"
            Bruc_MLST2.sh  "$r1" "$r2"

        elif [ "$tbFound" -eq 1 ]; then
            echo "Mycobacterium tuberculosis found in sample"
            spoligoSpacerFinder2.sh "$r1" "$r2"

        elif [ "$paraFound" -eq 1 ]; then
            echo "Paratuberculosis found in sample"
            #nothing yet

        else
            echo "Shouldn't get here" #because there's a check for that in oligo identifier2!
            exit 1
        fi

        #Run processZips2.sh
        processZips2.sh "$ID" | tee "${baseDir}"/"${n}"/tee_processZips_out.txt

    else
        echo "No Brucella, TB or paraTB signature found in "$n""
        exit 1
    fi
done


#
#  Created by Tod Stuber on 11/05/12, modified 1/22/15.
#
