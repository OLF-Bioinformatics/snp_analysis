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
    k=26 \
    in1="$r1" \
    in2="$r2" \
    ref="${dependents}"/oligos.fasta \
    stats=""${baseDir}"/"${n}"/"${n}"_oligos.txt"


# #Name   Reads   ReadsPct
# primertb7   97  0.00472%
# primertb3   90  0.00438%
# primertb2   87  0.00423%
# primertb4   79  0.00384%
# primertb157 78  0.00380%


#Brucella
ab1="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbab1" | cut -f 2)"
ab3="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbab3" | cut -f 2)"
ab5="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbab5" | cut -f 2)"
mel="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbmel" | cut -f 2)"
suis1="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbsuis1" | cut -f 2)"
suis2="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbsuis2" | cut -f 2)"
suis3="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbsuis3" | cut -f 2)"
ceti1="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbceti1" | cut -f 2)"
ceti2="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbceti2" | cut -f 2)"
canis4="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbcanis4" | cut -f 2)"
canis="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -Fw "pbcanis" | cut -f 2)"
ovis="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "pbovis" | cut -f 2)"

#Para
para="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "ppara" | cut -f 2)"

#TB
tbbov="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "primertbbov" | cut -f 2)"
tb2="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "primertb2" | cut -f 2)"
tb3="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "primertb3" | cut -f 2)"
tb4="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "primertb4" | cut -f 2)"
tb5="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "primertb5" | cut -f 2)"
tb6="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "primertb6" | cut -f 2)"
tb7="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "primertb7" | cut -f 2)"
tb157="$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -F "primertb157" | cut -f 2)"

#Set to zero unfound primers
[ -z "$ab1" ] && ab1=0
[ -z "$ab3" ] && ab3=0
[ -z "$ab5" ] && ab5=0
[ -z "$mel" ] && mel=0
[ -z "$suis1" ] && suis1=0
[ -z "$suis2" ] && suis2=0
[ -z "$suis3" ] && suis3=0
[ -z "$ceti1" ] && ceti1=0
[ -z "$ceti2" ] && ceti2=0
[ -z "$canis4" ] && canis4=0
[ -z "$canis" ] && canis=0
[ -z "$ovis" ] && ovis=0

[ -z "$tbbov" ] && tbbov=0
[ -z "$tb2" ] && tb2=0
[ -z "$tb3" ] && tb3=0
[ -z "$tb4" ] && tb4=0
[ -z "$tb5" ] && tb5=0
[ -z "$tb6" ] && tb6=0
[ -z "$tb7" ] && tb7=0
[ -z "$tb157" ] && tb157=0

[ -z "$para" ] && para=0

#Count number of hits
bruccounts=("$ab1" "$ab3" "$ab5" "$mel" "$suis1" "$suis2" "$suis3" "$ceti1" "$ceti2" "$canis4" "$canis" "$ovis")
tbcounts=("$tb157" "$tb7" "$tbbov" "$tb5" "$tb2" "$tb3" "$tb4" "$tb6")
paracounts="$para"

# echo "${bruccounts[@]}"
# echo "${tbcounts[@]}"
# echo "$paracounts"

#convert counts to binary (present or absent)
brucbinary=$(echo "${bruccounts[@]}" | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]")
tbbinary=$(echo "${tbcounts[@]}" | awk '{for(i=1;i<=NF;i++) if ($i >= 1) print 1; else print 0}' | tr -cd "[:print:]")
parabinary=$(echo "$paracounts" | awk '{for(i=1;i<=NF;i++) if ($i >= 10) print 1; else print 0}' | tr -cd "[:print:]")

# echo "$brucbinary"
# echo "$tbbinary"
# echo "$parabinary"


##################
#                #
#   Species ID   #
#                #
##################


#to count supporting reads
sumBruOligos=0
sumTbOligos=0
sumParaOligos=0

catch=""
ID=""

#Brucella

unset bruID
declare -A bruID=( \
    ["011111111111"]="ab1" \
    ["101111111111"]="ab1" \
    ["110111111111"]="ab1" \
    ["111011111111"]="mel" \
    ["111101111111"]="suis1" \
    ["111110111111"]="suis2" \
    ["111111011111"]="suis3" \
    ["111111101111"]="ceti1" \
    ["111111100111"]="ceti1" \
    ["111111110111"]="ceti2" \
    ["111111111011"]="suis4" \
    ["111111111001"]="canis" \
    ["111111111110"]="ovis" \
    )

wait

unset bruCatch
declare -A bruCatch=( \
    ["011111111111"]="Brucella abortus bv 1, 2 or 4" \
    ["101111111111"]="rucella abortus bv 3" \
    ["110111111111"]="Brucella abortus bv 5, 6 or 9" \
    ["111011111111"]="Brucella melitensis" \
    ["111101111111"]="Brucella suis bv 1" \
    ["111110111111"]="Brucella suis bv 2" \
    ["111111011111"]="Brucella suis bv 3" \
    ["111111101111"]="Brucella ceti bv 1" \
    ["111111100111"]="Brucella ceti bv 1" \
    ["111111110111"]="Brucella ceti bv 2" \
    ["111111111011"]="Brucella suis bv 4" \
    ["111111111001"]="Brucella canis" \
    ["111111111110"]="Brucella ovis" \
    )

wait

#Sample ID for log
echo -e "**********   "$n"   **********\n" | tee -a "${reports}"/oligo_identifier_log.txt

#count reads
checkBru=0
for t in "${bruccounts[@]}"; do
    let checkBru+=t
done

bruFound=0

if [ "$checkBru" -gt 0 ]; then
    echo "Brucella species found"
    bruFound=1

    #number of reads supporting the oligos found by bbduk for Brucella
    bruOligos=($(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -E "^pb" | cut -f 2)) #array

    for s in "${bruOligos[@]}"; do
        let sumBruOligos+=s
    done

    echo "Number of reads supporting the Brucella-specific oligos = "$sumBruOligos"" | tee -a "${reports}"/oligo_identifier_log.txt
    echo "$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -E "^pb" | cut -f 1,2)" | tee -a "${reports}"/oligo_identifier_log.txt
    
    catch="${bruCatch["$brucbinary"]}"
    ID="${bruID["$brucbinary"]}"

    #reporting
    echo "$catch" >> "${baseDir}"/"${n}"/tee_bruc_oligo_identifier_out2.txt
    echo "$brucbinary" >> "${baseDir}"/"${n}"/tee_bruc_oligo_identifier_out2.txt
    echo "${bruccounts[@]}" >> "${baseDir}"/"${n}"/tee_bruc_oligo_identifier_out2.txt

    #no ID
    if [ -z "$ID" ]; then
        catch="Odd isolate. Unexpected SNP pattern."
        echo -e ""$n" Cannot identify a reference genome for Brucella isolate "$n".\noligo_identifier2.sh stats:\n\tOligo counts: "$bruccounts"\n\tBinary: "$brucbinary"" \
            >> "${reports}"/dailyReport.txt
    fi
fi


#Mycobacterium


unset tbID
declare -A tbID=( \
    ["11101111"]="TB1" \
    ["01100111"]="TB2" \
    ["01101011"]="TB3" \
    ["11101011"]="TB3" \
    ["01101111"]="TB4a" \
    ["01101101"]="TB4b" \
    ["11101101"]="TB4b" \
    ["11111111"]="TB5" \
    ["11001111"]="TB6" \
    ["10101110"]="TB7" \
    ["11001110"]="TBBOV" \
    ["11011110"]="TBBOV" \
    ["11001100"]="TBBOV" \
    ["11111110"]="TBBOV" \
    )

wait

unset tbCatch
declare -A tbCatch=( \
    ["11101111"]="Mycobacterium bovis group 1" \
    ["01100111"]="Mycobacterium bovis group 2" \
    ["01101011"]="Mycobacterium bovis group 3" \
    ["11101011"]="Mycobacterium bovis group 3" \
    ["01101111"]="Mycobacterium bovis group 4a" \
    ["01101101"]="Mycobacterium bovis group 4b" \
    ["11101101"]="Mycobacterium bovis group 4b" \
    ["11111111"]="Mycobacterium bovis group 5" \
    ["11001111"]="Mycobacterium bovis group 6" \
    ["10101110"]="Mycobacterium bovis group 7" \
    ["11001110"]="Mycobacterium bovis (TBBOV)" \
    ["11011110"]="Mycobacterium bovis (TBBOV)" \
    ["11001100"]="Mycobacterium bovis (TBBOV)" \
    ["11111110"]="Mycobacterium bovis (TBBOV)" \
    )

wait

#count reads
checkTb=0
for t in "${tbcounts[@]}"; do
    let checkTb+=t
done

tbFound=0

if [ "$checkTb" -gt 0 ]; then
    echo "Mycobacterium bovis complex species found"
    tbFound=1

    #number of reads supporting the oligos found by bbduk for TB
    tbOligos=($(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -E "^primertb" | cut -f 2)) #array

    for s in "${tbOligos[@]}"; do
        let sumTbOligos+=s
    done

    echo "Number of reads supporting the TB-specific oligos = "$sumTbOligos"" | tee -a "${reports}"/oligo_identifier_log.txt
    echo "$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -E "^primertb" | cut -f 1,2)" | tee -a "${reports}"/oligo_identifier_log.txt

    catch="${tbCatch["$tbbinary"]}"
    ID="${tbID["$tbbinary"]}"

    #reporting
    echo "$catch" >> "${baseDir}"/"${n}"/tee_tb_oligo_identifier_out2.txt
    echo "$tbbinary" >> "${baseDir}"/"${n}"/tee_tb_oligo_identifier_out2.txt
    echo "${tbcounts[@]}" >> "${baseDir}"/"${n}"/tee_tb_oligo_identifier_out2.txt

    #no ID
    if [ -z "$ID" ]; then
        catch="Odd isolate. Unexpected SNP pattern."
        echo -e ""$n" Cannot identify a reference genome for Mycobacterium bovis isolate "$n".\noligo_identifier2.sh stats:\n\tOligo counts: "$bruccounts"\n\tBinary: "$brucbinary"" \
            >> "${reports}"/dailyReport.txt
    fi
fi


#Paratyberculosis

#count reads
checkPara=0
for t in "${paracounts[@]}"; do
    let checkPara+=t
done

paraFound=0

if [ "$checkPara" -gt 0 ]; then
    echo "Mycobacterium avium subsp. paratuberculosis found"
    paraFound=1

    #number of reads supporting the oligos found by bbduk for Para
    paraOligos=($(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -E "^pp" | cut -f 2)) #array

    for s in "${paraOligos[@]}"; do
        let sumParaOligos+=s
    done

    echo "Number of reads supporting the paratuberculosis-specific oligo = "$sumParaOligos"" | tee -a "${reports}"/oligo_identifier_log.txt
    echo "$(cat "${baseDir}"/"${n}"/"${n}"_oligos.txt | grep -E "^pp" | cut -f 1,2)" | tee -a "${reports}"/oligo_identifier_log.txt

    if [ "$parabinary" == 1 ]; then
        catch="Mycobacterium avium subsp. paratuberculosis"
        ID="para"
        echo "$catch" >> "${baseDir}"/"${n}"/tee_tb_oligo_identifier_out2.txt
        echo "$parabinary" >> "${baseDir}"/"${n}"/tee_tb_oligo_identifier_out2.txt
        echo "${paracounts[@]}" >> "${baseDir}"/"${n}"/tee_tb_oligo_identifier_out2.txt
    else
        echo "Could not find a match for "$n"" | tee -a "${reports}"/dailyReport.txt
        echo -e ""$n" Unable to find a reference" >> "${reports}"/dailyReport.txt
    fi
fi

#number of target found
targetFound=$(($bruFound+$tbFound+$paraFound))

#No target species found
if [ "$targetFound" == 0  ]; then 
    echo ""$n" needs special attention" >> "${reports}"/dailyReport.txt
    echo ""$n" was not identified as Brucella, Mycobabterium bovis complex or Mycobacterium avium subsp. paratuberculosis." \
    | mail -s ""$n" Needs special attention" marc-olivier.duceppe@inspection.gc.ca
    exit 1
fi  

#More that one target species found
if [ "$targetFound" -gt 1 ]; then
    echo "More that one target species found "$n""

    unset readArray
    declare -A readArray=( \
        ["Brucella"]="$sumBruOligos" \
        ["TB"]="$sumTbOligos" \
        ["Para"]="$sumParaOligos" \
    )
    wait

    # echo "${readArray[@]}"

    #total number of reads matching the oligos for all target species
    total=0
    for t in "${readArray[@]}"; do
        let total+=t
    done
    # echo "$total"

    echo "Estimated fraction of targets in "$n":"
    for k in "${!readArray[@]}"; do
        echo ""$k": $(bc <<< "scale=2; "${readArray["$k"]}"*100/"$total"")%" | sed 's/ \./ 0\./'
    done | sort -t " " -rn -k2,2

    if [ "$sumBruOligos" -gt "$sumTbOligos" ] && [ "$sumBruOligos" -gt "$sumTbOligos" ]; then
        bruFound=1
        tbFound=0
        paraFound=0
        echo "Although other targets were detected, only Brucella reads will be considered for sample "$n""

    elif [ "$sumTbOligos" -gt "$sumBruOligos" ] && [ "$sumTbOligos" -gt "$sumParaOligos" ]; then
        bruFound=0
        tbFound=1
        paraFound=0
        echo "Although other targets were detected, only TB reads will be considered for sample "$n""

    elif [ "$sumParaOligos" -gt "$sumBruOligos" ] [ "$sumParaOligos" -gt "$sumTbOligos" ]; then
        bruFound=0
        tbFound=0
        paraFound=1
        echo "Although other targets were detected, only ParaTB reads will be considered for sample "$n""

    else
        ID=""
        echo "Sample "$n" contain high levels of multiple targets and its analysis will be skipped."
        exit 1
    fi

fi

#Write results to file to be able to retrieve by the parent process
echo -e "ID: "$catch"\n" | tee -a "${reports}"/oligo_identifier_log.txt

echo -e "\
export bruFound="$bruFound"
export tbFound="$tbFound"
export paraFound="$paraFound"
export ID=\""$ID"\"\
" >> "${baseDir}"/"${n}"/variables.txt



#
#  Created by Stuber, Tod P - APHIS on 04/11/2014.
#
