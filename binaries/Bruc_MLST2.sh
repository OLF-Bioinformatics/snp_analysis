#!/bin/bash


#####################
#                   #
#   Sofware paths   #
#                   #
#####################


picardPath=""${HOME}"/prog/picard-tools-2.4.1/picard.jar"
GATKPath=""${HOME}"/prog/gatk/GenomeAnalysisTK.jar"

sendreport=""${reports}"/mlstCheck.txt"
mlst=""${dependents}"/Brucella_MLST"

#Output folder for mlst files
mlstResults=""${baseDir}"/"${n}"/mlst"

#create folder if doesn't exist
[ -d "$mlstResults" ] || mkdir -p "$mlstResults"


######################
#                    #
#     Resources      #
#                    #
######################


#computer performance
cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"


##########################
#                        #
#   Reference indexing   #
#                        #
##########################


#Reference
# Grab reads and reference and place them in variables
ref="${mlst}"/ST1-MLST.fasta
refName="$(basename "$ref")"
refNameNoExt="${refName%.*}"

#Reference sequence indexing
[ -s "${mlst}"/"${refName}".bwt ] && [ -s "${mlst}"/"${refName}".pac ] && [ -s "${mlst}"/"${refName}".ann ] && [ -s "${mlst}"/"${refName}".amb ] && [ -s "${mlst}"/"${refName}".sa ] || bwa index "$ref"
[ -s "${mlst}"/"${refName}".fai ] || samtools faidx "$ref"
[ -s "${mlst}"/"${refNameNoExt}".dict ] || java "$memJava" -jar "$picardPath" CreateSequenceDictionary REFERENCE="$ref" OUTPUT="${mlst}"/"${refNameNoExt}".dict #index if does not exist


#######################
#                     #
#      Alignment      #
#                     #
####################### 


#Align the reads on reference sequence
echo "Aligning reads on reference Brucella sequence, compressing alignment file (sam to bam) and sorting output bam (by chromosome and by position)"

rg="@RG\tID:"${n}"\tPL:ILLUMINA\tPU:"${n}"_RG1_UNIT1\tLB:"${n}"_LIB1\tSM:"${n}""

bwa mem -M -r 1 -t "$cpu" -R "$rg" "$ref" "$r1" "$r2" | \
    sambamba view -t "$cpu" -f bam -h -S /dev/stdin | \
    sambamba sort -t "$cpu" -o "${mlstResults}"/"${n}".sorted.bam /dev/stdin # mapped and unmapped reads

#Detect SNPs
java "$memJava" -jar "$GATKPath" -T UnifiedGenotyper \
    -R "$ref" \
    -I "${mlstResults}"/"${n}".sorted.bam \
    -o "${mlstResults}"/"${n}".vcf \
    -glm BOTH \
    -out_mode EMIT_ALL_SITES \
    -nct "$cpu"

#Cleanup
rm "${mlstResults}"/"${n}".sorted.bam*
rm "${mlstResults}"/"${n}".vcf.idx

# Position 1629 was too close to the end of glk sequence.  Reads would not assemble properly to call possilbe SNP, therefore 100 bases of the gene were added.
# Because of this all positions beyond this point are 100 more.  Same with position 1645 and 2693.

pos=('231' '297' '363' '398' '429' '523' '631' '730' '1247' '1296' '1342' '1381' '1648' '1685' '1741' '1754' '2165' '2224' '2227' '2297' '2300' '2344' '2352' '2403' '2530' '2557' '2578' '2629' '3045' '3054' '3118' '3295' '3328' '3388' '3966' '3969' '4167' '4271' '4296' '4893' '4996' '4998' '5058' '5248' '5672' '5737' '5928' '5963' '5984' '5987' '6025' '6045' '6498' '6499' '6572' '6627' '6715' '6735' '6745' '6785' '6810' '6828' '6845' '6864' '6875' '7382' '7432' '7464' '7594' '7660' '7756')
# echo "${pos[@]}"


# awk 'BEGIN{OFS="\t"} $2 ~ /^231$/ {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' "${mlstResults}"/"${n}".vcf > "${mlstResults}"/0031.txt

#Array to strore the Alleles at all the positions
alleles=()

for i in "${pos[@]}"; do
    #if present, print alternative allele, otherwise print reference allele
    alleles+=($(cat "${mlstResults}"/"${n}".vcf \
        | awk -v pos="$i" 'BEGIN{OFS="\t"}
            {
                if ($2 == pos && $5 == ".")
                    { print $4 }
                else if ($2 == pos && $5 != ".")
                    { print $5 }
            }'))
done

snps=$(echo "${alleles[@]}" | tr -d " ")

echo "This is the string of SNPs"
echo "$snps"


#################################################

#check if new MLST SNPs present in sample

cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 201 && $2 <= 789) print $0 }' > "${mlstResults}"/"${i}".section1.txt
cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 1190 && $2 <= 1754) print $0 }' > "${mlstResults}"/"${i}".section2.txt
cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 2155 && $2 <= 2629) print $0 }' > "${mlstResults}"/"${i}".section3.txt
cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 3030 && $2 <= 3499) print $0 }' > "${mlstResults}"/"${i}".section4.txt
cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 3900 && $2 <= 4368) print $0 }' > "${mlstResults}"/"${i}".section5.txt
cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 4769 && $2 <= 5254) print $0 }' > "${mlstResults}"/"${i}".section6.txt
cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 5655 && $2 <= 6076) print $0 }' > "${mlstResults}"/"${i}".section7.txt
cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 6477 && $2 <= 6966) print $0 }' > "${mlstResults}"/"${i}".section8.txt
cat "${mlstResults}"/"${n}".vcf | grep -v "#" | awk '{ if ($2 >= 7367 && $2 <= 7796) print $0 }' > "${mlstResults}"/"${i}".section9.txt

cat "${mlstResults}"/"${i}".section*.txt > "${mlstResults}"/"${n}".sectionall.txt
rm "${mlstResults}"/"${i}".section?.txt

cat "${mlstResults}"/"${n}".sectionall.txt | \
    awk '$2 !~ /^231$|^297$|^363$|^398$|^429$|^523$|^631$|^730$|^1247$|^1296$|^1342$|^1381$|^1648$|^1685$|^1741$|^1754$|^2165$|^2224$|^2227$|^2297$|^2300$|^2344$|^2352$|^2403$|^2530$|^2557$|^2578$|^2629$|^3045$|^3054$|^3118$|^3295$|^3328$|^3388$|^3966$|^3969$|^4167$|^4271$|^4296$|^4893$|^4996$|^4998$|^5058$|^5248$|^5672$|^5737$|^5928$|^5963$|^5984$|^5987$|^6025$|^6045$|^6498$|^6499$|^6572$|^6627$|^6715$|^6735$|^6745$|^6785$|^6810$|^6828$|^6845$|^6864$|^6875$|^7382$|^7432$|^7464$|^7594$|^7660$|^7756$/ {print $0}' \
    > "${mlstResults}"/"${n}".target.vcf

#remove temporary file
rm "${mlstResults}"/"${n}".sectionall.txt



# Find possible SNPs
cat "${mlstResults}"/"${n}".target.vcf \
    | awk '$5 != "." {snps[$2]} END{for(i in snps) print i}' \
    | sort -n \
    > "${mlstResults}"/possibleSNPs.txt

for l in $(cat "${mlstResults}"/possibleSNPs.txt); do
    for i in "${mlstResults}"/*.target.vcf; do
        awk -v x="$l" 'BEGIN{OFS="\t"} $2 ~ "^"x"$" {if ( $5 == "." ) print FILENAME, $4; else print FILENAME, $5}' "$i"
    done > "${mlstResults}"/"${l}".txt
done

if [ -s "${mlstResults}"/possibleSNPs.txt ]; then
    echo "There may be a new MLST SNP to be using."
    echo "See reference at position:"
    echo ""$n" may have an unused SNP:" >> "$sendreport"
    cat "${mlstResults}"/possibleSNPs.txt | tee -a "$sendreport"
else
    rm "${mlstResults}"/possibleSNPs.txt
fi

#remove vcf file. Not needed anymore
rm "${mlstResults}"/*.vcf



###################################################

unset mlstArray
declare -A mlstArray=( \
    ["CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"]="MLST type 01" \
    ["CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"]="MLST type 02" \
    ["CTCCCGTGGCGACCCGAGCGAAGCGGGAAGGCCACGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"]="MLST type 03" \
    ["CTCCCGGGGCGACCCGAGCGAAGCGGGAAGGCCAAGGCGCGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"]="MLST type 04" \
    ["CTCCCGGGGCGACCCGATCGAAGCGGGAAGGCCACGGCGAGTGAGCAGCCGGGCATCTGTCCCGCGGGGTA"]="MLST type 05" \
    ["TTCCTGGGGCAACCCGAGCGAGGCAGGGAGGCCGCGGCTCGTGAGCGGTCGGGCATCTGTCCCGCGGGGTA"]="MLST type 06" \
    ["CTTCCTGGCCGAGCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"]="MLST type 07" \
    ["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCCGTCTCGCGGTGC"]="MLST type 08" \
    ["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGTGGTGCT"]="MLST type 09" \
    ["CTTCCTGGCCGACCCGAGTGAAGGGGGGGGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"]="MLST type 10" \
    ["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCGCGGCTGGGTACCTGTCTCGCGGTGCT"]="MLST type 11" \
    ["CTTCCTGGCCGACCCGAGTGAAGGGGGGAGGCCACGGCGCGTGCTCGGCTGGGTACCTGTCTCGCGGTGCT"]="MLST type 12" \
    ["CCCCCGGGCCGACTCGAGCGAAGCGAAGAGGCCACGGCGCGTGAGTGACCAGGCACCTATCCCACGGGGTA"]="MLST type 13" \
    ["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAGTGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA"]="MLST type 14" \
    ["CCCCCGGGCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA"]="MLST type 15" \
    ["CCCCCGGCCCGACCCGGGCGAAGCGGGGAGGCTACGGTGCGTGAGTGGCCAGGCACCTGTCCCGCGAGGTA"]="MLST type 16" \
    ["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGGTA"]="MLST type 17" \
    ["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACCTGTCCCGCAGGCTA"]="MLST type 18" \
    ["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGACACGGCGCGTGAGTGGCCAGGCACCTGTCCCGCGGGGTA"]="MLST type 19" \
    ["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGTCCCGCAGGGTA"]="MLST type 20" \
    ["CCCCCGGGCCGGCCCAAGCGAAGCGGGGAGGCTACAATGCGTGAGTGGCCAGGCACATGCCCCGCAGGGTA"]="MLST type 21" \
    ["CCCCCGGGCCGACCCGAGCGAGGCGGGGAGGCCACGGCGCGGGAGTGGCCAGACACCTGTCCTGCGGGGTA"]="MLST type 22" \
    ["CCCCCGGGCTGACCCGAGCGAAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"]="MLST type 23" \
    ["CCCCCGGGCTGACCCGAGCGGAACGGGGAAGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"]="MLST type 23x" \
    ["CCCCCGGGCCGACCCGAGCAAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"]="MLST type 24" \
    ["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAGGCACCTGTCCCGCGGGGTA"]="MLST type 25" \
    ["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGGCCACGGCGCGTAAGTGGCCAAGCACCTGTTCCGCGGGGTA"]="MLST type 26" \
    ["CCCCCGGGCCGACCCGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCACCTGTCCCGCGGGGTA"]="MLST type 27" \
    ["CCCTCGGGCCGACCTGAGCGAAGCGGGGAGACCACGGCGCATAAGTGGCCAGGCTCCTGTCCCGCGGGGTA"]="MLST type 28" \
)

typing="${mlstArray["$snps"]}"


#Reporting
if [ -n "$typing" ]; then
    echo ""$n" --> $typing" | tee -a "$sendreport"
else
    echo ""$n" --> No MLST type found" | tee -a "$sendreport"
fi

#
#  Created by Stuber, Tod P - APHIS on 2014-09-02.
#

