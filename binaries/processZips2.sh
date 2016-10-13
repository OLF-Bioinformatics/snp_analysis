#!/bin/bash


######################
#                    #
#    User Defined    #
#                    #
######################


#where programs are all installed
prog=""${HOME}"/prog"
scripts=""${HOME}"/scripts"
gatk=""${prog}"/gatk/GenomeAnalysisTK.jar"
picard=""${prog}"/picard-tools-2.4.1/picard.jar"
igvtools=""${prog}"/IGVTools/igvtools.jar"


#######################
#                     #
#    Data Stucture    #
#                     #
#######################


#Directory structure
sampleDir=""${baseDir}"/"${n}""
logs=""${sampleDir}"/logs"
fastq=""${sampleDir}"/fastq"
trimmed=""${sampleDir}"/trimmed"
merged=""${sampleDir}"/merged"
aligned=""${sampleDir}"/aligned"
realigned=""${sampleDir}"/realigned"
variant=""${sampleDir}"/variant"
tree=""${sampleDir}"/tree"

#Create required directories if don't exist
[ -d "$logs" ] || mkdir -p "$logs"
[ -d "$fastq" ] || mkdir -p "$fastq"
[ -d "$trimmed" ] || mkdir -p "$trimmed"
[ -d "$merged" ] || mkdir -p "$merged"
[ -d "$aligned" ] || mkdir -p "$aligned"
[ -d "realigned" ] || mkdir -p "$realigned"
[ -d "$variant" ] || mkdir -p "$variant"
[ -d "$tree" ] || mkdir -p "$tree"


#####################
#                   #
#     Resources     #
#                   #
#####################


#computer performance
cpu=$(nproc) #total number of cores
mem=$(($(grep MemTotal /proc/meminfo | awk '{print $2}')*85/100000000)) #85% of total memory in GB
memJava="-Xmx"$mem"g"


#######################
#                     #
#   Initiating logs   #
#                     #
#######################


#Date
echo -e "$(date)" | tee "${logs}"/log.txt
echo -e "User: $(whoami)" | tee -a "${logs}"/log.txt
echo -e "Processors: "$cpu"" | tee -a "${logs}"/log.txt
echo -e "Memory: "$mem"G" | tee -a "${logs}"/log.txt


####################
#                  #
#     Trimming     #
#                  #
####################


#Trim the reads with bbmap tool kit (bbduk plugin)
#about twice as fast as trimmomatic
start=$(date +%s)

[ -e "${logs}"/trimming.txt ] && rm "${logs}"/trimming.txt
echo -e "Quality trimming sample "$n"..." | tee -a "${logs}"/log.txt

bbduk.sh "$memJava" \
    in1="$r1" \
    in2="$r2" \
    ref="${prog}"/bbmap/resources/nextera.fa.gz \
    ktrim=r k=23 mink=11 hdist=1 tbo tpe \
    qtrim=lr trimq=10 \
    minlen=64 \
    out1="${trimmed}"/"${n}"_Trimmed_1P.fastq.gz \
    out2="${trimmed}"/"${n}"_Trimmed_2P.fastq.gz \
    pigz=t \
    unpigz=t \
    2> >(tee -a "${logs}"/trimming.txt)

wait

end=$(date +%s)
elapsed=$(($end - $start))
printf "Trimming finished in %dh:%dm:%ds\n" \
    $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a "${logs}"/log.txt


###################
#                 #
#     Merging     #
#                 #
###################


#paired-end read merging
start=$(date +%s)

t1=$(find "$trimmed" -type f | grep -F "fastq.gz" | grep -F "_1P")
t2=$(sed "s/_1P/_2P/" <<< "$t1")

[ -e "${logs}"/merging.txt ] && rm "${logs}"/merging.txt
echo -e "Merging paried-end reads for "$n"..." | tee -a "${logs}"/log.txt

bbmerge.sh "$memJava" \
    in1="$t1" \
    in2="$t2" \
    out="${merged}"/"${n}"_merged.fastq.gz \
    outu1="${merged}"/"${n}"_unmerged_1P.fastq.gz \
    outu2="${merged}"/"${n}"_unmerged_2P.fastq.gz \
    pigz=t \
    ziplevel=5 \
    unpigz=t \
    2> >(tee -a "${logs}"/merging.txt)

wait

end=$(date +%s)
elapsed=$(($end - $start))
printf "Merging finished in %dh:%dm:%ds\n" \
    $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a "${logs}"/log.txt

#remove trimmed files
rm "${trimmed}"/*


###########################
#                         #
#   Reference selection   #
#                         #
###########################


unset bruRefGenome
declare -A bruRefGenome=( \
    ["ab1"]=""${dependents}"/Brucella_abortus/NC_00693c.fasta" \
    ["mel"]=""${dependents}"/Brucella_melitensis/BmelitensisM5-90.fasta" \
    ["suis1"]=""${dependents}"/Brucella_suis_bv1/NC_01725c.fasta" \
    ["suis2"]=""${dependents}"/Brucella_suis_bv2/Bsuisbv2-94-11.fasta" \
    ["suis3"]=""${dependents}"/Brucella_suis_bv3/B-REF-BS3-686.fasta" \
    ["suis4"]=""${dependents}"/Brucella_suis_bv4/B-REF-BS4-40.fasta" \
    ["canis"]=""${dependents}"/Brucella_canis/BcanisATCC23365.fasta" \
    ["ceti1"]=""${dependents}"/Brucella_ceti-grp1/Bceti1Cudo.fasta" \
    ["ceti2"]=""${dependents}"/Brucella_ceti-grp2/Bceti2-TE10759.fasta" \
    ["ovis"]=""${dependents}"/Brucella_ovis/BovisATCC25840.fasta" \
    ["TB1"]=""${dependents}"/TB1/NC_017528.fasta" \
    ["TB2"]=""${dependents}"/TB2/NC_021251.fasta" \
    ["TB3"]=""${dependents}"/TB3/NC_021193it3-readreference.fasta" \
    ["TB4a"]=""${dependents}"/TB4a/NC002755.fasta" \
    ["TB4a"]=""${dependents}"/TB4b/NC018143.fasta" \
    ["TB5"]=""${dependents}"/TB5/APKD01000001.fasta" \
    ["TB6"]=""${dependents}"/TB6/NC_015758.fasta" \
    ["TBBOV"]=""${dependents}"/Mycobacterium_bovis/NC_002945.fasta" \
    ["para"]=""${dependents}"/paraTB/NC_002944.fasta" \
)

wait


unset bruRefVCF
declare -A bruRefVCF=( \
    ["ab1"]=""${dependents}"/Brucella_abortus/NC_00693cHighestQualitySNPs.vcf" \
    ["mel"]=""${dependents}"/Brucella_melitensis/melHighestQualitySNPs.vcf" \
    ["suis1"]=""${dependents}"/Brucella_suis_bv1/NC_01725cHighestQualitySNPs.vcf" \
    ["suis2"]=""${dependents}"/Brucella_suis_bv2/suis2HighestQualitySNPs.vcf" \
    ["suis3"]=""${dependents}"/Brucella_suis_bv3/suis3HighestQualitySNPs.vcf" \
    ["suis4"]=""${dependents}"/Brucella_suis_bv4/suis4HighestQualitySNPs.vcf" \
    ["canis"]=""${dependents}"/Brucella_canis/canisHighestQualitySNPs.vcf" \
    ["ceti1"]=""${dependents}"/Brucella_ceti-grp1/ceti1HighestQualitySNPs.vcf" \
    ["ceti2"]=""${dependents}"/Brucella_ceti-grp2/ceti2HighestQualitySNPs.vcf" \
    ["ovis"]=""${dependents}"/Brucella_ovis/BovisATCC25840HighestQualitySNPs.vcf" \
    ["TB1"]=""${dependents}"/TB1/HQ-NC_017528.vcf" \
    ["TB2"]=""${dependents}"/TB2/HQ-NC021251.vcf" \
    ["TB3"]=""${dependents}"/TB3/13-7575-highqualitysnps.vcf" \
    ["TB4a"]=""${dependents}"/TB4a/HQ-NC002755.vcf" \
    ["TB4a"]=""${dependents}"/TB4b/HQ-NC018143.vcf" \
    ["TB5"]=""${dependents}"/TB5/HQ-16-2185-11.vcf" \
    ["TB6"]=""${dependents}"/TB6/HQ-NC015758.vcf" \
    ["TBBOV"]=""${dependents}"/Mycobacterium_bovis/HighestQualitySNPs.vcf" \
    ["para"]=""${dependents}"/paraTB/HQ-NC002944.vcf" \
)

wait

#results from  oligo_identifier2.sh
source "${sampleDir}"/variables.txt

wait

genome="${bruRefGenome["$ID"]}"
hqs="${bruRefVCF["$ID"]}"

if [ -z "$genome" ] || [ -z "$hqs" ]; then
    echo "Could not assign reference genome and VCF files to "$n". Aborting SNP calling."
    exit 1
fi


##########################
#                        #
#   Reference indexing   #
#                        #
##########################


#Create dictionary file for the reference genome (needed by gatk)
if [ -e "${genome%.*}".dict ]; then #check if indexing already done
    echo -e "Reference genome $(basename "$genome") already indexed. Skipping this step."
else
    echo -e "Indexing reference genome $(basename "$genome") with Picard tools..."
    java "$memJava" -jar "$picard" CreateSequenceDictionary \
        R="$genome" \
        O="${genome%.*}".dict #have to strip off the ".fasta" extension and replace it by ".dict"
fi

wait

#Indexing reference VCF file (High Quality SNPs)
if [ -e "${hqs}".tbi ]; then #check if indexing already done
    echo -e "High quality SNPs for $(basename "$genome") already indexed. Skipping this step."
else
    echo -e "Indexing High quality SNPs for $(basename "$genome")..."
    tabix "$hqs" #careful, the vcf file must be compressed with bgzip, not pigz
    if [ -e "${hqs}".tbi ]; then #if didn't work
        echo -e "The VCF file must be compressed with bgzip to be compatible with tabix indexing."
        exit 1
    fi
fi

wait

#index reference genome for samtools
if [ -e "${genome}".fai ]; then #check if indexing already done
    echo -e "Reference genome $(basename "$genome") already indexed. Skipping this step."
else
    echo -e "Indexing reference genome $(basename "$genome") with Samtools..."
    samtools faidx "$genome"
fi

wait

#index reference genome for bwa
if [ -e "${genome}".sa ] && [ -e "${genome}".amb ] \
&& [ -e "${genome}".ann ] && [ -e "${genome}".pac ] \
&& [ -e "${genome}".bwt ]; then #check if indexing already done
    echo -e "Reference genome $(basename "$genome") already indexed. Skipping this step."
else
    echo -e "Indexing reference genome $(basename "$genome") with bwa..."
    bwa index "$genome"
fi

wait


#######################
#                     #
#      Alignment      #
#                     #
####################### 


#Map reads to referenced genome
start=$(date +%s)

echo -e "Aligning "$n" pre-processed reads on $(basename "$genome")..." | tee -a "${logs}"/log.txt

#unmerged
mu1=$(find "$merged" -type f | grep -F "fastq.gz" | grep -F "_1P")
mu2=$(sed "s/_1P/_2P/" <<< "$mu1")

#merged
m=$(find "$merged" -type f | grep -F "fastq.gz" | grep -F "_merged")

#read group
rg="@RG\tID:"${n}"\tCN:OLF\tLB:Nextera\tPL:ILLUMINA\tPM:MiSeq\tSM:"$n""

#map merged reads (single end)

#Align trimmed reads to reference genome (SAM output with Read Group information),
#convert SAM to BAM, filter out unmapped, sort, remove duplicates and index BAM
#https://github.com/lomereiter/sambamba/wiki/%5Bsambamba-view%5D-Filter-expression-syntax
bwa mem -t "$cpu" -r 1 -a -M -R "$rg" "$genome" "$m" | \
    sambamba view -t "$cpu" -f bam -h -F "not (unmapped)" -S /dev/stdin | \
    sambamba sort -t "$cpu" -o "${aligned}"/"${n}"_merged.bam /dev/stdin

wait

sambamba markdup -r -t "$cpu" "${aligned}"/"${n}"_merged.bam "${aligned}"/"${n}"_merged_nodup.bam

wait

#map umerged (paried-end)
bwa mem -t "$cpu" -r 1 -a -M -R "$rg" "$genome" "$mu1" "$mu2" | \
    sambamba view -t "$cpu" -f bam -h -F "not (unmapped)" -S /dev/stdin | \
    sambamba sort -t "$cpu" -o "${aligned}"/"${n}"_unmerged.bam /dev/stdin

wait

sambamba markdup -r -t "$cpu" "${aligned}"/"${n}"_unmerged.bam "${aligned}"/"${n}"_unmerged_nodup.bam

wait

#merge bam files
# Usage: sambamba-merge [options] <output.bam> <input1.bam> <input2.bam> [...]
sambamba merge -t "$cpu" \
    "${aligned}"/"${n}"_all.bam \
    "${aligned}"/"${n}"_merged_nodup.bam \
    "${aligned}"/"${n}"_unmerged_nodup.bam

wait

#elapsed time
end=$(date +%s)
elapsed=$(($end - $start))
printf "Mapping with bwa finished in %dh:%dm:%ds\n" \
    $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a "${logs}"/log.txt

#remove temp bam files
find "$aligned" -maxdepth 1 -type f | grep -vF "all" | xargs rm

#Remove merged fastq files
rm "${merged}"/*


#########################
#                       #
#      Realignment      #
#                       #
######################### 


start=$(date +%s)
echo -e "Running recommended workflows for variant discovery analysis with GATK..." | tee -a "${logs}"/log.txt

#Realign indels - step 1
java "$memJava" -jar "$gatk" -T RealignerTargetCreator \
    -R "$genome" \
    -I "${aligned}"/"${n}"_all.bam \
    -filterNoBases \
    -o "${realigned}"/"${n}".intervals \
    -nt "$cpu"

wait

#Realign indels - step 2
java "$memJava" -jar "$gatk" -T IndelRealigner \
    -R "$genome" \
    -I "${aligned}"/"${n}"_all.bam \
    -filterNoBases \
    -targetIntervals "${realigned}"/"${n}".intervals \
    -o "${realigned}"/"${n}"_realigned.bam

wait

#remove old aligned files
rm "${aligned}"/*

#Recalibrate bases
#known sites obtaines by running bwa and mpileup on all the samples
#then, filtering with vcftools in such a way to keep only the very high quality SNPs
# MQ 60, minQ 999, max-missing 1, maf 0.1
java "$memJava" -jar "$gatk" -T BaseRecalibrator \
    -R "$genome" \
    -I "${realigned}"/"${n}"_realigned.bam \
    -filterNoBases \
    -knownSites "$hqs" \
    -o "${realigned}"/"${n}"_recal_data.table \
    --maximum_cycle_value 1000 \
    -nct "$cpu"

wait

#Print reads
java "$memJava" -jar "$gatk" -T PrintReads \
    -R "$genome" \
    -I "${realigned}"/"${n}"_realigned.bam \
    -filterNoBases \
    -BQSR "${realigned}"/"${n}"_recal_data.table \
    -o "${realigned}"/"${n}"_realigned_recalibrated.bam \
    -nct "$cpu"

wait

#Collect Depth of coverage info for every 
java "$memJava" -jar "$gatk" -T DepthOfCoverage \
    -R "$genome" \
    -I "${realigned}"/"${n}"_realigned_recalibrated.bam  \
    -o "${realigned}"/"${n}"_coverage.txt \
    --omitIntervalStatistics \
    --omitLocusTable \
    --omitPerSampleStats \
    -nt "$cpu"

wait

# #Coverage stats by bbmap (20s vs 8min!) -> gives different values than GATK
# pileup.sh "$memJava" \
#     in="${realigned}"/"${n}"_realigned_recalibrated.bam \
#     out="${sampleDir}"/"${n}"_averageCoverage.txt \
#     basecov="${sampleDir}"/"${n}"_perBaseCoverage.txt


#######################
#                     #
#   Variant Calling   #
#                     #
#######################


#The vcf and the bam files output here are the ones used in IGV for manual SNP curation
#$n.hapreadyAll.vcf
java "$memJava" -jar "$gatk" -T HaplotypeCaller \
    -R "$genome" \
    -I "${realigned}"/"${n}"_realigned_recalibrated.bam \
    --bamOutput "${variant}"/"${n}"_haplotypes.bam \
    --dontUseSoftClippedBases \
    --allowNonUniqueKmersInRef \
    -o "${variant}"/"${n}".vcf

wait

#index variant file for IGV
java "$memJava" -jar "$igvtools" index \
    "${variant}"/"${n}".vcf

wait

#elapsed time
end=$(date +%s)
elapsed=$(($end - $start))
printf "GATK best practices workflow finished in %dh:%dm:%ds\n" \
    $(($elapsed/3600)) $(($elapsed%3600/60)) $(($elapsed%60)) | tee -a "${logs}"/log.txt


#########################
#                       #
#   Variant filtering   #
#                       #
#########################


#Keep only the SNPs
cat "${variant}"/"${n}".vcf | awk '/#/ || $4 ~ /^[ATGC]$/ && $5 ~ /^[ATGC]$/ {print $0}' > "${variant}"/"${n}"_snps.vcf


#########################
#                       #
#   VCF manipulations   #
#                       #
#########################


#Split header lines from position calls
cat "${variant}"/"${n}"_snps.vcf | grep "#" > "${variant}"/"${n}".header
cat "${variant}"/"${n}"_snps.vcf | grep -v "#" > "${variant}"/"${n}".body

#SNP positons that will be used
cat "${variant}"/"${n}".body | awk '{print $1 "%" $2}' > "${variant}"/"${n}".calledSNPpositions

#Zero coverage positions
cat "${realigned}"/"${n}"_coverage.txt | awk 'BEGIN {FS="[:\t]"} $3 == 0 {print $1 "%" $2}' > "${variant}"/"${n}".zeroCoveragePositions

#Remove zero coverage positions that will be are in $n.hapreadyOnlySNPs.vcf
cat "${variant}"/"${n}".calledSNPpositions "${variant}"/"${n}".zeroCoveragePositions \
    | sort | uniq -d \
    > "${variant}"/"${n}".duplicates

cat "${variant}"/"${n}".zeroCoveragePositions "${variant}"/"${n}".duplicates \
    | sort | uniq -u \
    > "${variant}"/"${n}".keepTheseZeroCovPositions

zeroposition=$(cat "${variant}"/"${n}".keepTheseZeroCovPositions | grep -c ".*")
refsize=$(cat "$genome" | wc -m | awk '{print $1}')

#Fromat $n.keepTheseZeroCovPositions to VCF
cat "${variant}"/"${n}".keepTheseZeroCovPositions \
    | sed 's/%/ /' \
    | awk 'BEGIN{OFS="\t"}{print $1, $2, ".", ".", ".", ".", ".", ".", "GT", "./."}' \
    > "${variant}"/"${n}".vcfFormated

cat "${variant}"/"${n}".body "${variant}"/"${n}".vcfFormated \
    | awk 'BEGIN{OFS="\t"}{if ($4 == ".") print $1, $2, $3, "N", $5, $6, $7, $8, $9, $10; else print $0}' \
    > "${variant}"/"${n}".SNPsMapzeroNoHeader.vcf

cat "${variant}"/"${n}".header "${variant}"/"${n}".SNPsMapzeroNoHeader.vcf > "${variant}"/"${n}".unsortSNPsZeroCoverage.vcf

java "$memJava" -jar "$igvtools" sort \
    "${variant}"/"${n}".unsortSNPsZeroCoverage.vcf \
    "${variant}"/"${n}".SNPsZeroCoverage.vcf

java "$memJava" -jar "$igvtools" index \
    "${variant}"/"${n}".SNPsZeroCoverage.vcf



# Emit all sites to VCF, not just the SNPs and indels.
# This allows making a UnifiedGenotyper VCF similar to what was used before using the Haplotypecaller.
java "$memJava" -jar "$gatk" -T UnifiedGenotyper \
    -R "$genome" \
    -out_mode EMIT_ALL_SITES \
    -I "${realigned}"/"${n}"_realigned_recalibrated.bam \
    -o "${variant}"/"${n}".allsites.vcf \
    -nt "$cpu"

wait

# This removes all positions same as the reference. 
# These positions are found by removing rows were column (field) 8 begins with AN=2.
# This decreases the size of the VCF considerably.
# The final VCF contains all SNP, indel or zero mapped/coverage positions
cat "${variant}"/"${n}".allsites.vcf \
    | awk ' $0 ~ /#/ || $8 !~ /^AN=2;/ {print $0}' \
    > "${variant}"/"${n}".diffsites.vcf

java "$memJava" -jar "$igvtools" index \
    "${variant}"/"${n}".diffsites.vcf

wait

#clean up
rm igv.log
# find "${variant}" -type f | grep -v "SNPsZeroCoverage" # | xargs rm

#################
#               #
#    Metrics    #
#               #
#################


#Quality Score Distribution
echo "***Quality Score Distribution"
java "$memJava" -jar "$picard" QualityScoreDistribution \
    REFERENCE_SEQUENCE="$genome" \
    INPUT="${realigned}"/"${n}"_realigned_recalibrated.bam \
    CHART_OUTPUT="${variant}"/"${n}".QualityScorceDistribution.pdf \
    OUTPUT="${variant}"/"${n}".QualityScoreDistribution \
    ASSUME_SORTED=true

wait

#Mean Quality by Cycle
echo "***Mean Quality by Cycle"
java "$memJava" -jar "$picard" CollectMultipleMetrics \
    REFERENCE_SEQUENCE="$genome" \
    INPUT="${realigned}"/"${n}"_realigned_recalibrated.bam \
    OUTPUT="${variant}"/"${n}".Quality_by_cycle \
    PROGRAM=MeanQualityByCycle \
    ASSUME_SORTED=true

wait

#Collect Alignment Summary Metrics
echo "***Collect Alignment Summary Metrics"
java "$memJava" -jar "$picard" CollectAlignmentSummaryMetrics \
    REFERENCE_SEQUENCE="$genome" \
    INPUT="${realigned}"/"${n}"_realigned_recalibrated.bam \
    OUTPUT="${variant}"/"${n}".AlignmentMetrics \
    ASSUME_SORTED=true

wait

#Collect GC Bias Error
echo "***Collect GC Bias Error"
java "$memJava" -jar "$picard" CollectGcBiasMetrics \
    REFERENCE_SEQUENCE="$genome" \
    INPUT="${realigned}"/"${n}"_realigned_recalibrated.bam \
    OUTPUT="${variant}"/"${n}".CollectGcBiasMetrics \
    CHART_OUTPUT="${variant}"/"${n}".GC.pdf \
    SUMMARY_OUTPUT="${variant}"/"${n}".GC.summary.txt \
    ASSUME_SORTED=true

wait

#Collect Insert Size Metrics
echo "***Collect Insert Size Metrics"
java "$memJava" -jar "$picard" CollectInsertSizeMetrics \
    REFERENCE_SEQUENCE="$genome" \
    INPUT="${realigned}"/"${n}"_realigned_recalibrated.bam \
    HISTOGRAM_FILE="${variant}"/"${n}".InsertSize.pdf \
    OUTPUT="${variant}"/"${n}".CollectInsertSizeMetrics \
    ASSUME_SORTED=true

wait

#Move to qualityvalues subfolder
qual=""${sampleDir}"/qualityvalues"
mkdir "$qual"
mv "${variant}"/*.pdf "$qual"


#################
#               #
#     Stats     #
#               #
#################


#create a new stats file
echo "fastq.gz file sizes:" > "${sampleDir}"/"${n}".stats.txt
for f in $(find -L "$sampleDir" -type f | grep -F ".fastq.gz" | grep -F ""${n}"_R"); do
    #stat -> file size in bytes
    #awk -> convert ot MB
    # rest -> remove decimal (not rounded)
    stat -Lc %s "$f" \
        | awk '{ foo = $1 / 1024 / 1024 ; print foo ".MB" }' \
        | cut -d "." -f 1,3 | tr -d "." \
        >> "${sampleDir}"/"${n}".stats.txt
done

wait

# echo "Unmapped fastq file sizes:" | tee -a "${sampleDir}"/"${n}".stats.txt
# du -sh ../unmappedReads/*.gz  | sed 's%../unmappedReads%%' | tee -a "${sampleDir}"/"${n}".stats.txt

# echo "Unmapped contig count:" | tee -a "${sampleDir}"/"${n}".stats.txt
# grep -c ">" ./spades_output/scaffolds.fasta | tee -a "${sampleDir}"/"${n}".stats.txt
# echo "" | tee -a "${sampleDir}"/"${n}".stats.txt

#TOTAL_READS
cat "${variant}"/"${n}".AlignmentMetrics \
    | sed -n 7,8p \
    | awk '{print $2}' \
    >> "${sampleDir}"/"${n}".stats.txt

#PF_READS
cat "${variant}"/"${n}".AlignmentMetrics \
    | sed -n 7,8p \
    | awk '{print $3}' \
    >> "${sampleDir}"/"${n}".stats.txt

#PF_ALIGNED_BASES
cat "${variant}"/"${n}".AlignmentMetrics \
    | sed -n 7,8p \
    | awk '{print $8}' \
    >> "${sampleDir}"/"${n}".stats.txt

#Average depth of coverage
aveCoverage=$(cat "${realigned}"/"${n}"_coverage.txt \
    | awk 'NR > 1' \
    | awk '{sum+=$4} END { print sum/NR"X"}')

wait

echo -e "\nAverage depth of coverage: "$aveCoverage"" | tee -a "${sampleDir}"/"${n}".stats.txt

#genome coverage
percGenomeMissing=$(awk -v x="$zeroposition" -v y="$refsize" 'BEGIN { print(x/y)*100}')
percGenomeCoverage="$(echo "100 - $percGenomeMissing" | bc)"
echo -e "\nPercent of reference with coverage: "$percGenomeCoverage"%" | tee -a "${sampleDir}"/"${n}".stats.txt

#cat ${n}.stats.txt | grep -v "Failed" | grep -v "Duplicates" | grep -v "Proper-pairs" | tee -a "${sampleDir}"/"${n}".stats.txt

#Mean insert size
echo -e "\nMean_Insert_Size\tStandard_Deviation:" >> "${sampleDir}"/"${n}".stats.txt
cat "${variant}"/"${n}".CollectInsertSizeMetrics \
    | awk 'BEGIN {OFS="\t"} { print $5,$6 }' \
    | awk 'FNR == 8 {print $0}' \
    >> "${sampleDir}"/"${n}".stats.txt

#Mean read lengh
echo -e "\nMean_Read_Length:" | tee -a "${sampleDir}"/"${n}".stats.txt
cat "${variant}"/"${n}".AlignmentMetrics \
    | awk 'BEGIN {OFS="\t"} { print $16 }' \
    | awk 'FNR == 10 {print $0}' \
    >> "${sampleDir}"/"${n}".stats.txt

#Add SNP call numbers to stats.txt file
echo -e "\nSNP and zero coverage positions:" | tee -a "${sampleDir}"/"${n}".stats.txt
cat "${variant}"/"${n}".SNPsZeroCoverage.vcf \
    | egrep -v "#" \
    | grep -c ".*" \
    | tee -a "${sampleDir}"/"${n}".stats.txt

#AC=2
echo -e "\nSNPs of AC2 and QUAL > 300:" | tee -a "${sampleDir}"/"${n}".stats.txt
cat "${variant}"/"${n}".SNPsZeroCoverage.vcf \
    | egrep -v "#" \
    | egrep "AC=2" \
    | awk '$6 > 300' \
    | grep -c ".*" \
    | tee -a "${sampleDir}"/"${n}".stats.txt

#Read counts
readcount=$(cat "${variant}"/"${n}".AlignmentMetrics \
    | sed -n 8p \
    | awk '{print $3}')

#Put in the coverageReport
echo -e ""$n"\t"$ID"\t"$readcount"\t"$aveCoverage"\t"$percGenomeCoverage"%" | tee -a "${reports}"/coverageReport.txt "${reports}"/dailyReport.txt


# if [ -f "${sampleDir}"/*identifier_out.txt ];then
#     cat "${sampleDir}"/*out1.txt ../*out2.txt > ../${n}-identification.txt
#     rm "${sampleDir}"/*identifier_out*
# fi


#cleanup
find "${realigned}" -type f | grep -v "recalibrated" | xargs rm
# mv ${n}.Metrics_summary.txt "$qual"
# mv ${n}.stats.txt "$qual"
# rm ${n}.Quality_by_cycle.insert_size_metrics
# rm ${n}.AlignmentMetrics
# mv "${startingdir}"/fastq ${startingdir}/spoligo
# rm "${startingdir}"/spoligo/*fastq
# rm -r "${startingdir}"/temp
# ln "${qual}"/"${n}".stats.txt ./stats-"${n}".txt


#Make dailyStats.txt for each stats.txt made for each isolate.
echo -e "\n\n\n" >> "${reports}"/dailyStats.txt
echo "ADD_MARKER" >> "${reports}"/dailyStats.txt
echo "" >> "${reports}"/dailyStats.txt
echo "<------- "$n" "$ID" ------->" >> "${reports}"/dailyStats.txt
cat "${sampleDir}"/"${n}".stats.txt >> "${reports}"/dailyStats.txt

echo "**************************** END "$n" ****************************"
