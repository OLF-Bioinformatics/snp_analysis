#!/bin/bash


######################
#                    #
#    User Defined    #
#                    #
######################


#where fastq reads are
export fastqPath="/media/3tb_hdd/data/Mycobaterium_bovis"

#Analysis root directory
export baseDir=""${HOME}"/analyses/mbovis_script1" #make variable global (for called scripts)

#script dependenties
export dependents=""${HOME}"/prog/snp_analysis/script_dependents"


#####################
#                   #
#   Data Stucture   #
#                   #
#####################


export reports=""${baseDir}"/reports"

[ -d "$baseDir" ] || mkdir -p "$baseDir"
[ -d "$reports" ] || mkdir -p "$reports"


##########################
#                        #
#   Initiating reports   #
#                        #
##########################


#Date
starttime=$(date +%s)
echo -e "Start time: "$(date)"" | tee "${reports}"/dailyTime.txt

#Run the whole thing
echo "Identifying species (TB, Brucella and paratuberculosis) and running processZips.sh..."

loopFiles2.sh
wait

#Runtime statistics
echo "End Time: $(date)" | tee -a "${reports}"/dailyTime.txt
endtime=$(date +%s)
runtime=$(($endtime - $starttime))
printf 'Total runtime: %dh:%dm:%ds\n' $(($runtime/3600)) $(($runtime%3600/60)) $(($runtime%60)) | tee -a "${reports}"/dailyTime.txt

echo "e-mailing files"

cat "${reports}"/dailyTime.txt > "${reports}"/email_processZips.txt
echo "" >> "${reports}"/email_processZips.txt
echo "ADD_MARKER" >> "${reports}"/email_processZips.txt
echo "" >> "${reports}"/dailyReport.txt
cat "${reports}"/dailyReport.txt >> "${reports}"/email_processZips.txt

cat "${reports}"/spoligoCheck.txt >> "${reports}"/email_processZips.txt
cat "${reports}"/mlstCheck.txt >> "${reports}"/email_processZips.txt

echo "ADD_MARKER" >> "${reports}"/email_processZips.txt

cat "${reports}"/dailyStats.txt >> "${reports}"/email_processZips.txt
echo "" >> "${reports}"/email_processZips.txt

cat "${reports}"/email_processZips.txt \
	| grep -v "*" \
	| grep -v "Stats for BAM file" \
    | sed 's/ADD_MARKER/******************************************/g' > "${reports}"/email_processZips2.txt

# rm "${reports}"/email_processZips.txt

if [ $1 == "me" ]; then
    email_list="marc-olivier.duceppe@inspection.gc.ca"
else
    email_list="marc-olivier.duceppe@inspection.gc.ca olga.andrievskaea@inspection.gc.ca susan.nadin-davis@inspection.gc.ca"
fi

cat "${reports}"/email_processZips2.txt | mail -s "WGS results" "$email_list"

echo $(date) >> "${reports}"/mlstCheck_all.txt
cat "${reports}"/mlstCheck.txt >> "${reports}"/mlstCheck_all.txt



 # Created by Tod Stuber on 11/09/12.

