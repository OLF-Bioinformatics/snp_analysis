#!/bin/bash


######################
#                    #
#    User Defined    #
#                    #
######################


#script dependenties
export dependents=""${HOME}"/prog/snp_analysis/script_dependents"


#################
#               #
#    Options    #
#               #
#################


function displayHelp()
{
echo -e "\
Usage: email_loopFiles2.sh [flag] -b analysisFolder -v vcfSourceFolder

Mandatory flags:

    -b         base directory. Folder in which the analysis will take place
    -f         Fastq files root folder. Will look recursively for \".fastq.gz\" files

Optional flags:

    -m         mail just \"M\"e
"
}

#Colored error message
BLUE='\033[1;34m'
NC='\033[0m' # No Color

baseDir=""
fastqPath=""
mflag=0

#Hardcoded paths for debug
# export fastqPath="/media/3tb_hdd/data/Mycobaterium_bovis/outbreak2016"
# export baseDir=""${HOME}"/analyses/mbovis_script1_2017" #make variable global (for called scripts)

options=':mb:f:h'

while getopts "$options" opt; do
    case "$opt" in
        m)
            mflag=1
            ;;
        b)
            export baseDir="$OPTARG"  # make variable global (for called scripts)
            ;;
        f)
            export fastqPath="$OPTARG"  # make variable global (for called scripts)
            ;;
        h)
            displayHelp
            exit 1
            ;;
        \?)
            printf ""${BLUE}"Invalid option: -"$OPTARG"\n\n"${NC}"" >&2
            # echo "Invalid option: -"$OPTARG"" >&2
            displayHelp
            exit 1
            ;;
        :)
            printf ""${BLUE}"Option -"$OPTARG" requires an argument.\n\n"${NC}"" >&2
            # echo "Option -"$OPTARG" requires an argument." >&2
            displayHelp
            exit 1
            ;;
    esac
done

shift $(($OPTIND - 1))

# Exit if option flags b or v are missing
if [[ -z "$baseDir" ]] || [[ -z "$fastqPath" ]]; then
    echo "Both \"-b\" and \"-f\" options are mandatory"
    displayHelp
    exit 1
fi


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

if [ "$mflag" -eq 1 ]; then
    email_list="marc-olivier.duceppe@inspection.gc.ca"
else
    email_list="marc-olivier.duceppe@inspection.gc.ca,olga.andrievskaea@inspection.gc.ca,susan.nadin-davis@inspection.gc.ca"
fi

cat "${reports}"/email_processZips2.txt | mail -s "WGS results" -t "$email_list"
# mail -s "WGS results" -t "$email_list" < "${reports}"/email_processZips2.txt

# if [ -e "${reports}"/mlstCheck.txt ]; then
#     echo $(date) >> "${reports}"/mlstCheck_all.txt
#     cat "${reports}"/mlstCheck.txt >> "${reports}"/mlstCheck_all.txt
# fi



 # Created by Tod Stuber on 11/09/12.

