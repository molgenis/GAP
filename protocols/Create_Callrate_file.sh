#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string Project
#string pythonVersion
#string beadArrayVersion
#string intermediateDir
#string bpmFile
#string projectRawTmpDataDir
#string gapVersion
#list Sample_ID
#list SentrixBarcode_A
#list SentrixPosition_A
#string logsDir

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list



set -e
set -u

python "${EBROOTGAP}/scripts/Make_Callrate_Report.py" "${bpmFile}" "${projectRawTmpDataDir}" "${intermediateDir}/callratedata_project.txt"

#Create header for callrate report
echo -en "Sample ID\tCall Rate\tGender" > "${intermediateDir}/callrate_header.txt"


#Add header to callrate report to create final results
#(cat ${intermediateDir}/callrate_header.txt; printf '\n'; cat ${intermediateDir}/callratedata_project.txt > "${intermediateDir}/Callrates_${Project}.txt")

printf "\n" | cat "${intermediateDir}/callrate_header.txt" - "${intermediateDir}/callratedata_project.txt" > "${intermediateDir}/Callrates_${Project}.txt"

#add gender [replace M/F/U with Male/Female/Unknown ]

perl -pi -e 's|M|Male|g' "${intermediateDir}/Callrates_${Project}.txt"
perl -pi -e 's|F|Female|g' "${intermediateDir}/Callrates_${Project}.txt"
perl -pi -e 's|U|Unknown|g' "${intermediateDir}/Callrates_${Project}.txt"

#Replace barcode with sampleid

barcodelist=()

n_elements=${Sample_ID[@]}
max_index=${#Sample_ID[@]}-1
for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
do
echo "	perl -pi -e \"s|${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}|${Sample_ID[samplenumber]}|\" \"${intermediateDir}/Callrates_${Project}.txt\""
	perl -pi -e "s|${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}|${Sample_ID[samplenumber]}|" "${intermediateDir}/Callrates_${Project}.txt"
done

