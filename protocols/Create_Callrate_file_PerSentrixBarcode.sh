#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string Project
#string pythonVersion
#string beadArrayVersion
#string bpmFile
#string projectRawTmpDataDir
#list Sample_ID
#string SentrixBarcode_A
#list SentrixPosition_A
#string pipeline
#string CallrateDir
#string gapVersion
#string logsDir
#string intermediateDir

#Function to check if array contains value
array_contains () {
	local array="$1[@]"
	local seeking=$2
	local in=1
	for element in "${!array-}"; do
		if [[ "$element" == "$seeking" ]]; then
			in=0
			break
		fi
	done
	return $in
}

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list

set -e
set -u

mkdir -p "${CallrateDir}/"

makeTmpDir "${CallrateDir}/"
tmpCallrateDir="${MC_tmpFile}"


INPUTARRAYS=()

for array in "${SentrixBarcode_A}_${SentrixPosition_A[@]}"
do
	array_contains INPUTARRAYS "${array}" || INPUTARRAYS+=("$array")    # If GTCfile does not exist in array add it
done


if [ "${pipeline}" == 'research' ]
then
	for i in "${INPUTARRAYS[@]}"
	do
		python "${EBROOTGAP}/scripts/Make_Callrate_Report.py" "${bpmFile}" "${projectRawTmpDataDir}/${i}/" "${tmpCallrateDir}/${i}_callratedata_project.txt"
	done


	rm -f "${tmpCallrateDir}/callratedata_project.txt"
	for j in "${INPUTARRAYS[@]}"
	do
		echo "cat ${tmpCallrateDir}/${j}_callratedata_project.txt >> ${tmpCallrateDir}/callratedata_project.txt"
		cat "${tmpCallrateDir}/${j}_callratedata_project.txt" >> "${tmpCallrateDir}/callratedata_project.txt"
	done
else
	python "${EBROOTGAP}/scripts/Make_Callrate_Report_PerSentrixBarcode.py" "${bpmFile}" "${projectRawTmpDataDir}" "${SentrixBarcode_A}" "${tmpCallrateDir}/callratedata_project.txt"
fi

#Create header for callrate report
echo -en "Sample ID\tCall Rate\tGender" > "${tmpCallrateDir}/callrate_header.txt"

printf "\n" | cat "${tmpCallrateDir}/callrate_header.txt" - "${tmpCallrateDir}/callratedata_project.txt" > "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"


#add gender [replace M/F/U with Male/Female/Unknown ]

perl -pi -e 's|M|Male|g' "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"
perl -pi -e 's|F|Female|g' "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"
perl -pi -e 's|U|Unknown|g' "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"

#Replace barcode with sampleid

#barcodelist=()

#n_elements=${Sample_ID[@]}
max_index=${#Sample_ID[@]}-1
for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
do
	echo "	perl -pi -e \"s|${SentrixBarcode_A}_${SentrixPosition_A[samplenumber]}|${Sample_ID[samplenumber]}|\" \"${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt\""
	perl -pi -e "s|${SentrixBarcode_A}_${SentrixPosition_A[samplenumber]}|${Sample_ID[samplenumber]}|" "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"
done

# Move results from tmp to intermediateDir

echo "mv ${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt ${CallrateDir}/Callrates_${SentrixBarcode_A}.txt"
mv "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt" "${CallrateDir}/Callrates_${SentrixBarcode_A}.txt"



