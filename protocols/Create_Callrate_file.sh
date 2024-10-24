#MOLGENIS walltime=05:59:00 mem=10gb ppn=1

#string Project
#string pythonVersion
#string beadArrayVersion
#string intermediateDir
#string bpmFile
#string projectRawTmpDataDir
#list Sample_ID
#list SentrixBarcode_A
#list SentrixPosition_A
#string logsDir
#string diagnosticOutputFolder
#string resultDir
#string pipeline
#string CallrateDir
#string gapVersion

set -e
set -u
set -o pipefail

#Function to check if array contains value
array_contains () {
	local array="${1}[@]"
	local seeking="${2}"
	local in='no'
	for element in "${!array-}"; do
		if [[ "${element}" == "${seeking}" ]]; then
			in='yes'
			break
		fi
	done
	echo "${in}"
}

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list

mkdir -p "${CallrateDir}/"

makeTmpDir "${CallrateDir}/"
tmpCallrateDir="${MC_tmpFile}"

mkdir -p "${diagnosticOutputFolder}/${Project}"

INPUTARRAYS=()
for array in "${SentrixBarcode_A}_${SentrixPosition_A}"[@]
do
	element_exists="$(set -e; array_contains INPUTARRAYS "${array}")"
	if [[ "${element_exists}" == 'no' ]]; then
		# If file does not exist in array add it.
		INPUTARRAYS+=("${array}")
	fi
done

if [[ "${pipeline}" == 'research' ]]
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
	python "${EBROOTGAP}/scripts/Make_Callrate_Report.py" "${bpmFile}" "${projectRawTmpDataDir}" "${tmpCallrateDir}/callratedata_project.txt"
fi

#Create header for callrate report
echo -en "Sample ID\tCall Rate\tGender" > "${tmpCallrateDir}/callrate_header.txt"


#Add header to callrate report to create final results
#(cat ${intermediateDir}/callrate_header.txt; printf '\n'; cat ${intermediateDir}/callratedata_project.txt > "${intermediateDir}/Callrates_${Project}.txt")

printf "\n" | cat "${tmpCallrateDir}/callrate_header.txt" - "${tmpCallrateDir}/callratedata_project.txt" > "${tmpCallrateDir}/Callrates_${Project}.txt"

#add gender [replace M/F/U with Male/Female/Unknown ]

perl -pi -e 's|M|Male|g' "${tmpCallrateDir}/Callrates_${Project}.txt"
perl -pi -e 's|F|Female|g' "${tmpCallrateDir}/Callrates_${Project}.txt"
perl -pi -e 's|U|Unknown|g' "${tmpCallrateDir}/Callrates_${Project}.txt"

#Replace barcode with sampleid


max_index=${#Sample_ID[@]}-1

for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
do
echo "	perl -pi -e \"s|${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}|${Sample_ID[samplenumber]}|\" \"${tmpCallrateDir}/Callrates_${Project}.txt\""
	perl -pi -e "s|${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}|${Sample_ID[samplenumber]}|" "${tmpCallrateDir}/Callrates_${Project}.txt"
done

# Move results from tmp to intermediateDir


echo "mv ${tmpCallrateDir}/Callrates_${Project}.txt ${CallrateDir}/Callrates_${Project}.txt"
mv "${tmpCallrateDir}/Callrates_${Project}.txt" "${CallrateDir}/Callrates_${Project}.txt"


#Put results in resultsfolder

rsync -a "${CallrateDir}/Callrates_${Project}.txt" "${resultDir}"
rsync -a "${CallrateDir}/Callrates_${Project}.txt" "${diagnosticOutputFolder}/${Project}"
