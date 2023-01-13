#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string Project
#list Sample_ID
#list SentrixBarcode_A
#list SentrixPosition_A
#string resultDir
#string callrateDir
#string gapVersion
#string logsDir
#string intermediateDir
#string diagnosticOutputFolder

set -e
set -u
set -o pipefail

#Function to check if array contains value
array_contains () {
	local array="$1[@]"
	local seeking=$2
	local in=1
	for element in "${!array-}"; do
		if [[ "${element}" == "${seeking}" ]]; then
			in=0
			break
		fi
	done
	return "${in}"
}

module load "${gapVersion}"
module list

INPUTARRAYS=()

for array in "${SentrixBarcode_A[@]}"
do
	array_contains INPUTARRAYS "${array}" || INPUTARRAYS+=("${array}")    # Make a list of unique SentrixBarcode_A per project.
done

## Merge all Callrate files from different SentrixBarcode_A to one project Callrate file.
echo -e "Sample ID\tCall Rate\tGender" > "${callrateDir}/Callrates_${Project}.txt"

for i in "${INPUTARRAYS[@]}"
do
	echo "${callrateDir}/Callrates_${i}.txt"
	awk 'FNR>1' "${callrateDir}/Callrates_${i}.txt" >> "${callrateDir}/Callrates_${Project}.txt"
done


#Put results in resultsfolder
rsync -a "${callrateDir}/Callrates_${Project}.txt" "${resultDir}"

#rsync to DiagnosticOutput folder, can be removed when data is stored on shared isilon storage.
rsync -a "${callrateDir}/Callrates_${Project}.txt" "${diagnosticOutputFolder}/${Project}/"