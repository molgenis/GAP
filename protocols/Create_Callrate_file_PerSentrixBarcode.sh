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
#string callrateDir
#string gapVersion
#string logsDir
#string intermediateDir

#Function to check if array contains value

module load "${gapVersion}"
module list

set -e
set -u
set -o pipefail

mkdir -p "${callrateDir}/"

makeTmpDir "${callrateDir}/"
tmpCallrateDir="${MC_tmpFile}"

python "${EBROOTGAP}/scripts/Make_Callrate_Report_PerSentrixBarcode.py" "${bpmFile}" "${projectRawTmpDataDir}" "${SentrixBarcode_A}" "${tmpCallrateDir}/callratedata_project.txt"

#Create header for callrate report
echo -en "Sample ID\tCall Rate\tGender" > "${tmpCallrateDir}/callrate_header.txt"

printf "\n" | cat "${tmpCallrateDir}/callrate_header.txt" - "${tmpCallrateDir}/callratedata_project.txt" > "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"


#add gender [replace M/F/U with Male/Female/Unknown ]

perl -pi -e 's|M|Male|g' "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"
perl -pi -e 's|F|Female|g' "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"
perl -pi -e 's|U|Unknown|g' "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"

#Replace barcode with sampleid

max_index=${#Sample_ID[@]}-1
for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
do
	echo "	perl -pi -e \"s|${SentrixBarcode_A}_${SentrixPosition_A[samplenumber]}|${Sample_ID[samplenumber]}|\" \"${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt\""
	perl -pi -e "s|${SentrixBarcode_A}_${SentrixPosition_A[samplenumber]}|${Sample_ID[samplenumber]}|" "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt"
done

# Move results from tmp to intermediateDir

echo "mv ${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt ${callrateDir}/Callrates_${SentrixBarcode_A}.txt"
mv "${tmpCallrateDir}/Callrates_${SentrixBarcode_A}.txt" "${callrateDir}/Callrates_${SentrixBarcode_A}.txt"
