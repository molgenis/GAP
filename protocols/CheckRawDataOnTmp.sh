#!/bin/bash

#list SentrixBarcode_A
#list SentrixPosition_A
#string GTCtmpDataDir
#string intermediateDir
#string ateambotUser
#string Project
#string logsDir
#list Sample_ID

array_contains () {
	local array="$1[@]"
	local seeking="${2}"
	local in=1
	for element in "${!array-}"; do
		if [[ "${element}" == "${seeking}" ]]; then
			in=0
			break
		fi
	done
	return "${in}"
}

array_contains_missingSamples () {
	local array="$1[@]"
	local seeking="${2}"
	local in=1
	missing="false"
	for element in "${!array-}"; do
		if [[ "${element}" == "${seeking}" ]]; then
			in=0
		missing="true"
		continue
		fi
	done
}

allRawDataAvailable='true'

mkdir -p "${logsDir}/${Project}/"
max_index=${#Sample_ID[@]}-1
## Removing previously created data.requested file
rm -f "${logsDir}/${Project}/${Project}.data.requested"

declare -a arrayUniqueMissingSampleNames

for sample in "${SentrixBarcode_A[@]}"
do
	array_contains arrayUniqueMissingSampleNames "${sample}" || arrayUniqueMissingSampleNames+=("${sample}")
done
	
declare -a arrayMissingSampleNames

for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
do	
	dataProcessingStarted='false'
	TMPDATADIR="${GTCtmpDataDir}/${SentrixBarcode_A[samplenumber]}/"
	mkdir -vp "${TMPDATADIR}"
	
	for i in "${arrayUniqueMissingSampleNames[@]}"
	do
		if [[ -f "${TMPDATADIR}/missingIDATs.txt" ]]
		then
			arrayRejected=()
			while read line
			do
				missingPosition=$(echo "${line}" | awk 'BEGIN {FS=":"}{print $2}')
				missingSampleName=$(echo "${line}" | awk 'BEGIN {FS=":"}{print $1}')
				arrayMissingPosition+=("${missingPosition}")
				arrayMissingSampleNames+=("${missingSampleName}")
			done<"${TMPDATADIR}/missingIDATs.txt"
		fi
	done
	
	## Check if data copying already started, so only checking if file exists (when it is not started, it will create an ${logsDir}/${project}/${project}.data.requested)
	if [[ -f "${logsDir}/${Project}/${Project}.data.started" ]]
	then
		dataProcessingStarted='true'
	fi
	
	array_contains_missingSamples arrayMissingPosition "${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}"
	if [[ "${missing}" == "false" ]]
	then
		if [[ ! -f "${TMPDATADIR}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc" ]]
		then
			if [[ "${dataProcessingStarted}" == 'false' ]]
			then
				echo "${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc" >> "${logsDir}/${Project}/${Project}.data.requested"
			fi
			echo "${TMPDATADIR}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc is missing"
			allRawDataAvailable='false'
		else
			echo "${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc is available"
		fi
	else
		echo -e "Sample is not in original data: ${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc, skipped"
	fi
done

if [[ "${allRawDataAvailable}" == 'true' ]]
then
	echo "rawdata already available"
	if [[  -f "${logsDir}/${Project}/${Project}.data.started" ]]
	then
		mv "${logsDir}/${Project}/${Project}.data."{started,finished}
	else
		
		touch "${logsDir}/${Project}/${Project}.data.finished"
	fi
else
	echo "all Data is not yet available, exiting"
	rm -f "${logsDir}/${Project}/${Project}.data.finished"
	trap - EXIT
	exit 0
fi
