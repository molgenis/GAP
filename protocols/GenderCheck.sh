#!/bin/bash

#MOLGENIS walltime=01:59:00 mem=1gb ppn=1

#string callrateDir
#string Project
#string intermediateDir
#string genderCheckDir
#string logsDir


mkdir -p "${genderCheckDir}"

# remove header from the callrate file before comparing the genders
sed -e 1d "${callrateDir}/Callrates_${Project}.txt" > "${genderCheckDir}/callratedata_${Project}.txt"
rm -f "${logsDir}/${Project}/run01.pipeline.gendercheckfailed"

# reading callrate file and comparing the gender column
while read -r line
do
	genderIDColumn1=$(echo "${line}" | awk 'BEGIN { FS = "\t" } ; { print $1 }')
	genderColumn1=$(echo "${genderIDColumn1}" | awk 'BEGIN { FS = "_" } ; { print $4 }')
	genderColumn2=$(echo "${line}" | awk 'BEGIN { FS = "\t" } ; { print $3 }')

	if [[ "${genderColumn1}" != "${genderColumn2}" ]]
	then
		printf 'Gender check for sample %s failed: gender from adlas is %s, but callculated gender from array is %s.\r\n' \
			"${genderIDColumn1}" \
			"${genderColumn1}" \
			"${genderColumn2}" \
			>> "${logsDir}/${Project}/run01.pipeline.gendercheckfailed"
	fi

done <"${genderCheckDir}/callratedata_${Project}.txt"