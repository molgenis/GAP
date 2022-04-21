#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=1

#string callrateDir
#string Project
#string intermediateDir
#string genderCheckDir
#string logsDir


mkdir -p "${genderCheckDir}"

# remove header from the callrate file before comparing the genders
sed -e 1d "${callrateDir}/Callrates_${Project}.txt" > "${genderCheckDir}/callratedata_${Project}.txt"

# reading callrate file and comparing the gender column
while read line
do
#       echo "${line}"
        genderIDColumn1=$(echo "${line}" | awk 'BEGIN { FS = "\t" } ; { print $1 }')
#       echo "genderIDColumn1 ${genderIDColumn1}..."
        genderColumn1=$(echo "${genderIDColumn1}" | awk 'BEGIN { FS = "_" } ; { print $4 }')
#       echo "genderColumn1:${genderColumn1}..."
        genderColumn2=$(echo "${line}" | awk 'BEGIN { FS = "\t" } ; { print $3 }')
#       echo "genderColumn2:${genderColumn2}..."

        if [[ "${genderColumn1}" == "${genderColumn2}" ]];
        then
                echo -e "gender check for sample ${genderIDColumn1} passed" >> "${logsDir}/${Project}/run01_pipeline.gendercheckpassed"
        else
                echo -e "gender check for sample ${genderIDColumn1} failed: gender from adlas:${genderColumn1}, gender callculated:${genderColumn2}" >> "${logsDir}/${Project}/run01.pipeline.gendercheckfailed"
        fi

#       genderID=$(awk 'BEGIN { FS = "_" } ; { print $3 }' "${genderIDColumn}")
#       echo "genderID ${genderID}"
#       genderCall=$(awk 'BEGIN { print $3}' "${line}")
#       echo "genderCall ${genderCall}"

done < "${genderCheckDir}/callratedata_${Project}.txt"