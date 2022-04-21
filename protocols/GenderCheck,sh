#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=1

#string CallrateDir
#string Project
#string intermediateDir
#string GenderCheckDir
#string logsDir


mkdir -p "${GenderCheckDir}"

# remove header from the callrate file before comparing the genders
sed -e 1d "${CallrateDir}/Callrates_${Project}.txt" > "${GenderCheckDir}/callratedata_${Project}.txt"

# reading callrate file and comparing the gender column
while read line
do
#       echo "${line}"
        GenderIDColum1=$(echo "${line}" | awk 'BEGIN { FS = "\t" } ; { print $1 }')
#       echo "GenderIDColum1 ${GenderIDColum1}..."
        GenderColum1=$(echo "${GenderIDColum1}" | awk 'BEGIN { FS = "_" } ; { print $4 }')
#       echo "GenderColum1:${GenderColum1}..."
        GenderColum2=$(echo "${line}" | awk 'BEGIN { FS = "\t" } ; { print $3 }')
#       echo "GenderColum2:${GenderColum2}..."

        if [[ "${GenderColum1}" == "${GenderColum2}" ]];
        then
                echo -e "gender check for sample ${GenderIDColum1} passed" >> "${logsDir}/${Project}/run01_pipeline.gendercheckpassed"
        else
                echo -e "gender check for sample ${GenderIDColum1} failed: gender from adlas:${GenderColum1}, gender callculated:${GenderColum2}" >> "${logsDir}/${Project}/run01.pipeline.gendercheckfailed"
        fi

#       GenderID=$(awk 'BEGIN { FS = "_" } ; { print $3 }' "${GenderIDColum}")
#       echo "GenderID ${GenderID}"
#       GenderCall=$(awk 'BEGIN { print $3}' "${line}")
#       echo "GenderCall ${GenderCall}"

done < "${GenderCheckDir}/callratedata_${Project}.txt"