#!/bin/bash

inputDir="${1}"
finalReport="${2}"


set -e
set -u

inputReports=("${inputDir}/"*)
tmpFinalReport="${finalReport}.tmp"

first="true"

for i in "${inputReports[@]}"
do
	if [[ ${first} == "true" ]]
	then
		cat ${i} > ${tmpFinalReport}
		first='false'
		headerNumber=$(( $(head -30 "${i}" |grep -n "\[Data\]" | grep -Eo '^[^:]+')+2))
		echo "headerNumber:${headerNumber}"
	else
		cat ${i} | tail -n+${headerNumber} >> ${tmpFinalReport}
	fi
done

echo "mv ${tmpFinalReport} ${finalReport}"
mv "${tmpFinalReport}" "${finalReport}"
