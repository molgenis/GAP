#MOLGENIS walltime=02:00:00 mem=2gb ppn=1

#string intermediateDir
#string logsDir
#string Project
#list arrayFinalReport
#string finalReport

set -e
set -u

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

makeTmpDir "${finalReport}"
tmpFinalReport="${MC_tmpFile}"

INPUTREPORTS=()

for file in "${arrayFinalReport[@]}"
do
	array_contains INPUTREPORTS "${file}" || INPUTREPORTS+=("$file")    # If bamFile does not exist in array add it
done

first="true"

for i in "${INPUTREPORTS[@]}"
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
