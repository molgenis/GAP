#MOLGENIS walltime=02:00:00 mem=2gb ppn=1

#string intermediateDir
#string logsDir
#string Project
#list arrayFinalReport
#string finalReport

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

makeTmpDir "${finalReport}"
tmpFinalReport="${MC_tmpFile}"

INPUTREPORTS=()

for file in "${arrayFinalReport[@]}"
do
	element_exists="$(set -e; array_contains INPUTREPORTS "${file}")"
	if [[ "${element_exists}" == 'no' ]]; then
		# If file does not exist in array add it.
		INPUTREPORTS+=("${file}")
	fi
done

first="true"

for i in "${INPUTREPORTS[@]}"
do
	if [[ ${first} == "true" ]]
	then
		cat "${i}" > "${tmpFinalReport}"
		first='false'
		found="$(head -30 "${i}" | grep -n "\[Data\]" | grep -Eo '^[^:]+')"
		headerNumber=$(( ${found} + 2))
		echo "headerNumber:${headerNumber}"
	else
		tail -n+"${headerNumber}" "${i}" >> "${tmpFinalReport}"
	fi
done

echo "mv ${tmpFinalReport} ${finalReport}"
mv "${tmpFinalReport}" "${finalReport}"
