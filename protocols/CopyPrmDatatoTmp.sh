#list SentrixBarcode_A
#list SentrixPosition_A
#string GTCprmDataDir
#string GTCtmpDataDir
#string prmHost
#string ateambotUser
#string Project
#string logsDir

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

INPUTREPORTS=()

for file in "${SentrixBarcode_A[@]}"
do
	element_exists="$(set -e; array_contains INPUTREPORTS "${file}")"
	if [[ "${element_exists}" == 'no' ]]; then
		# If file does not exist in array add it.
		INPUTREPORTS+=("${file}")
	fi
done

for i in "${INPUTREPORTS[@]}"
do
	mkdir -vp "${GTCtmpDataDir}/${i}"
	GTC_DIR="${GTCprmDataDir}/${i}"

	if [[ "${prmHost}" == "localhost" ]]
	then
		rsync --verbose --recursive --links --no-perms --times --group --no-owner --devices --specials --checksum \
		"${GTC_DIR}" \
		"${GTCtmpDataDir}"

	else
		rsync --verbose --recursive --links --no-perms --times --group --no-owner --devices --specials --checksum \
		"${ateambotUser}@${prmHost}:${GTC_DIR}" \
		"${GTCtmpDataDir}"
	fi
done
