#MOLGENIS walltime=02:00:00 mem=4gb
#list SentrixBarcode_A,SentrixPosition_A
#string intermediateDir
#string resultDir
#string computeVersion
#string Project
#string projectJobsDir
#string projectRawTmpDataDir
#string intermediateDir
#string genScripts
#string pipeline
#string runID
#string logsDir
#string perlVersion
#string group
#string gapVersion
#string workDir
#string workflowpath
#string host

umask 0007

module load ${computeVersion}
module list

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

#Create ProjectDirs
mkdir -p -m 2770 "${intermediateDir}"
mkdir -p -m 2770 "${resultDir}"
mkdir -p -m 2770 "${projectJobsDir}"
mkdir -p -m 2770 "${projectRawTmpDataDir}"

rocketPoint=$(pwd)

declare -a arrayUniqueMissingSampleNames

for sample in "${SentrixBarcode_A[@]}"
do
	array_contains arrayUniqueMissingSampleNames "${sample}" || arrayUniqueMissingSampleNames+=("${sample}")
done
	
declare -a arrayMissingSampleNames
cd "${projectRawTmpDataDir}"
for i in "${arrayUniqueMissingSampleNames[@]}"
do
	if [[ -f "../../../../../../rawdata/array/GTC/${i}/missingIDATs.txt" ]]
	then
		arrayRejected=()
		while read line
		do
			missingPosition=$(echo "${line}" | awk 'BEGIN {FS=":"}{print $2}')
			missingSampleName=$(echo "${line}" | awk 'BEGIN {FS=":"}{print $1}')
			arrayMissingPosition+=("${missingPosition}")
			arrayMissingSampleNames+=("${missingSampleName}")
		done<"../../../../../../rawdata/array/GTC/${i}/missingIDATs.txt"
	fi
done

max_index=${#SentrixBarcode_A[@]}-1


if [ "${pipeline}" == 'diagnostics' ] 
then
	for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
	do
		array_contains_missingSamples arrayMissingPosition "${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}"
		if [ "${missing}" == "false" ]
		then
			ln -sf "../../../../../../rawdata/array/GTC/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc" \
				"${projectRawTmpDataDir}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc"

			ln -sf "../../../../../../rawdata/array/GTC/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc.md5" \
			"${projectRawTmpDataDir}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc.md5"
		else
			echo -e "\n SAMPLE IS MISSING: ${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}"
		fi
	done
else
	for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
	do
		mkdir -p "${SentrixBarcode_A[samplenumber]}"
		ln -sf "../../../../../../rawdata/array/GTC/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc" \
			"${projectRawTmpDataDir}/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc"

		ln -sf "../../../../../../rawdata/array/GTC/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc.md5" \
		"${projectRawTmpDataDir}/${SentrixBarcode_A[samplenumber]}/${SentrixBarcode_A[samplenumber]}_${SentrixPosition_A[samplenumber]}.gtc.md5"
	done
fi

sampleSheetCsv="${genScripts}/${Project}.csv"
perl -pi -e 's/\r(?!\n)//g' "${sampleSheetCsv}"

cd "${rocketPoint}"

if [[ "${#arrayMissingSampleNames[@]:-0}" -ne "0" ]]
then
	teller=1
	size=${#arrayMissingSampleNames[@]}
	for m in "${arrayMissingSampleNames[@]}"
	do
		if [[ "${teller}" -lt "${size}" ]]
		then
			missingIDATsGrepCommand+="${m}|"
		elif [ "${teller}" == ${size} ]
		then
			echo "last line"
			missingIDATsGrepCommand+="${m}"
		fi
		teller=$((teller+1))
	done
		echo "missingIDATsGrepCommand: ${missingIDATsGrepCommand}"
	if grep -E "${missingIDATsGrepCommand}" "${sampleSheetCsv}" > "${Project}.removedSamples.csv"
	then
		grep -E -v "${missingIDATsGrepCommand}" "${sampleSheetCsv}" > "${intermediateDir}/${Project}.filteredSamplesheet.csv"
		sampleSheetCsv="${intermediateDir}/${Project}.filteredSamplesheet.csv"
	fi

fi

#Copying samplesheet to project jobs,results folder

cp "${sampleSheetCsv}" "${projectJobsDir}/${Project}.csv"
cp "${sampleSheetCsv}" "${resultDir}/${Project}.csv"

#
# Execute MOLGENIS/compute to create job scripts to analyse this project.
#
module load "${gapVersion}"
cd "${rocketPoint}"

perl "${EBROOTGAP}/scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/parameters_${host}.csv" > "${rocketPoint}/parameters_host_converted.csv"
perl "${EBROOTGAP}/scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/parameters_${group}.csv" > "${rocketPoint}/parameters_group_converted.csv"
perl "${EBROOTGAP}/scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/parameters_${pipeline}.csv" > "${rocketPoint}/parameters_converted.csv"

sh "${EBROOTMOLGENISMINCOMPUTE}/molgenis_compute.sh" \
-p "${genScripts}/parameters_converted.csv" \
-p "${genScripts}/parameters_group_converted.csv" \
-p "${genScripts}/parameters_host_converted.csv" \
-p "${sampleSheetCsv}" \
-p "${EBROOTGAP}/chromosomes.homo_sapiens.csv" \
-rundir "${projectJobsDir}" \
-w "${workflowpath}" \
--header "${EBROOTGAP}/templates/slurm/header.ftl" \
--submit "${EBROOTGAP}/templates/slurm/submit.ftl" \
--footer "${EBROOTGAP}/templates/slurm/footer.ftl" \
-b slurm \
-g \
-weave \
-runid "${runID}" \
-o "gapVersion=${gapVersion};\
runID=${runID}"


sampleSize=$(cat "${genScripts}/${Project}.csv" |  wc -l)

if [ -f "${Project}.removedSamples.csv" ]
then
	echo -e "\n################### THE FOLLOWING LINES ARE MISSING THE IDATS ###############\n"
	cat "${Project}.removedSamples.csv"
	cat "${Project}.removedSamples.csv" > "${logsDir}/${Project}/${Project}.pipeline.missingsamples"
fi


if [ "${pipeline}" == 'research' ] && [ "${sampleSize}" -gt 1000 ]
then
	echo "Samplesize is ${sampleSize}"
	ml "${perlVersion}"
	perl ${EBROOTGAP}/scripts/RemoveDuplicatesCompute.pl "${projectJobsDir}/"*"_mergeFinalReports_0.sh"
fi
