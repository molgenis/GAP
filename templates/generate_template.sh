#!/bin/bash

module load Molgenis-Compute/v17.08.1-Java-1.8.0_74
module load GAP/v2.2.1
module list

host=$(hostname -s)
environmentParameters="parameters_${host}"

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to copy (sync) data from a succesfully finished analysis project from tmp to prm storage.
Usage:
	$(basename $0) OPTIONS
Options:
	-h   Show this help.
	-a   sampleType (default=GAP)
	-p   project
	-g   group (default=basename of ../../../ )
	-f   filePrefix (default=basename of this directory)
	-r   runID (default=run01)
	-t   tmpDirectory (default=basename of ../../ )
	-x   excludeGTCFiles
	-w   workdir (default=/groups/\${group}/\${tmpDirectory})

===============================================================================================================
EOH
	trap - EXIT
	exit 0
}

while getopts "t:g:w:f:r:h:x:" opt;
do
	case $opt in h)showHelp;; t)tmpDirectory="${OPTARG}";; g)group="${OPTARG}";; w)workDir="${OPTARG}";; f)filePrefix="${OPTARG}";; p)project="${OPTARG}";; r)runID="${OPTARG}";; x)excludeGTCsFile="${OPTARG}";;
	esac
done

if [[ -z "${tmpDirectory:-}" ]]; then tmpDirectory=$(basename $(cd ../../ && pwd )) ; fi ; echo "tmpDirectory=${tmpDirectory}"
if [[ -z "${group:-}" ]]; then group=$(basename $(cd ../../../ && pwd )) ; fi ; echo "group=${group}"
if [[ -z "${workDir:-}" ]]; then workDir="/groups/${group}/${tmpDirectory}" ; fi ; echo "workDir=${workDir}"
if [[ -z "${filePrefix:-}" ]]; then filePrefix=$(basename $(pwd )) ; fi ; echo "filePrefix=${filePrefix}"
if [[ -z "${runID:-}" ]]; then runID="run01" ; fi ; echo "runID=${runID}"
if [[ -z  "${excludeGTCsFile}" ]];then excludeGTCsFile="FALSE" ; fi ; echo "excludeGTCsFile=${excludeGTCsFile}"
genScripts="${workDir}/generatedscripts/${filePrefix}/"
samplesheet="${genScripts}/${filePrefix}.csv" ; mac2unix "${samplesheet}"

### Which pipeline to run
samplesHeetColumnNames=()
sampleSheetColumnOffsets=()
IFS="${SAMPLESHEET_SEP}" sampleSheetColumnNames=($(head -1 "${_samplesheet}"))
for (( _offset = 0 ; _offset < ${#sampleSheetColumnNames[@]:-0} ; _offset++ ))
do
	_sampleSheetColumnOffsets["${sampleSheetColumnNames[${_offset}]}"]="${_offset}"
done
if [[ ! -z "${sampleSheetColumnOffsets['pipeline']+isset}" ]]; then
	pipelineFieldIndex=$((${sampleSheetColumnOffsets['pipeline']} + 1))
	IFS=$'\n' pipeline=($(tail -n +2 "${sampleSheet}" | cut -d "${SAMPLESHEET_SEP}" -f ${pipelineFieldIndex} | head -1))
else
	pipeline="diagnostics"
fi

host=$(hostname -s)
echo "${host}"

projectDir="${workDir}/projects/${filePrefix}/${runID}/jobs/"

mkdir -p -m 2770 "${workDir}/projects/"
mkdir -p -m 2770 "${workDir}/projects/${filePrefix}/"
mkdir -p -m 2770 "${workDir}/projects/${filePrefix}/${runID}/"
mkdir -p -m 2770 "${workDir}/projects/${filePrefix}/${runID}/jobs/"

#samplesheet="${genScripts}/${filePrefix}.csv" ; mac2unix "${samplesheet}"

perl "${EBROOTGAP}/scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/parameters_${host}.csv" > "${genScripts}/parameters_host_converted.csv"
perl "${EBROOTGAP}/scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/${pipeline}_parameters.csv" > "${genScripts}/parameters_converted.csv"

sh "${EBROOTMOLGENISMINCOMPUTE}/molgenis_compute.sh" \
-p "${genScripts}/parameters_converted.csv" \
-p "${genScripts}/parameters_host_converted.csv" \
-p "${samplesheet}" \
-w "${EBROOTGAP}/Prepare_${pipeline}_workflow.csv" \
-weave \
--generate \
-rundir "${genScripts}/scripts" \
--runid "${runID}" \
-o "outputdir=scripts/jobs;\
mainParameters=${genScripts}/parameters_converted.csv;\
samplesheet=${samplesheet};\
Project=${filePrefix};\
pipeline=${pipeline};\
runID=${runID};\
excludeGTCsFile=${excludeGTCsFile:-};"
