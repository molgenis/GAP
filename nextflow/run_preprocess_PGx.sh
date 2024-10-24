#!/bin/bash

set -e
set -u

function showHelp() {
	#
	# Display commandline help on STDOUT.
	#
	cat <<EOH
===============================================================================================================
Script to preprocess GAP data to PGX input
Usage:
	$(basename "${0}") OPTIONS
Options:
	-h   Show this help.

    Required:
	-s   samplesheet [path to samplesheet (including GDIO controls)]
	-f   gtcDir [path to gtcDir]
	
    Optional:
	-n   nextflow [path to nextflow workflow]
	-c   recalculateControls [default=no] means that all samples in the 'pipeline' column tagged as research will be included in a later step of the preprocessing]
	     [yes = all samples from samplesheet will be preprocessed]
===============================================================================================================
EOH
	trap - EXIT
	exit 0
}


while getopts "s:c:g:n:h" opt;
do
	case $opt in h)showHelp;; s)samplesheet="${OPTARG}";; g)gtcDir="${OPTARG}";; n)nextflow="${OPTARG}" ;; c)calculateControls="${OPTARG}"
esac
done

if [[ -z "${samplesheet:-}" &&  -z "${gtcDir:-}" ]]
then
	showHelp
	echo "samplesheet and gtcDir not defined" 
fi
if [[ -z "${nextflow:-}" ]]; then nextflow="/groups/umcg-pgx/tmp07/umcg-gvdvries/GAP/nextflow/main.nf" ; echo "nextflow workflow will be set to: /groups/umcg-pgx/tmp07/umcg-gvdvries/GAP/nextflow/main.nf"; fi
if [[ -z "${calculateControls:-}" ]] ; then calculateControls="no"; echo "controls will not be recalculated" ;fi


if [[ "${calculateControls}" == "yes" ]]
then
	newSamplesheet="${samplesheet}"
else

	declare -a _sampleSheetColumnNames=()
	declare -A _sampleSheetColumnOffsets=()

	workDir=$(dirname "${samplesheet}")
	newSamplesheet="${workDir}/diagnosticSamplesOnly.csv"
	IFS="," read -r -a _sampleSheetColumnNames <<< "$(head -1 "${samplesheet}")"

	for (( _offset = 0 ; _offset < ${#_sampleSheetColumnNames[@]} ; _offset++ ))
	do
		_sampleSheetColumnOffsets["${_sampleSheetColumnNames[${_offset}]}"]="${_offset}"
	done

	if [[ -n "${_sampleSheetColumnOffsets["pipeline"]+isset}" ]] 
	then
		echo "column [pipeline] is found in the samplesheet"
		_pipelineFieldIndex=$((${_sampleSheetColumnOffsets["pipeline"]} + 1))
		head -1 "${samplesheet}" > "${newSamplesheet}"
		awk -v pipeline="${_pipelineFieldIndex}" 'BEGIN {FS=","}{if(NR>1 &&$pipeline != "research"){print $0}}' "${samplesheet}" >> "${newSamplesheet}"
	fi

	echo "New samplesheet with only non-research samples is this: ${newSamplesheet}"
fi
echo "${newSamplesheet}" 
ml nextflow

nextflow run -profile slurm --samplesheet "${newSamplesheet}" --gtcDir "${gtcDir}" "${nextflow}" --recalculate "${calculateControls}" -resume
