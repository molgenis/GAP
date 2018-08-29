#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string pythonVersion
#string beadArrayVersion
#string gapVersion
#string bpmFile
#string projectRawTmpDataDir
#string intermediateDir
#list SentrixBarcode_A
#list SentrixPosition_A
#string concordanceInputDir
#list Sample_ID
#string logsDir
#string Project


set -e
set -u


module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"


python "${EBROOTGAP}/scripts/gtc_final_report_diagnostics.py" "${bpmFile}" "${projectRawTmpDataDir}" "${intermediateDir}"


for ((i=0;i<${#Sample_ID[@]}-1;i++))
do
	barcodelist+=("${Sample_ID[${i}]}:${SentrixBarcode_A[${i}]}_${SentrixPosition_A[${i}]}")
done


for input_file in "${intermediateDir}/"*".gtc.txt"
do
	barcode_position="${input_file%%.*}"
	barcode="${barcode_position%_*}"
	position="${barcode_position##_*}"

	log_ratios=$(awk '{if ($3 != "X" && $3 != "Y" && $3 != "XY" ) print $6}' "${input_file}")
	sd=$(echo "${log_ratios[@]}" | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}')

        echo "$(basename ${input_file}): ${sd}"


	for barcode_combined in ${barcodelist[@]}
	do
		sample_id=$(echo "${barcode_combined}" | awk 'BEGIN {FS=":"}{print $1}')
		sentrix_barcode=$(echo "${barcode_combined}" | awk 'BEGIN {FS="_"}{print $1}')
		sentrix_position=$(echo "${barcode_combined}" | awk 'BEGIN {FS="_"}{print $2}')

		if  [[ "${sd}" < 0.2 && "${barcode}" == "${sentrix_barcode}" && "${position}" == "${sentrix_position}" ]]
		then
			echo "mv ${intermediateDir}/concordance_${input_file} ${concordanceInputDir}"
                        mv "${intermediateDir}/concordance_${input_file}" "${concordanceInputDir}/concordance_${sample_id}.txt"
		fi

	done
done
