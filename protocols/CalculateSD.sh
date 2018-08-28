#!/bin/bash


#string pythonVersion
#string beadArrayVersion
#string gapVersion
#string bpmFile
#string projectRawTmpDataDir
#string intermediateDir
#string concordanceInputDir
#string Sample_ID
#string SentrixBarcode_A
#string SentrixPosition_A

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
	_barcode_position="${input_file%%.*}"
	_barcode="${barcode_position%_*}"
	_position="${barcode_position##_*}"

	_log_ratios=$(awk '{if ($3 != "X" && $3 != "Y" && $3 != "XY" ) print $6}' "${input_file}")
	_sd=$(echo "${_log_ratios[@]}" | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}')

        echo "$(basename ${input_file}): ${_sd}"


	for barcode_combined in ${barcodelist[@]}
	do
		_sample_id=$(echo "${barcode_combined}" | awk 'BEGIN {FS=":"}{print $1}')
		_sentrix_barcode=$(echo "${barcode_combined}" | awk 'BEGIN {FS="_"}{print $1}')
		_sentrix_position=$(echo "${barcode_combined}" | awk 'BEGIN {FS="_"}{print $2}')

		if  [[ "${sd}" < 0.2 && "${_barcode}" == "${_sentrix_barcode}" && "${_position}" == "${_sentrix_position}" ]]
		then
			echo "mv ${intermediateDir}/concordance_${input_file} ${concordanceInputDir}"
                        mv "${intermediateDir}/concordance_${input_file}" "${concordanceInputDir}/${_sample_id}"
		fi

	done
done
