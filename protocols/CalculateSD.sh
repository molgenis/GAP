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
#string tmpName
#string logsDir
#string Project
#string GTCtoVCF
#string fastaFile
#string GTCtmpDataDir
#string tmpTmpdir

set -e
set -u


module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module load "${GTCtoVCF}"

python "${EBROOTGAP}/scripts/gtc_final_report_diagnostics.py" "${bpmFile}" "${projectRawTmpDataDir}" "${intermediateDir}"

python "${EBROOTGTCTOVCF}"/gtc_to_vcf.py \
--gtc-paths "${GTCtmpDataDir}"/"${SentrixBarcode_A}"/ \
--output-vcf-path "${intermediateDir}"/ \
--manifest-file "${bpmFile}" \
--genome-fasta-file "${fastaFile}" \
--skip-indels \
--log-file "${tmpTmpdir}"/"${SentrixBarcode_A}"_GTCtoVCF.log.txt

#Replace barcode with sampleid


for ((i=0;i<${#Sample_ID[@]}-1;i++))
do
    barcodelist+=("${Sample_ID[${i}]}:${SentrixBarcode_A[${i}]}_${SentrixPosition_A[${i}]}")
done


for input_file in $(find "${intermediateDir}"/[0-9]*_R*.gtc.txt)
do

        log_ratios=$(awk '{if ($3 != "X" && $3 != "Y" && $3 != "XY" ) print $6}' "${input_file}")
        sd=$(echo "${log_ratios[@]}" | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}')

        echo "$(basename ${input_file}): ${sd}"

        for sample_barcode in ${barcodelist[@]}
        do
        filename=$(basename ${input_file})
        barcode_position=${filename%%.*}
        barcode=${barcode_position%_*}
        position=${barcode_position##*_}

                sample_id=$(echo "${sample_barcode}" | awk 'BEGIN {FS=":"}{print $1}')
                barcode_combined=$(echo "${sample_barcode}" | awk 'BEGIN {FS=":"} {print $2}')
                sentrix_barcode=$(echo "${barcode_combined}" | awk 'BEGIN {FS="_"}{print $1}')
                sentrix_position=$(echo "${barcode_combined}" | awk 'BEGIN {FS="_"}{print $2}')


                if  [[ "${sd}" < 0.2 && "${barcode}" == "${sentrix_barcode}" && "${position}" == "${sentrix_position}" ]]
                then
                    echo "${intermediateDir}/${input_file} ${concordanceInputDir}/${sample_id}.vcf"
                    mv "${intermediateDir}/${filename}.vcf" "${concordanceInputDir}/${sample_id}.vcf"
                    awk '{OFS="\t"}{if ($0 ~ "#CHROM" ){ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"'$sample_id'"} else {print $0}}' "${concordanceInputDir}/${sample_id}.vcf" > "${concordanceInputDir}/${sample_id}.FINAL.vcf"
                    bgzip -c "${concordanceInputDir}/${sample_id}.FINAL.vcf" > "${concordanceInputDir}/${sample_id}.FINAL.vcf.gz"
                    tabix -p vcf "${concordanceInputDir}/${sample_id}.FINAL.vcf.gz"
                fi
        done
done
