#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=30gb ppn=6

#string pythonVersion
#string beadArrayVersion
#string gapVersion
#string bpmFile
#string projectRawTmpDataDir
#string intermediateDir
#string SentrixBarcode_A
#string SentrixPosition_A
#string concordanceInputDir
#string Sample_ID
#string tmpName
#string logsDir
#string Project
#string GTCtoVCFVersion
#string fastaFile
#string GTCtmpDataDir
#string tmpTmpdir
#string HTSlibVersion

set -e
set -u


module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module load "${GTCtoVCFVersion}"
module load "${HTSlibVersion}"

# python "${EBROOTGAP}/scripts/gtc_final_report_diagnostics.py" "${bpmFile}" "${projectRawTmpDataDir}" "${intermediateDir}"

python "${EBROOTGTCTOVCF}/gtc_to_vcf.py" \
--gtc-paths "${GTCtmpDataDir}/${SentrixBarcode_A}/${SentrixBarcode_A}_${SentrixPosition_A}.gtc" \
--output-vcf-path "${intermediateDir}/" \
--manifest-file "${bpmFile}" \
--genome-fasta-file "${fastaFile}" \
--skip-indels \
--log-file "${tmpTmpdir}/${SentrixBarcode_A}_${SentrixPosition_A}_GTCtoVCF.log.txt"

log_ratios=$(awk '{if ($3 != "X" && $3 != "Y" && $3 != "XY" ) print $6}' "${intermediateDir}/${SentrixBarcode_A}_${SentrixPosition_A}.gtc.txt")
sd=$(echo "${log_ratios[@]}" | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}')

echo "${Sample_ID}": "${sd}"

if  [[ "${sd}" < 0.2 ]]
	then
	echo "move VCF to concordancedir when standard deviation <0.20 ."
	mv "${intermediateDir}/${SentrixBarcode_A}_${SentrixPosition_A}.vcf" "${intermediateDir}/${Sample_ID}.vcf"
	awk '{OFS="\t"}{if ($0 ~ "#CHROM" ){ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"'${Sample_ID}'"} else {print $0}}' "${intermediateDir}/${Sample_ID}.vcf" > "${intermediateDir}/${Sample_ID}.FINAL.vcf"
#	cp "${intermediateDir}/${Sample_ID}.FINAL.vcf" "${concordanceInputDir}/"
fi
