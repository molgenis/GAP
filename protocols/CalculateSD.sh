#MOLGENIS walltime=05:59:00 mem=8gb ppn=1

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
#string resultDir
#string calculateSDDir
#string pennCNV_reportDir

mkdir -p "${calculateSDDir}/"

makeTmpDir "${calculateSDDir}/"
tmpCalculateSDDir="${MC_tmpFile}"

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module load "${GTCtoVCFVersion}"
module load "${HTSlibVersion}"

#Convert gtc file to a VCF

python "${EBROOTGTCTOVCF}/gtc_to_vcf.py" \
--gtc-paths "${GTCtmpDataDir}/${SentrixBarcode_A}/${SentrixBarcode_A}_${SentrixPosition_A}.gtc" \
--output-vcf-path "${tmpCalculateSDDir}/" \
--manifest-file "${bpmFile}" \
--genome-fasta-file "${fastaFile}" \
--skip-indels \
--log-file "${tmpCalculateSDDir}/${SentrixBarcode_A}_${SentrixPosition_A}_GTCtoVCF.log.txt"

log_ratios=$(awk '{if (NR>10){if ($3 != "X" && $3 != "Y" && $3 != "XY" ){ print $5}}}' "${resultDir}/PennCNV_reports/${Sample_ID}.txt")s

sd=$(echo "${log_ratios[@]}" | awk '{sum+=$1; sumsq+=$1*$1}END{print sqrt(sumsq/NR - (sum/NR)**2)}')
echo "${sd}" > "${tmpCalculateSDDir}/${Sample_ID}_${SentrixBarcode_A}_${SentrixPosition_A}.FINAL.vcf.sd"

mv "${tmpCalculateSDDir}/${SentrixBarcode_A}_${SentrixPosition_A}.vcf" "${tmpCalculateSDDir}/${Sample_ID}_${SentrixBarcode_A}_${SentrixPosition_A}.vcf"
awk '{OFS="\t"}{if ($0 ~ "#CHROM" ){ print $1,$2,$3,$4,$5,$6,$7,$8,$9,"'"${Sample_ID}"_"${SentrixBarcode_A}"_"${SentrixPosition_A}"'"} else {print $0}}' "${tmpCalculateSDDir}/${Sample_ID}_${SentrixBarcode_A}_${SentrixPosition_A}.vcf" > "${tmpCalculateSDDir}/${Sample_ID}_${SentrixBarcode_A}_${SentrixPosition_A}.FINAL.vcf"

#Copy VCF to resultsdir
mkdir -p "${resultDir}/vcf/"

#Move vcf and sd values to intermediateDir
echo "moving ${tmpCalculateSDDir}/${Sample_ID}_${SentrixBarcode_A}_${SentrixPosition_A}.FINAL.vcf.sd ${resultDir}/vcf/"
mv "${tmpCalculateSDDir}/${Sample_ID}_${SentrixBarcode_A}_${SentrixPosition_A}.FINAL.vcf.sd" "${resultDir}/vcf/"
echo "moving ${tmpCalculateSDDir}/${Sample_ID}_${SentrixBarcode_A}_${SentrixPosition_A}.FINAL.vcf to ${resultDir}/vcf/"
mv "${tmpCalculateSDDir}/${Sample_ID}_${SentrixBarcode_A}_${SentrixPosition_A}.FINAL.vcf" "${resultDir}/vcf/"
