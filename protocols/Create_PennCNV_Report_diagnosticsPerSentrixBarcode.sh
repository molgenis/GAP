#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string pythonVersion
#string beadArrayVersion
#string bpmFile
#string projectRawTmpDataDir
#string intermediateDir
#string tmpTmpdir
#string tmpDir
#string workDir
#string tmpName
#string Project
#string logsDir
#string SentrixBarcode_A
#list SentrixPosition_A
#string PennCNV_reportDir
#list Sample_ID
#string gapVersion

set -e
set -u

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list


mkdir -p "${PennCNV_reportDir}"

makeTmpDir "${PennCNV_reportDir}"
tmpPennCNV_reportDir="${MC_tmpFile}"

python "${EBROOTGAP}/scripts/Make_PennCNV_report_diagnostics.py" "${bpmFile}" "${projectRawTmpDataDir}" "${tmpPennCNV_reportDir}" "${SentrixBarcode_A}"

barcodelist=()

n_elements=${Sample_ID[@]}
max_index=${#Sample_ID[@]}-1
for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
do
    barcodelist+=("${Sample_ID[samplenumber]}:${SentrixBarcode_A}_${SentrixPosition_A[samplenumber]}")
done

for i in ${barcodelist[@]}
do
        echo "processing $i"
        barcodeCombined=$(echo ${i} | awk 'BEGIN {FS=":"}{print $2}')
        echo "${barcodeCombined}"
        echo "mv ${tmpPennCNV_reportDir}/${barcodeCombined}.gtc.txt ${PennCNV_reportDir}"
        mv "${tmpPennCNV_reportDir}/${barcodeCombined}.gtc.txt" "${PennCNV_reportDir}"
done

