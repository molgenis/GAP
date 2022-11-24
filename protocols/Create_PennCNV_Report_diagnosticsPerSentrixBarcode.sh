#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string pythonVersion
#string beadArrayVersion
#string bpmFile
#string projectRawTmpDataDir
#string Project
#string SentrixBarcode_A
#list SentrixPosition_A
#string pennCNV_reportDir
#list Sample_ID
#string gapVersion
#string resultDir
#string logsDir
#string intermediateDir

module purge
module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list


mkdir -p "${pennCNV_reportDir}"
mkdir -p "${resultDir}/PennCNV_reports/"

makeTmpDir "${pennCNV_reportDir}"
tmpPennCNV_reportDir="${MC_tmpFile}"

## Make a list of all samples to be processed per SentrixBarcode

samplelist=()

max_index=${#Sample_ID[@]}-1
for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
do
	samplelist+=("${Sample_ID[samplenumber]}:${SentrixBarcode_A}_${SentrixPosition_A[samplenumber]}")
done

## Process all samples in the samplelist. A PennCNV report per sample is made, compatible with downstream diagnostics. 

for i in "${samplelist[@]}"
do
	python "${EBROOTGAP}/scripts/Make_PennCNV_report_diagnosticsPerSentrixBarcode.py" "${bpmFile}" "${projectRawTmpDataDir}" "${tmpPennCNV_reportDir}" "${i}"
	echo "processing $i"
	barcodeCombined=$(echo "${i}" | awk 'BEGIN {FS=":"}{print $1}')
	echo "${barcodeCombined}"
	echo "mv ${tmpPennCNV_reportDir}/${barcodeCombined}.txt ${resultDir}/PennCNV_reports/"
	mv "${tmpPennCNV_reportDir}/${barcodeCombined}.txt" "${resultDir}/PennCNV_reports/"
done

