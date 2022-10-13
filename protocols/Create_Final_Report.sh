#!/bin/bash

#MOLGENIS walltime=02:00:00 mem=2gb ppn=1

#string pythonVersion
#string beadArrayVersion
#string gapVersion
#string bpmFile
#string projectRawTmpDataDir
#string intermediateDir
#string tmpTmpdir
#string tmpDir
#string workDir
#string tmpName
#string Project
#string logsDir
#string arrayFinalReport
#string samplesheet
#string SentrixBarcode_A

set -e
set -u

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list

makeTmpDir "${arrayFinalReport}"
tmpArrayFinalReport="${MC_tmpFile}"

python "${EBROOTGAP}/scripts/gtc_final_report.py" \
--manifest "${bpmFile}" \
--samplesheet "${samplesheet}" \
--gtc_directory "${projectRawTmpDataDir}/${SentrixBarcode_A}/" \
--output_file "${tmpArrayFinalReport}"


echo "mv ${tmpArrayFinalReport} ${arrayFinalReport}"
mv "${tmpArrayFinalReport}" "${arrayFinalReport}"
