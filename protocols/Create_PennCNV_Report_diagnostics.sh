#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

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
#list SentrixBarcode_A
#list SentrixPosition_A
#string PennCNV_reportDir

set -e
set -u

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list


mkdir -p "${PennCNV_reportDir}"

makeTmpDir "${PennCNV_reportDir}"
tmpPennCNV_reportDir="${MC_tmpFile}"


python "${EBROOTGAP}/scripts/Make_PennCNV_report_diagnostics.py" "${bpmFile}" "${projectRawTmpDataDir}" "${tmpPennCNV_reportDir}"

echo "mv ${tmpPennCNV_reportDir}/${SentrixBarcode_A}_${SentrixPosition_A}.gtc.txt ${PennCNV_reportDir}/"
mv "${tmpPennCNV_reportDir}/${SentrixBarcode_A}_${SentrixPosition_A}.gtc.txt" "${PennCNV_reportDir}/"
