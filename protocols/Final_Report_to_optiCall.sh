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
#string finalReport
#string samplesheet
#string optiCallDir

set -e
set -u

module load "${gapVersion}"
module list

mkdir -p "${optiCallDir}"

cd "${optiCallDir}"
bash ${EBROOTGAP}/scripts/GS_to_Opticall.sh -i "${finalReport}" -o "${optiCallDir}"
cd -

