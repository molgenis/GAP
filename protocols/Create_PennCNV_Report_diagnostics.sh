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

set -e
set -u

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list

python "${EBROOTGAP}/scripts/Make_PennCNV_report_diagnostics.py" "${bpmFile}" "${projectRawTmpDataDir}" "${intermediateDir}"
