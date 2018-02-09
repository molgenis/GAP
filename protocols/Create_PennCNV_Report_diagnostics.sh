#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string pythonVersion
#string beadArrayVersion
#string GSA
#string bpmFile
#string inputDir
#string intermediateDir
#string tmpTmpdir
#string tmpDir
#string workDir
#string tmpName

set -e
set -u

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module list

mkdir -p "${intermediateDir}"

python "${GSA}/Scripts/Make_PennCNV_report_diagnostics.py" "${bpmFile}" "${inputDir}" "${intermediateDir}"
