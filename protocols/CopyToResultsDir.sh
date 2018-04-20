#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6

#string intermediateDir
#string resultDir
#string Project
#string diagnosticOutputFolder

set -e
set -u


#Copying Diagnostics outputfiles to resultsDir


rsync -a "${intermediateDir}/${Project}_PennCNV.txt" "${resultDir}"
rsync -a "${intermediateDir}/Callrates_${Project}.txt" "${resultDir}"

#Copying files to Diagnostics output folder so DARWIN can further process the output files

mkdir -p "${diagnosticOutputFolder}/${Project}"

rsync -a "${intermediateDir}/${Project}_PennCNV.txt" "${diagnosticOutputFolder}/${Project}"
rsync -a "${intermediateDir}/Callrates_${Project}.txt" "${diagnosticOutputFolder}/${Project}"
