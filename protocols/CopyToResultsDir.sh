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


#Touch file for DARWIN so they know pipeline is finished and can start proceeding step to put data in SNP Module...
touch "${diagnosticOutputFolder}/${Project}".finished

# Touch log file for GAP_Automated for starting copying project data to PRM

if [ -f "${logsDir}//${Project}/${runID}.pipeline.started" ]
then
	mv "${logsDir}/${project}/${runID}.pipeline".{started,finished}
else
	touch "${logsDir}/${project}/${runID}.pipeline.finished"
fi
rm -f "${logsDir}/${project}/${runID}.pipeline.failed"
echo "${logsDir}/${project}/${runID}.pipeline.finished is created"

if [ ! -d "${logsDir}/${project}/" ]
then
	mkdir -p "${logsDir}/${project}/"
fi

touch pipeline.finished