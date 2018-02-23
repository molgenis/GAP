#!/bin/bash

#MOLGENIS walltime=05:59:00 mem=10gb ppn=6



#string intermediateDir
#string resultDir
#string Project


set -e
set -u


#Copying Diagnostics outputfiles to resultsDir


mkdir -p ${resultDir}

rsync -a "${intermediateDir}/${Project}_PennCNV.txt" "${resultDir}"
rsync -a "${intermediateDir}/Callrates_${Project}.txt" "${resultDir}"
