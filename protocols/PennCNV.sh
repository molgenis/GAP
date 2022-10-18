#!/bin/bash

#MOLGENIS walltime=23:59:00 mem=20gb ppn=6

#string pennCNVVersion
#string intermediateDir
#string tmpTmpdir
#string tmpDir
#string workDir
#string tmpName
#string Project
#string logsDir
#string finalReport
#string samplesheet
#string pennCNVInputDir
#string pennCNVVersion
#string hmmFile
#string pfbFile
#string pennCNVDir

set -e
set -u

module load "${pennCNVVersion}"
module list

makeTmpDir "${pennCNVDir}"
tmpPennCNVDir="${MC_tmpFile}"

cd "${pennCNVInputDir}"

find "$PWD" -type f > "${tmpPennCNVDir}/${Project}.penncnv.list"

perl "${EBROOTPENNCNV}/detect_cnv.pl" \
-test -hmm "${EBROOTPENNCNV}/lib/${hmmFile}" \
-pfb "${pfbFile}" \
-list "${tmpPennCNVDir}/${Project}.penncnv.list" \
-log "${tmpPennCNVDir}/${Project}.penncnv.log" \
-out "${tmpPennCNVDir}/${Project}.rawcnv"

cd -

echo "mv ${tmpPennCNVDir}/${Project}.penncnv.list ${pennCNVDir}/${Project}.penncnv.list"
echo "mv ${tmpPennCNVDir}/${Project}.penncnv.log ${pennCNVDir}/${Project}.penncnv.log"
echo "mv ${tmpPennCNVDir}/${Project}.rawcnv ${pennCNVDir}/${Project}.rawcnv"

mv "${tmpPennCNVDir}/${Project}.penncnv.list" "${pennCNVDir}/${Project}.penncnv.list"
mv "${tmpPennCNVDir}/${Project}.penncnv.log" "${pennCNVDir}/${Project}.penncnv.log"
mv "${tmpPennCNVDir}/${Project}.rawcnv" "${pennCNVDir}/${Project}.rawcnv"

