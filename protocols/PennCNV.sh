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

cd "${pennCNVInputDir}"

find $PWD -type f > ${pennCNVDir}/${Project}.penncnv.list

perl ${EBROOTPENNCNV}/detect_cnv.pl \
-test -hmm ${EBROOTPENNCNV}/lib/${hmmFile} \
-pfb ${pfbFile} \
-list ${pennCNVDir}/${Project}.penncnv.list \
-log ${pennCNVDir}/${Project}.penncnv.log \
-out ${pennCNVDir}/${Project}.rawcnv

cd -
