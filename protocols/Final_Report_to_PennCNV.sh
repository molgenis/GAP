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

set -e
set -u
set -o pipefail

module load "${pennCNVVersion}"
module list

makeTmpDir "${pennCNVInputDir}"
tmpPennCNVInputDir="${MC_tmpFile}"

mkdir -p "${pennCNVInputDir}"

perl "${EBROOTPENNCNV}/split_illumina_report.pl" \
--prefix "${tmpPennCNVInputDir}/" \
--suffix '.txt' \
"${finalReport}"

echo "mv ${tmpPennCNVInputDir}/ ${pennCNVInputDir}"
mv "${tmpPennCNVInputDir}" "${pennCNVInputDir}"
