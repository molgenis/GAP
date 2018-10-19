#MOLGENIS walltime=59:59:00 mem=20gb ppn=6

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

mkdir -p "${zCallDir}"

cd "${optiCallDir}"
bash ${EBROOTGAP}/scripts/GS_to_zCall.sh -i "${finalReport}" -o "${zCallDir}"
cd -

