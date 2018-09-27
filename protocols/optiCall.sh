#MOLGENIS walltime=01:59:00 mem=2gb ppn=6

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
#string optiCallVersion
#string chr

set -e
set -u

module load "${optiCallVersion}"
module list

${EBROOTOPTICALL}/opticall \
-in ${optiCallDir}/${chr} \
-out ${optiCallDir}/${chr}
