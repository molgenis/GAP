#MOLGENIS walltime=01:59:00 mem=1gb ppn=1

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

set -e
set -u

module load "${pythonVersion}"
module load "${beadArrayVersion}"
module load "${gapVersion}"
module list

python ${EBROOTGAP}/scripts/gtc_final_report.py \
--manifest "${bpmFile}" \
--samplesheet "${samplesheet}" \
--gtc_directory "${projectRawTmpDataDir}" \
--output_file "${finalReport}"
