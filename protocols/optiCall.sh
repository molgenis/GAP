#MOLGENIS walltime=23:59:00 mem=8gb ppn=2

#string gapVersion
#string Project
#string logsDir
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
