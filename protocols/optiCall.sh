#MOLGENIS walltime=01:59:00 mem=2gb ppn=6

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
