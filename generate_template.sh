\#!/bin/bash

module load Molgenis-Compute/v17.08.1-Java-1.8.0_74
module list

project=PROJECTNAME
runID=XX
host=$(hostname -s)
echo ${host}

projectDir=/groups/umcg-gaf/tmp03/projects/${project}/run${runID}/jobs/
genScripts=/groups/umcg-gaf/tmp03/generatedscripts/${project}
EBROOT_GAP=/home/umcg-mbenjamins/Github/GAP/


mkdir -p ${projectDir}

samplesheet="${genScripts}/${project}.csv" ; mac2unix "${samplesheet}"

perl "${EBROOT_GAP}/Scripts/convertParametersGitToMolgenis.pl" "${EBROOT_GAP}/parameters_${host}.csv" > "${genScripts}/parameters_host_converted.csv"
perl "${EBROOT_GAP}/Scripts/convertParametersGitToMolgenis.pl" "${EBROOT_GAP}/parameters.csv" > "${genScripts}/parameters_converted.csv"


sh "${EBROOTMOLGENISMINCOMPUTE}/molgenis_compute.sh" \
-p "${genScripts}/parameters_converted.csv" \
-p "${genScripts}/parameters_host_converted.csv" \
-p "${genScripts}/${project}.csv" \
-w "${EBROOT_GAP}/Prepare_GAP_workflow.csv" \
-rundir "${genScripts}/scripts" \
--runid "${runID}" \
outputdir="scripts/jobs;mainParameters=${genScripts}/parameters_converted.csv"\
-weave \
--generate
