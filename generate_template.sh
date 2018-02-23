#!/bin/bash

module load Molgenis-Compute/v17.08.1-Java-1.8.0_74
module list

project=${PROJECT}
runID=${RUN}
host=$(hostname -s)
echo ${host}

projectDir=${PROJECTDIR}
genScripts=${GENERATED_SCRIPTSDir}

mkdir -p ${projectDir}

samplesheet="${genScripts}/${project}.csv" ; mac2unix "${samplesheet}"

perl "/groups/umcg-gaf/tmp03/umcg-mbenjamins/GSA/Scripts/convertParametersGitToMolgenis.pl" "/groups/umcg-gaf/tmp03/umcg-mbenjamins/GSA/parameters_${host}.csv" > "${genScripts}/parameters_host_converted.csv"
perl "/groups/umcg-gaf/tmp03/umcg-mbenjamins/GSA/Scripts/convertParametersGitToMolgenis.pl" "/groups/umcg-gaf/tmp03/umcg-mbenjamins/GSA/parameters.csv" > "${genScripts}/parameters_converted.csv"


sh "${EBROOTMOLGENISMINCOMPUTE}/molgenis_compute.sh" \
-p "${genScripts}/parameters_converted.csv" \
-p "${genScripts}/parameters_host_converted.csv" \
-p "${genScripts}/${project}.csv" \
-w "/groups/umcg-gaf/tmp03/umcg-mbenjamins/GSA/workflow.csv" \
-rundir "${projectDir}" \
-b slurm \
-weave \
--generate
