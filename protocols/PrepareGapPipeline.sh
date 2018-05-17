#list SentrixBarcode_A
#list SentrixPosition_A
#string projectRawTmpDataDir
#string intermediateDir
#string resultDir
#string computeVersion
#string Project
#string projectJobsDir
#string projectRawTmpDataDir
#string genScripts
#string gapVersion
#string pipeline

umask 0007

module load ${computeVersion}
module load ${gapVersion}
module list


#Create ProjectDirs
mkdir -p "${intermediateDir}"
mkdir -p "${resultDir}"
mkdir -p "${projectJobsDir}"
mkdir -p "${projectRawTmpDataDir}"


#Create Symlinks

rocketPoint=$(pwd)
host=$(hostname -s)

cd "${projectRawTmpDataDir}"

max_index=${#SentrixPosition_A[@]}-1

for i in ${SentrixBarcode_A[@]}
do
	for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
	do
		ln -sf "../../../../../rawdata/array/GTC/${i}/${i}_${SentrixPosition_A[samplenumber]}.gtc" \
		"${projectRawTmpDataDir}/${i}_${SentrixPosition_A[samplenumber]}.gtc"

		ln -sf "../../../../../rawdata/array/GTC/${i}/${i}_${SentrixPosition_A[samplenumber]}.gtc.md5" \
		"${projectRawTmpDataDir}/${i}_${SentrixPosition_A[samplenumber]}.gtc.md5"
	done
done

#
# Execute MOLGENIS/compute to create job scripts to analyse this project.
#

cd "${rocketPoint}"


perl "${EBROOTGAP}/Scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/parameters_${host}.csv" > "${rocketPoint}/parameters_host_converted.csv"
perl "${EBROOTGAP}/Scripts/convertParametersGitToMolgenis.pl" "${EBROOTGAP}/${pipeline}_parameters.csv" > "${rocketPoint}/parameters_converted.csv"


sh "${EBROOTMOLGENISMINCOMPUTE}/molgenis_compute.sh" \
-p "${genScripts}/parameters_converted.csv" \
-p "${genScripts}/parameters_host_converted.csv" \
-p "${genScripts}/${Project}.csv" \
-rundir "${projectJobsDir}" \
-w "${EBROOTGAP}/diagnostics_workflow.csv" \
-b slurm \
-g \
-weave \
-runid "${runID}"
