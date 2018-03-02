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

umask 0007

module load ${computeVersion}
module list


#Create ProjectDirs
mkdir -p "${intermediateDir}"
mkdir -p "${resultDir}"
mkdir -p "${projectJobsDir}"
mkdir -p "${projectRawTmpDataDir}"


#Create Symlinks

rocketPoint=$(pwd)
host=$(hostname -s)
EBROOT_GAP=/home/umcg-mbenjamins/Github/GAP/

cd "${projectRawTmpDataDir}"

max_index=${#SentrixPosition_A[@]}-1

for i in ${SentrixBarcode_A[@]}
do
	for ((samplenumber = 0; samplenumber <= max_index; samplenumber++))
	do
		ln -sf "../../../../../rawdata/array/${i}/${i}_${SentrixPosition_A[samplenumber]}.gtc" \
		"${projectRawTmpDataDir}/${i}_${SentrixPosition_A[samplenumber]}.gtc"


		ln -sf "../../../../../rawdata/array/${i}/${i}_${SentrixPosition_A[samplenumber]}.gtc.md5" \
		"${projectRawTmpDataDir}/${i}_${SentrixPosition_A[samplenumber]}.gtc.md5"
	done
done

#
# Execute MOLGENIS/compute to create job scripts to analyse this project.
#

cd "${rocketPoint}"


perl "${EBROOT_GAP}/Scripts/convertParametersGitToMolgenis.pl" "${EBROOT_GAP}/parameters_${host}.csv" > "${rocketPoint}/parameters_host_converted.csv"
perl "${EBROOT_GAP}/Scripts/convertParametersGitToMolgenis.pl" "${EBROOT_GAP}/parameters.csv" > "${rocketPoint}/parameters_converted.csv"


sh "${EBROOTMOLGENISMINCOMPUTE}/molgenis_compute.sh" \
-p "${genScripts}/parameters_converted.csv" \
-p "${genScripts}/parameters_host_converted.csv" \
-p "${genScripts}/${Project}.csv" \
-rundir "${projectJobsDir}" \
-w "${EBROOT_GAP}/GAP_workflow.csv" \
-b slurm \
-g \
-weave \
-runID "run${runid}"
