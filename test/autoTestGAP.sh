set -e
set -u

function preparePipeline(){

	local _projectName="NIST_TRIO"
	rm -f "${tmpfolder}/logs/${_projectName}/run01.pipeline.finished"
	echo "TMPFOLDER: ${tmpfolder}"
	pwd
	rsync -r --verbose --recursive --links --no-perms --times --group --no-owner --devices --specials ${pipelinefolder}/test/rawdata/203693990030 ${tmpfolder}/rawdata/array/GTC/

	rm -rf "${tmpfolder}/generatedscripts/${_projectName}/"

	rm -rf "${tmpfolder}/projects/GAP/${_projectName}/"

	rm -rf "${tmpfolder}/tmp/GAP/${_projectName}/"

	rm -rf "${tmpfolder}/logs/${_projectName}/"

	genFolder="${tmpfolder}/generatedscripts/${_projectName}/"
	mkdir "${genFolder}"
	mkdir "${tmpfolder}/logs/${_projectName}/"

	cp "${pipelinefolder}/templates/generate_template.sh" "${genFolder}/generate_template.sh"

	EBROOTGAP="${pipelinefolder}"
	echo "${EBROOTGAP}" 
	## Grep used version of molgenis compute out of the parameters file

	cp "${pipelinefolder}/test/${_projectName}.csv" "${genFolder}"
	perl -pi -e "s|/groups/umcg-gsad/tmp01/|${tmpfolder}/|g" "${genFolder}/${_projectName}.csv"
	
	cd "${genFolder}"
	perl -pi -e 's|workflow=\${EBROOTGAP}/workflow_diagnostics.csv|workflow=\${EBROOTGAP}/test_workflow.csv|' "${genFolder}/generate_template.sh"
	perl -pi -e 's|\${EBROOTGAP}|\'${pipelinefolder}'|g'  "${genFolder}/generate_template.sh"
	perl -pi -e 's|\${runID}|run01|' "${genFolder}/generate_template.sh"
	
	sed -i 's|gapVersion.*|gapVersion=GAP/betaAutotest;\\|' "${genFolder}/generate_template.sh"		

	# changing values in the diagnostics parameters file
	echo "bedToolsVersion,BEDTools/2.30.0-GCCcore-11.3.0" >> "${pipelinefolder}/parameters_diagnostics.csv"
	echo "ngsUtilsVersion,ngs-utils/22.10.1" >> "${pipelinefolder}/parameters_diagnostics.csv"

	module load  Molgenis-Compute/v19.01.1-Java-11-LTS
	pwd
	bash generate_template.sh -p "${_projectName}"
	cd scripts

	perl -pi -e 's|\${EBROOTGAP}|\'${pipelinefolder}'|g'  "PrepareGapPipeline_0.sh"

	bash submit.sh

	cd "${tmpfolder}/projects/GAP/${_projectName}/run01/jobs/"
	bash submit.sh --qos=regular
}

function checkIfFinished(){
	local _projectName="NIST_TRIO"
	count=0
	minutes=0
	while [ ! -f "${tmpfolder}/projects/GAP/${_projectName}/run01/jobs/s05_autoTestGAPResults_0.sh.finished" ]
	do

		echo "${_projectName} is not finished in ${minutes} minutes, sleeping for 2 minutes"
		sleep 120
		minutes=$((minutes+2))

		count=$((count+2))
		if [ "${count}" -eq 30 ]
		then
			echo "the test was not finished within 30 minutes, let's kill it"
			echo -e "\n"
			for i in $(ls "${tmpfolder}/projects/GAP/${_projectName}/run01/jobs/"*.sh)
			do
				if [ ! -f "${i}.finished" ]
				then
					echo "$(basename ${i}) is not finished"
				fi
			done
			exit 1
		fi
	done
	echo ""
	echo "${_projectName} test succeeded!"
	echo ""
}

tmpdirectory="tmp08"
groupName="umcg-gsad"

pipelinefolder="/groups/${groupName}/${tmpdirectory}/tmp/GAP/betaAutotest/"
tmpfolder="/groups/${groupName}/${tmpdirectory}"

if [ -d "${pipelinefolder}" ]
then
	rm -rf "${pipelinefolder}"
	echo "removed ${pipelinefolder}"
	mkdir "${pipelinefolder}"
fi
cd "${pipelinefolder}"

echo "pr number: $1"

PULLREQUEST="${1}"

git clone https://github.com/molgenis/GAP.git
cd GAP
git fetch --tags --progress https://github.com/molgenis/GAP/ +refs/pull/*:refs/remotes/origin/pr/*
COMMIT=$(git rev-parse refs/remotes/origin/pr/$PULLREQUEST/merge^{commit})
echo "checkout commit: COMMIT"
pwd
git checkout -f "${COMMIT}"

mv * ../
cd ..
rm -rf GAP/

### create testworkflow
cp "${pipelinefolder}/workflow_diagnostics.csv" test_workflow.csv
tail -1 workflow_diagnostics.csv | perl -p -e 's|,|\t|g' | awk '{print "s05_autoTestGAPResults,test/protocols/autoTestGAPResults.sh,"$1}' >> test_workflow.csv


cd "${pipelinefolder}"
pwd
mkdir -p /home/umcg-molgenis/GAP/vcf/
mkdir -p /home/umcg-molgenis/GAP/PennCNV_reports/

cp "${pipelinefolder}/test/results/vcf/"*".vcf"* "/home/umcg-molgenis/GAP/vcf/"
cp "${pipelinefolder}/test/results/PennCNV_reports/"*"_TRUE.txt" "/home/umcg-molgenis/GAP/PennCNV_reports/"
cp "${pipelinefolder}/test/autoTestArray.bed" "/home/umcg-molgenis/GAP/"

cp "${pipelinefolder}/test/results/Callrates_NIST_TRIO_TRUE.txt" "/home/umcg-molgenis/GAP/"
cp "${pipelinefolder}/test/results/NIST_TRIO_PennCNV_TRUE.txt" "/home/umcg-molgenis/GAP/"

preparePipeline
checkIfFinished
