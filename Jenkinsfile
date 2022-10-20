node {
	stage ('Checkout') {
	checkout scm
	}
        stage ('Automated test') {
		sh '''
		echo "ghprbPullId=$ghprbPullId"
        echo "sha1=$sha1"

		echo "Login to Gearshift"
         	sudo ssh -tt airlock+gearshift 'bash -s << 'ENDSSH'
		echo "Starting automated test"
		sh /home/umcg-molgenis/test_GAP.sh ${env.CHANGE_ID}
		ENDSSH'
		'''
	}
	stage('ShellCheck') {
		sh "check/shellcheck.sh"			
	}
	stage('IndentationCheck') {
		sh "check/indentationcheck.sh"
	}
	
	try {
		echo "TEST SUCCESSFULL"
	} finally {
		recordIssues (enabledForFailure: true, failOnError: true, qualityGates: [[threshold: 1, type: 'TOTAL', unstable: false]], tools: [checkStyle(name: 'ShellCheck')], trendChartType: 'NONE')
	}
}
