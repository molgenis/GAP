node {
	checkout scm
	
	stage ('Automated test') {
		sh "test/autoTestGAP.sh"
	}
	stage ('Automated test') {
		sh "test/autoTestGAP.sh"
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
