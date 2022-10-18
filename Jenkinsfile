node {
	checkout scm
	
	def remote = [:]
  	remote.name = 'Gearshift'
  	remote.host = 'gearshift'
  	remote.user = 'umcg-molgenis'
  	stage('Remote SSH') {
    	sshCommand remote: remote, command: "echo moi"

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
