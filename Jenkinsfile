node {
	stage ('Checkout') {
	checkout scm
	}
        sshagent(credentials : ['umcg-molgenis']) {
            sh 'ssh -A -o StrictHostKeyChecking=no umcg-molgenis@airlock.hpc.rug.nl uptime'
	    sh 'ssh umcg-molgenis@gearshift uptime'
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
