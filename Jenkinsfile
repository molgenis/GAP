node {
	stage ('Checkout') {
		checkout scm
	}
    stage ('Automated test') {
        
        echo "Copy test from repo to molgenis home on Gearshift"
        sh "sudo scp test/autoTestGAP.sh airlock+gearshift:/home/umcg-molgenis/"
        
        echo "Login to Gearshift"
        sh '''sudo ssh -T airlock+gearshift <<ENDSSH
		echo "Starting automated test"
		sh autoTestGAP.sh '''+env.CHANGE_ID+'''
ENDSSH
'''
	}
	stage('ShellCheck') {
		sh "check/shellcheck.sh"			
	}
	stage('IndentationCheck') {
		sh "check/indentationcheck.sh"
	}	
	post {
		always {
		recordIssues (enabledForFailure: true, failOnError: true, qualityGates: [[threshold: 1, type: 'TOTAL', unstable: false]], tools: [checkStyle(name: 'ShellCheck')], trendChartType: 'NONE')
		}
	}
}
