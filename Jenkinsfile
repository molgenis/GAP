pipeline {
	agent any

stages {
	stage('Build PR') {
    		node {
        		checkout scm
    			}
		steps {
			sh "test/autoTestGAP.sh"
		}
		}
	stage('ShellCheck') {
		steps {
			sh "check/shellcheck.sh"
			}
		}
	stage('IndentationCheck') {
		steps {
			sh "check/indentationcheck.sh"
			}
		}
	}
post {
	always {
		script {
			recordIssues (enabledForFailure: true, failOnError: true, qualityGates: [[threshold: 1, type: 'TOTAL', unstable: false]], tools: [checkStyle(name: 'ShellCheck')], trendChartType: 'NONE')
		}
	}
}
}
