pipeline {
	agent any

stages {
	stage('Build PR') {
    		steps {
        		checkout scm
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
}
