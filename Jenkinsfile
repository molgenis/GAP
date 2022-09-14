pipeline {
	agent any

stages {
	stage('Build PR') {
	when {
		changeRequest()
	}
    	steps {
        	checkout scm
    		}
	}
}
}
