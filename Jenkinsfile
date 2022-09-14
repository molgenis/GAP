pipeline {
	agent any

	stage('Build PR') {
	when {
		changeRequest()
	}
    	steps {
        	checkout scm
    		}
	}
}
