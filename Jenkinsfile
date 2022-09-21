pipeline {
	agent any

stages {
	stage('Build PR') {
	when {
		changeRequest()
	}
    	node {
        	checkout scm
    	}
	}
}
}
