pipeline {
	agent any
	
	stages {
		stage('Code Checkout') {
            		steps {
                		checkout([
                    		$class: 'GitSCM', 
                    		branches: [[name: '*/test-jenkins']], 
                    		userRemoteConfigs: [[url: 'https://github.com/molgenis/GAP.git']]
                		])
			}
		}
	}
}
