pipeline {
	agent any
	
	stages {
		stage ('Checkout SCM'){
			steps {
    				checkout([$class: 'GitSCM'],
					  [branches: [name: 'jenkins-test']]
            			])
			}
			
		}
		
	}
}
