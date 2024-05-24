node {
	stage ('Checkout') {
		checkout scm
	}
    	stage ('Automated test') {
        
        echo "Copy test from repo to molgenis home on Talos"
        sh "sudo scp test/autoTestGAP.sh reception+talos:/home/umcg-molgenis/"
        
        echo "Login to Talos"
	    
	sh '''
            sudo ssh -tt reception+talos 'exec bash -l << 'ENDSSH'
	    	echo "Starting automated test"
		bash /home/umcg-molgenis/autoTestGAP.sh '''+env.CHANGE_ID+'''
ENDSSH'
        '''	
	}
}
