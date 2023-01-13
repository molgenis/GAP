node {
	stage ('Checkout') {
		checkout scm
	}
    stage ('Automated test') {
        
        echo "Copy test from repo to molgenis home on Gearshift"
        sh "sudo scp test/autoTestGAP.sh airlock+gearshift:/home/umcg-molgenis/"
        
        echo "Login to Gearshift"
	    
	sh '''
            sudo ssh -tt airlock+gearshift 'exec bash -l << 'ENDSSH'
	    	echo "Starting automated test"
		bash /home/umcg-molgenis/autoTestGAP.sh '''+env.CHANGE_ID+'''
ENDSSH'
        '''	
	}
}
