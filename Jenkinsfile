node {
    stage 'Build and Test'
        echo "building..."
        sh 'echo ${PATH}'
        sh 'echo ${WORKSPACE}'
        sh 'LD_LIBRARY_PATH="" /illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh install.py ${WORKSPACE}/install --setup illumina --python system --python-interpreter /illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh'
    
    stage 'Notify'
        emailext body: 'The build was run.', recipientProviders: [[$class: 'CulpritsRecipientProvider'], [$class: 'DevelopersRecipientProvider']], subject: 'Build notification'
}