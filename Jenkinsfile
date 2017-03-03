node {
    stage 'Checkout'

    deleteDir()

    checkout scm

    stage 'Build / Install'

    sh 'cd ${WORKSPACE} && /illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh install.py ${WORKSPACE}/install --setup illumina --python system --python-interpreter /illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh --with-rtgtools --no-tests'

    stage 'Test'

    sh 'cd ${WORKSPACE}/install && echo "PYTHON=/illumina/sync/software/groups/hap.py/latest/python-ve/bin/python-wrapper.sh ${WORKSPACE}/src/sh/run_tests.sh" | qsub -l excl -cwd -sync y -N jnks-hap.py'
}
