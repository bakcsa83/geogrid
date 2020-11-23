pipeline {
    agent any

    stages {
        stage('Build') {
            steps {
                sh 'mvn package -Dmaven.test.skip=true -U'
            }
        }
        stage('Test') {
            steps {
                sh 'mvn test -DargLine="-Xmx6G"'
            }
        }
    }
    post {
            always {
                script {
                    emailext subject: '$DEFAULT_SUBJECT',
                             body: '$DEFAULT_CONTENT',
                             recipientProviders: [
                                [$class: 'CulpritsRecipientProvider'],
                                [$class: 'DevelopersRecipientProvider'],
                                [$class: 'RequesterRecipientProvider']
                             ],
                             replyTo: '$DEFAULT_REPLYTO',
                             to: '$DEFAULT_RECIPIENTS'
                }
            }
        }
}