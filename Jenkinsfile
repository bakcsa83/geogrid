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
}