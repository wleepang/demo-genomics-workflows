#!/bin/bash

CWD=$(pwd)

# core requirements
sudo yum install -y jq

# nextflow requires java 8
sudo yum install -y java-1.8.0-openjdk

# The Cloud9 AMI is Amazon Linux 1, which defaults to Java 7
# change the default to Java 8 and set JAVA_HOME accordingly
sudo alternatives --set java /usr/lib/jvm/jre-1.8.0-openjdk.x86_64/bin/java
export JAVA_HOME=/usr/lib/jvm/jre-1.8.0-openjdk.x86_64

# do this so that any new shells spawned use the correct
# java version
echo "export JAVA_HOME=$JAVA_HOME" >> ~/.bash_profile

# get and install nextflow
mkdir -p ~/bin
cd ~/bin
curl -s https://get.nextflow.io | bash

# create a nextflow config for aws-batch
CONFIG_DIR=~/environment/config
mkdir -p $CONFIG_DIR

# TODO: use cloudformation outputs to get these values
DEFAULT_JOB_QUEUE=$(aws batch describe-job-queues | jq -r .jobQueues[].jobQueueName | grep default)
WORK_BUCKET=$(aws s3 ls | cut -d " " -f 3 | grep genomics-workflows)

cat <<EOF > $CONFIG_DIR/batch.config
workDir = "s3://${WORK_BUCKET}"
process.executor = "awsbatch"
process.queue = "${DEFAULT_JOB_QUEUE}"
executor.awscli = "/home/ec2-user/miniconda/bin/aws"
EOF

cd $CWD

# create a public key for ssh
ssh-keygen -f $HOME/.ssh/id_rsa -t rsa -N ''