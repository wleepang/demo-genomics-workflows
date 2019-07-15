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

# get and install nextflow
mkdir ~/bin
cd ~/bin
curl -s https://get.nextflow.io | bash

# create a nextflow config for aws-batch
CONFIG_DIR=~/environment/config
mkdir $CONFIG_DIR
DEFAULT_JOB_QUEUE=$(aws batch describe-job-queues | jq -r .jobQueues[].jobQueueName | grep default)

cat <<EOF > $CONFIG_DIR/batch.config
process.executor = "awsbatch"
process.queue = "${DEFAULT_JOB_QUEUE}"
executor.awscli = "/home/ec2-user/miniconda/bin/aws"
EOF

cd $CWD
