#!/bin/bash

WORKFLOW_NAME=$1  # custom name for workflow
OVERRIDES=$2      # path to json file for job overrides, e.g. file://path/to/overrides.json

AWS_REGION=$(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone| sed -r "s/(.*?)[a-z]/\1/")

# get the name of the high-priority queue
HIGHPRIORITY_JOB_QUEUE=$(aws --region $AWS_REGION batch describe-job-queues | jq -r .jobQueues[].jobQueueName | grep highpriority)

# submits nextflow workflow to AWS Batch
aws batch submit-job \
    --region $AWS_REGION \
    --job-definition nextflow \
    --job-name nf-workflow-${WORKFLOW_NAME} \
    --job-queue ${HIGHPRIORITY_JOB_QUEUE} \
    --container-overrides ${OVERRIDES}

