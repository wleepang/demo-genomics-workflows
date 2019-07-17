#!/bin/bash

WORKFLOW_NAME=$1  # custom name for workflow
OVERRIDES=$2      # path to json file for job overrides, e.g. file://path/to/overrides.json

# get the name of the high-priority queue
HIGHPRIORITY_JOB_QUEUE=$(aws batch describe-job-queues | jq -r .jobQueues[].jobQueueName | grep highpriority)

# submits nextflow workflow to AWS Batch
aws batch submit-job \
    --job-definition nextflow \
    --job-name nf-workflow-${WORKFLOW_NAME} \
    --job-queue ${HIGHPRIORITY_JOB_QUEUE} \
    --container-overrides ${OVERRIDES}

