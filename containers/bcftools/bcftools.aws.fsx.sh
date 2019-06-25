#!/bin/bash
set -e

## WORKFLOW_NAME and EXECUTION_ID are environment variables that need to be set
## in the AWS Batch Job Submission.  This is used as the unique identifier for
## the workflow so that any outputs generated are isolated from other
## concurrently running workflows.
if [ -z "$WORKFLOW_NAME" ]; then
    WORKFLOW_NAME="workflow"
fi

if [ -z "$EXECUTION_ID" ]; then
    EXECUTION_ID="0"
fi

WORKFLOW_ID="$WORKFLOW_NAME/$EXECUTION_ID"

COMMAND=$1
REFERENCE_NAME=$2
SAMPLE_ID=$3
INPUT_PATH=${4:-"./working/$WORKFLOW_ID"}
OUTPUT_PATH=${5:-"./working/$WORKFLOW_ID"}
REFERENCE_PATH=./reference

if [ ! -d "$OUTPUT_PATH" ]; then
    mkdir -p $OUTPUT_PATH
fi

function call() {
    bcftools call \
        -m \
        --threads 16 \
        -t chr21 \
        -o $OUTPUT_PATH/${SAMPLE_ID}.vcf \
        $INPUT_PATH/${SAMPLE_ID}.mpileup.vcf.gz
}

function mpileup() {
    bcftools mpileup \
        --threads 16 \
        -Oz \
        -r chr21 \
        -f $REFERENCE_PATH/${REFERENCE_NAME}.fasta \
        $INPUT_PATH/${SAMPLE_ID}.bam \
        > $OUTPUT_PATH/${SAMPLE_ID}.mpileup.vcf.gz
}


$COMMAND
