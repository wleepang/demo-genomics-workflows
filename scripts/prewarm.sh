#!/bin/bash
AWS_REGION=$(curl -s http://169.254.169.254/latest/meta-data/placement/availability-zone| sed -r "s/(.*?)[a-z]/\1/")

environments=$(aws --region $AWS_REGION batch describe-compute-environments)
ce_ondemand=$(echo $environments | jq -r .computeEnvironments[].computeEnvironmentName | grep ondemand)
ce_spot=$(echo $environments | jq -r .computeEnvironments[].computeEnvironmentName | grep spot)

aws batch \
    update-compute-environment \
    --region $AWS_REGION \
    --compute-environment $ce_ondemand \
    --compute-resources minvCpus=2,maxvCpus=100,desiredvCpus=2

aws batch \
    update-compute-environment \
    --region $AWS_REGION \
    --compute-environment $ce_spot \
    --compute-resources minvCpus=16,maxvCpus=100,desiredvCpus=16
