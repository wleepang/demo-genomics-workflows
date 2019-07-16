#!/bin/bash

BUCKET=$1
PREFIX="meats.txt"

mkdir -p .tmp && cd $_
curl -s "https://baconipsum.com/api/?type=meat&paras=1&format=text" > meats.txt
aws --region us-east-1 s3 cp ./meats.txt s3://${BUCKET}/${PREFIX}
cd .. && rm -rf .tmp
