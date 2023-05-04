#! /bin/bash

# 0 reference


# 1 set key variables
## path of yaml configuration
yaml=$1
## the task folder in project
dir=$(yq e '.dir' ${yaml})
## the input folder in task
input=$(yq e '.input' ${yaml})
## the output folder in task
output=$(yq e '.output' ${yaml})

## the path of md5 file
store_md5=$(yq e '.store_md5' ${yaml})
## the path of raw data
raw_data=$(yq e '.raw_data' ${yaml})

## the install path of softwares
fastqc=$(yq e '.fastqc' ${yaml})
multiqc=$(yq e '.multiqc' ${yaml})