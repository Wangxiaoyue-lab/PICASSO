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
## the suffix of fq files
fq1=$(yq e '.fq1' ${yaml})
fq2=$(yq e '.fq2' ${yaml})

## the install path of softwares
## No installation required

# 2 check the directory link data
if [ ! -d ${dir}/${input} ]; then
    mkdir ${dir}/${input}
fi
if [ ! -d ${dir}/${output} ]; then
    mkdir ${dir}/${output}
fi
mkdir ${dir}/${output}/check_md5.txt
## raw_data
ln -s ${raw_data} ${dir}/${input}/fq
ln -s ${store_md5} ${dir}/${input}/store_md5

# 3 compute md5 value
ls ${dir}/${input}/fq | while read sample;
do
md5_fq1=$(md5sum ${dir}/${input}/fq/${sample}/*${fq1} | awk '{print $1}')
md5_fq2=$(md5sum ${dir}/${input}/fq/${sample}/*${fq2} | awk '{print $1}')
if grep -q "${md5_fq1}" "${dir}/${input}/store_md5";then
    echo "${sample}/*${fq1}: md5 match" >> ${dir}/${output}/check_md5.txt
else
    echo "${sample}/*${fq1}: md5 mismatch" >> ${dir}/${output}/check_md5.txt
fi
if grep -q "${md5_fq2}" "${dir}/${input}/store_md5";then
    echo "${sample}/*${fq2}: md5 match" >> ${dir}/${output}/check_md5.txt
else
    echo "${sample}/*${fq2}: md5 mismatch" >> ${dir}/${output}/check_md5.txt
fi
done

cat ${dir}/${output}/check_md5.txt | grep "mismatch"