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
suffix_fq1=$(yq e '.suffix_fq1' ${yaml})
suffix_fq2=$(yq e '.suffix_fq2' ${yaml})

## the install path of softwares
## No installation required

## some important parameters
thread=$(yq e '.thread' ${yaml})

# 2 check the directory and link data
if [ ! -d ${dir}/${input} ]; then
    mkdir ${dir}/${input}
fi
if [ ! -d ${dir}/${output} ]; then
    mkdir ${dir}/${output}
fi
touch ${dir}/${output}/check_md5.txt
touch ${dir}/${output}/md5_fq.txt

## raw_data
ln -s ${raw_data} ${dir}/${input}/fq
ln -s ${store_md5} ${dir}/${input}/store_md5

# 3 compute md5 value
check_md5() {
    sample="$1" # the name of sample folder
    fq="$2" # the path of all fqs
    suffix_fq1="$3" 
    suffix_fq2="$4"
    md5_input="$5" # the ref txt of md5 provided
    md5_output="$6" # the result of currect files
    check_result="$7" # the result of compare

    md5_fq1=$(md5sum ${fq}/${sample}/*${suffix_fq1} | awk '{print $1}')
    md5_fq2=$(md5sum ${fq}/${sample}/*${suffix_fq2} | awk '{print $1}')
    cat ${md5_fq1} >> ${md5_output}
    cat ${md5_fq2} >> ${md5_output}
    if grep -q "${md5_fq1}" "${md5_input}";then
        echo "${sample}/*${suffix_fq1}: md5 match" >> ${check_result}
    else
        echo "${sample}/*${suffix_fq1}: md5 mismatch" >> ${check_result}
    fi
    if grep -q "${md5_fq2}" "${md5_input}";then
        echo "${sample}/*${suffix_fq2}: md5 match" >> ${check_result}
    else
        echo "${sample}/*${suffix_fq2}: md5 mismatch" >> ${check_result}
    fi
}
export -f check_md5

find "${dir}/${input}/fq" --mindepth 1 --maxdepth 1 -type d -print0 | \
    xargs -0 -n 1 -P ${thread} -I {} \
    bash -c 'check_md5 "{}" "${dir}/${input}/fq" "${suffix_fq1}" "${suffix_fq2}" "${dir}/${input}/store_md5" "${dir}/${output}/md5_fq.txt" "${dir}/${output}/check_md5.txt"'


#ls ${dir}/${input}/fq | while read sample;
#do
#md5_fq1=$(md5sum ${dir}/${input}/fq/${sample}/*${fq1} | awk '{print $1}')
#md5_fq2=$(md5sum ${dir}/${input}/fq/${sample}/*${fq2} | awk '{print $1}')
#cat ${md5_fq1} >> ${dir}/${output}/md5_fq.txt
#cat ${md5_fq2} >> ${dir}/${output}/md5_fq.txt
#if grep -q "${md5_fq1}" "${dir}/${input}/store_md5";then
#    echo "${sample}/*${fq1}: md5 match" >> ${dir}/${output}/check_md5.txt
#else
#    echo "${sample}/*${fq1}: md5 mismatch" >> ${dir}/${output}/check_md5.txt
#fi
#if grep -q "${md5_fq2}" "${dir}/${input}/store_md5";then
#    echo "${sample}/*${fq2}: md5 match" >> ${dir}/${output}/check_md5.txt
#else
#    echo "${sample}/*${fq2}: md5 mismatch" >> ${dir}/${output}/check_md5.txt
#fi
#done

cat ${dir}/${output}/check_md5.txt | grep "mismatch"