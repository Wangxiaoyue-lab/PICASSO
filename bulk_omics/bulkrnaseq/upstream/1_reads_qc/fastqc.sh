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

## the path of fq data
fq=$(yq e '.fq' ${yaml})
## the suffix of fq files
fq1=$(yq e '.fq1' ${yaml})
fq2=$(yq e '.fq2' ${yaml})

## the install path of softwares
fastqc=$(yq e '.fastqc' ${yaml})
multiqc=$(yq e '.multiqc' ${yaml})

## some important parameters
thread=$(yq e '.thread' ${yaml})

# 2 check the directory link data
if [ ! -d ${dir}/${input} ]; then
    mkdir ${dir}/${input}
fi
if [ ! -d ${dir}/${output} ]; then
    mkdir ${dir}/${output}
fi
## fq data
ln -s ${fq} ${dir}/${input}/fq


# 3 fastqc report
process_fastqc() {
    sample="$1"
    fq="$2"
    suffix_fq1="$3" 
    suffix_fq2="$4"
    fastqc="$5"
    output="$6"

    if [ ! -d ${output}/${sample} ]; then
        mkdir ${output}/${sample}
    fi
    ${fastqc} -o ${output}/${sample} \
        ${fq}/${sample}/*${suffix_fq1}
    ${fastqc} -o ${output}/${sample} \
        ${fq}/${sample}/*${suffix_fq2}
}

find "${dir}/${input}/fq" --mindepth 1 --maxdepth 1 -type d -print0 | \
    xargs -0 -n 1 -P ${thread} -I {} \
    bash -c 'process_fastqc "{}" "${dir}/${input}/fq" "${suffix_fq1}" "${suffix_fq2}" "${dir}/${input}/store_md5" "${dir}/${output}/md5_fq.txt" "${dir}/${output}/check_md5.txt"'


 