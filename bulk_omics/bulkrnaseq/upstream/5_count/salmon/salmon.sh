#!/bin/bash

# 0 参考引用
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
# https://github.com/COMBINE-lab/salmon
# https://salmon.readthedocs.io/en/latest/


# 1 set key variables
## path of yaml configuration
yaml=$1
## the task folder in project
dir=$(yq e '.dir' ${yaml})
## the input folder in task
input=$(yq e '.input' ${yaml})
## the output folder in task
output=$(yq e '.output' ${yaml})
## the report folder in task
report=$(yq e '.report' ${yaml})

## the path of reference
ref=$(yq e '.ref' ${yaml})
###应该有fa、fai、gtf三个文件
ref_fa=${ref}/$(yq e '.ref_fa' ${yaml})
ref_fai=${ref}/$(yq e '.ref_fai' ${yaml})
ref_gtf=${ref}/$(yq e '.ref_gtf' ${yaml})

## the install path of softwares
salmon=$(yq e '.salmon' ${yaml})
gffread=$(yq e '.gffread' ${yaml})
multiqc=$(yq e '.multiqc' ${yaml})

## fq data which should be split to subfolders according to samples
fq=$(yq e '.fq' ${yaml})
## the suffix of fq files
suffix_fq1=$(yq e '.suffix_fq1' ${yaml})
suffix_fq2=$(yq e '.suffix_fq2' ${yaml})

## important parameters
thread=$(yq e '.thread' ${yaml})


# 2 check the directory and link data
if [ ! -d ${dir}/${input} ]; then
    mkdir ${dir}/${input}
fi
if [ ! -d ${dir}/${output} ]; then
    mkdir ${dir}/${output}
fi
if [ ! -d ${dir}/${report} ]; then
    mkdir ${dir}/${report}
fi
## fq data
ln -s ${fq} ${dir}/${input}/fq
## the ref of salmon 
mkdir ${dir}/${input}/salmon_index
ln -s ${ref_fa} ${dir}/${input}/salmon_index/ref_fa
ln -s ${ref_fai} ${dir}/${input}/salmon_index/ref_fai
ln -s ${ref_gtf} ${dir}/${input}/salmon_index/ref_gtf
## the report file


# 3 construct the index of salmon 
## extract the cRNA sequence
${gffread} \
    ${dir}/${input}/salmon_index/ref_gtf \
    -g ${dir}/${input}/salmon_index/ref_fa \
    -w ${dir}/${input}/salmon_index/ref_fa.tmp
## clean the cDNA fa
cut -f 1 -d ' ' ${dir}/${input}/salmon_index/ref_fa.tmp \
    > ${dir}/${input}/salmon_index/ref_fa_cDNA.fa
## extract the name of sequence
grep '^>' ${dir}/${input}/salmon_index/ref_fa | cut -d ' ' -f 1 | \
    sed 's/^>//g' > ${dir}/${input}/salmon_index/ref_fa.decoys.txt
## combine the cDNA sequence and genome sequence
cat ${dir}/${input}/salmon_index/ref_fa_cDNA.fa \
    ${dir}/${input}/salmon_index/ref_fa \
    > ${dir}/${input}/salmon_index/ref_fa_genome.fa
## construct the index
${salmon} index \
    -t ${dir}/${input}/salmon_index/ref_fa_genome.fa \
    -d ${dir}/${input}/salmon_index/ref_fa.decoys.txt \
    -i ${dir}/${input}/salmon_index/ref_fa_salmon_index \
    -p ${thread} -k 31 --gencode


# 4 quantifing by salmon
ls ${dir}/${input}/fq | while read sample;
do
if [ ! -d ${dir}/${output}/${sample} ]; then
    mkdir ${dir}/${output}/${sample}
fi
${salmon} quant \
    -i ${dir}/${input}/salmon_index/ref_fa_salmon_index \
    -l A \
    -1 ${dir}/${input}/fq/${sample}/*${suffix_fq1} \
    -2 ${dir}/${input}/fq/${sample}/*${suffix_fq2} \
    -p ${thread} \
    --gcBias \
    --validateMappings \
    -g ${dir}/${input}/salmon_index/ref_gtf \
    -o ${dir}/${output}/${sample}
done


# 5 multiqc
${multiqc} ${dir}/${output}

# 6 merge result
files=$( ls ${dir}/${input}/fq | paste -sd ',' )
cd ${dir}/${output}/
${salmon} quantmerge \
    --quants {${files}} \
    --names {${files}} \
    --column=NumReads \
    -o ${dir}/${output}/merge_salmon_count.txt 
${salmon} quantmerge \
    --quants {${files}} \
    --names {${files}} \
    --column=tpm \
    -o ${dir}/${output}/merge_salmon_tpm.txt 

# 7 report the task
Start=$(data +%s)
End=$(data +%s)
runtime$((End-Start))

${salmon} -v
${gffread} -v
${multiqc} -v