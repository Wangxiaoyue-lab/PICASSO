#!/bin/bash

# 0 参考引用
# https://combine-lab.github.io/alevin-tutorial/2019/selective-alignment/
# https://github.com/COMBINE-lab/salmon
# https://salmon.readthedocs.io/en/latest/


# 1 设置关键变量
##配置文件路径
yaml=$1
## project的task文件夹
dir=$(yq e '.dir' ${yaml})
## task里的input存储文件夹
input=$(yq e '.input' ${yaml})
## task里的output存储文件夹
output=$(yq e '.output' ${yaml})

## 参考基因组文件夹
ref=$(yq e '.ref' ${yaml})
###应该有fa、fai、gtf三个文件
ref_fa=${ref}/$(yq e '.ref_fa' ${yaml})
ref_fai=${ref}/$(yq e '.ref_fai' ${yaml})
ref_gtf=${ref}/$(yq e '.ref_gtf' ${yaml})

## 软件的安装目录
salmon=$(yq e '.salmon' ${yaml})
gffread=$(yq e '.gffread' ${yaml})
multiqc=$(yq e '.multiqc' ${yaml})


## fq数据 应该按sample分好文件夹
fq=$(yq e '.fq' ${yaml})
fq1=$(yq e '.fq1' ${yaml})
fq2=$(yq e '.fq2' ${yaml})



# 2 构建并链接数据
if [ ! -d ${dir}/${input} ]; then
    mkdir ${dir}/${input}
fi
if [ ! -d ${dir}/${output} ]; then
    mkdir ${dir}/${output}
fi

## fq数据
ln -s ${fq} ${dir}/${input}/fq

## salmon的ref
mkdir ${dir}/${input}/salmon_index
ln -s ${ref_fa} ${dir}/${input}/salmon_index/ref_fa
ln -s ${ref_fai} ${dir}/${input}/salmon_index/ref_fai
ln -s ${ref_gtf} ${dir}/${input}/salmon_index/ref_gtf


# 3 salmon的index
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
    -p 12 -k 31 --gencode



## 3 salmon的定量
ls ${dir}/${input}/fq | while read sample;
do
if [ ! -d ${dir}/${output}/${sample} ]; then
    mkdir ${dir}/${output}/${sample}
fi
${salmon} quant \
    -i ${dir}/${input}/salmon_index/ref_fa_salmon_index \
    -l A \
    -1 ${dir}/${input}/fq/${sample}/*${fq1} \
    -2 ${dir}/${input}/fq/${sample}/*${fq2} \
    -p 10 \
    --gcBias \
    --validateMappings \
    -g ${dir}/${input}/salmon_index/ref_gtf \
    -o ${dir}/${output}/${sample}
done


# 4 multiqc
${multiqc} ${dir}/${output}

# 5 merge result
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
