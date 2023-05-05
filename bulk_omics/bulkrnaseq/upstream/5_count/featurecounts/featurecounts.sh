#! /bin/bash
 
# 0 参考引用


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
featurecounts=$(yq e '.featurecounts' ${yaml})
multiqc=$(yq e '.multiqc' ${yaml})

## fq data which should be split to subfolders according to samples
fq=$(yq e '.fq' ${yaml})
## the suffix of fq files
suffix_fq1=$(yq e '.suffix_fq1' ${yaml})
suffix_fq2=$(yq e '.suffix_fq2' ${yaml})

## important parameters
index=$(yq e '.index' ${yaml})
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


# 3 construct the index of featurecounts


# 4 quantifing by featurecounts
ls ${dir}/${input}/fq | grep -v 'multiqc'| while read sample;
do
mkdir ${dir}/${output}/${sample}
${featurecounts} -T ${thread} \
	-t exon \
	-p \
	-g gene_id \
	--extraAttributes gene_name,gene_type \
	-a ${ref_gtf} \
	-o ${dir}/${output}/${sample}/${sample}.featurecounts.txt \
	${dir}/${input}/fq/${sample}/*sortedByCoord.out.bam
done


# 5 multiqc
${multiqc} ${dir}/${output}

# 6 merge result
