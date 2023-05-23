#!/bin/sh
###############################################################################
# 这个脚本用于运行 cellranger count 命令
# 它会遍历 样本名称
# 为每个样本创建一个包含 cellranger count 命令的脚本文件
# 并使用 sbatch 命令提交作业
###############################################################################

# software
installed_cellranger=//




# define refrence genome
ref=///

# the stucture of folder
project=///
task=//
rawdata=${project}/${task}/input/rawdata
#--rawdata
#----sample1
#------test_sample_S1_L001_R1_001.fastq.gz
#------test_sample_S1_L001_R2_001.fastq.gz
#------...
#----sample2
#------...
#------...
#----sample3
#------...
#------...

mkdir ${project}/${task}/output/cellranger_result\n

# print the commands to scripts to run parallel
ls ${rawdata} | while read sample;
do
    sample_name=$(ls ${rawdata}/${sample} | head -n 1 | sed 's/_S1.*//')
    echo -e "\
    #!/bin/sh\n
    #SBATCH --job-name=cr_${sample}\n
    #SBATCH -p normal\n
    #SBATCH --mem=100G\n
    #SBATCH -n 8\n
    #SBATCH -N 1\n
    module load apps/cellranger/7.1.0\n
    mkdir ${project}/${task}/output/cellranger_result/${sample}\n
    cd ${project}/${task}/output/cellranger_result/${sample}
    cellranger count --id=cellranger_result_${sample} \
    --transcriptome=${ref} \
    --fastqs=${rawdata}/${sample} \
    --sample=${sample_name} \
    --localcores=8 \
    --nosecondary" > ${project}/${task}/script/cellranger_run_${sample}.sh\n
    sbatch  ${project}/${task}/script/cellranger_run_${sample}.sh
done




