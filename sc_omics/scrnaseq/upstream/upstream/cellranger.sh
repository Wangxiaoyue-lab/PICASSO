#!/bin/sh
###############################################################################
# 这个脚本用于运行 cellranger count 命令
# 它会遍历 样本名称
# 为每个样本创建一个包含 cellranger count 命令的脚本文件
# 并使用 sbatch 命令提交作业
###############################################################################

# define refrence genome
ref=/public/home/luoliheng/refgenome/refdata-gex-GRCh38-2020-A

# loop through DE HE, you can loop through any thing you like
dir_index=0
for sample_name in DE HE; do
    dir_index=$((dir_index + 1))

    dir="/public/home/luoliheng/perturb/raw_data_$sample_name""\
/CP2023030700120/H101SC23031823/RSHR01204/X101SC23031823-Z01/X101SC23031823-Z01-J00${dir_index}""\
/Rawdata/IPSC_$sample_name"

    # print the commands to scripts to run parallel
    echo "\
#!/bin/sh
#SBATCH --job-name=cellranger_$sample_name
#SBATCH -p normal
#SBATCH --mem=100G
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --exclusive
module load apps/cellranger/7.1.0
cellranger count --id=count_$sample_name \
--transcriptome=$ref \
--fastqs=$dir \
--sample=IPSC_$sample_name-1 \
--localcores=8 \
--nosecondary" >$sample_name.sh
    sbatch $sample_name.sh
done
