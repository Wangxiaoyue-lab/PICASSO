#!/bin/bash
#SBATCH --job-name=combine_loom
#SBATCH -p normal
#SBATCH --mem=256G
#SBATCH -n 5
#SBATCH -N 1
#SBATCH --exclusive

# Master folder
dir=/public/home/luoliheng/yujia/counts_results

# Script path
py_path=/public/home/luoliheng/yujia/analysis/gz+jz/re_analysis_gz/velocyto/combine_loom.py

# Python path
python_path=/public/home/luoliheng/miniconda3/envs/py1/bin/python

# Combined.loom path
output_dir=/public/home/luoliheng/yujia/analysis/gz+jz/re_analysis_gz/velocyto/combined.loom

# Empty folder to paste
files=""

# Paste all the paths
for day in 15 18 45 53 54; do
    single_file=${dir}/count_gz-${day}/velocyto/count_gz-${day}.loom
    files="${files} ${single_file}"
done

### Main order
${python_path} ${py_path} ${files} -o ${output_dir}

echo "done!"
