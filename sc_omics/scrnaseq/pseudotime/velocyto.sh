#!/bin/bash

msk_gtf=/public/home/luoliheng/refgenome/mm10_rmsk.gtf
annotation=/public/home/luoliheng/refgenome/refdata-gex-mm10-2020-A/genes/genes.gtf

for day in 15 18 45 53 54; do
    file_path=/public/home/luoliheng/yujia/counts_results/count_gz-${day}
    job_name="velocyto_${day}"

    # Generate Slurm script for this job
    cat >${job_name}.sh <<EOL
#!/bin/bash
#SBATCH --job-name=${job_name}
#SBATCH -p normal
#SBATCH --mem=200G
#SBATCH -n 5
#SBATCH -N 1
#SBATCH --exclusive
module purge
module load apps/samtools/1.16.1-gnu485

/public/home/luoliheng/miniconda3/envs/py1/bin/velocyto run10x -m ${msk_gtf} \
    ${file_path} ${annotation} -vvv
EOL

    # Submit the job
    sbatch ${job_name}.sh
done
