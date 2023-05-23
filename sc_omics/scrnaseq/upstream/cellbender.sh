#! /bin/bash

# key variable
project=///
task=//
result_cellranger=///
installed_cellbender=///

mkdir ${project}/${task}/output/result_cellbender


# run cellbender
ls ${result_cellranger} | while read sample; 
do
    h5file=${result_cellranger}/${sample}/outs/raw_feature_bc_matrix.h5
    echo -e "#! /bin/bash\n
        #SBATCH --job-name=cb_${sample}\n
        #SBATCH -p normal\n
        #SBATCH --mem=100G\n
        #SBATCH -n 8\n
        #SBATCH -N 1\n
    ${installed_cellbender} remove-background \
                    --input ${h5file} \
                    --output ${project}/${task}/output/result_cellbender/cellbender_${sample}.h5 \
                    --expected-cells 5000 \
                    --total-droplets-included 15000 \
                    --fpr 0.01 \
                    --epochs 150
    " > ${project}/${task}/script/cellbender_run_${sample}.sh
   sbatch ${project}/${task}/script/cellbender_run_${sample}.sh
done
