#! /bin/bash
# this script is ran in SGE system, I haven't tried in slurm system
b=0
i="compute-0-0 compute-0-1 compute-0-16 compute-0-17 compute-0-18 compute-0-19 compute-0-2 compute-0-20 compute-0-21 compute-0-22 compute-0-26 compute-0-27 compute-0-28 compute-0-29  compute-0-30 compute-0-31 compute-0-32 compute-0-4 compute-0-5 compute-0-6 compute-0-7"
# compute-0-3 不足10核
cat /data/luoliheng/yujia/NAFLD-singleCell/cell_count/namelist | while read ID; do

    a=$(echo $i | awk -v b=$(($b + 1)) '{print $b}')
    b=$(($b + 1))
    h5files=/data/luoliheng/yujia/NAFLD-singleCell/cell_count/counts_results/count_$ID/outs/raw_feature_bc_matrix.h5

    echo "#! /bin/bash
#$ -cwd
#$ -S /bin/bash
#$ -V
#$ -l p=10
#$ -l h=$a

cellbender remove-background \
                 --input $h5files \
                 --output .\output$ID.h5 \
                 --expected-cells 5000 \
                 --total-droplets-included 15000 \
                 --fpr 0.01 \
                 --epochs 150
" >$ID-$a.sh
    qsub $ID-$a.sh
done
