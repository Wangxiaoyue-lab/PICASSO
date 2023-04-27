#!/bin/sh

grna=/public/home/luoliheng/perturb/crop_count/mini_library.txt
for sample_name in DE HE; do

    bamfiles=/public/home/luoliheng/perturb/cellranger/count_$sample_name/outs/possorted_genome_bam.bam

    echo "\
#!/bin/sh
#SBATCH --job-name=cellranger_$sample_name
#SBATCH -p normal
#SBATCH --mem=100G
#SBATCH -n 8
#SBATCH -N 1
#SBATCH --exclusive
/public/home/luoliheng/perturb/crop_count/cropseq_count.py --lib-grna $grna --file-type bam \
--files $bamfiles -m 0 \
-n $sample_name-sgRNA_m0" >$sample_name.sh
    sbatch $sample_name.sh

done
