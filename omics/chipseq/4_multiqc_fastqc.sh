#parameter:
#$1/$dir: the path of project (only have the child folder called 'rawdata')
##fq files should be renamed to fastq and not uncompressed
##fq files should be stored in pairs

day=$(date +'%Y-%m')
dir=$1 
shell_name=$0
multiqc_path=/public/home/caojun/anaconda3/envs/py3.8/bin

mkdir ${dir}/4_multiqc_fastqc_${day}
mkdir ${dir}/4_multiqc_fastqc_${day}/data
ln -s ${dir}/1_fastqc_${day}/output/ ${dir}/4_multiqc_fastqc_${day}/data/before
ln -s ${dir}/3_fastqc_${day}/output/ ${dir}/4_multiqc_fastqc_${day}/data/after
mkdir ${dir}/4_multiqc_fastqc_${day}/output
mkdir ${dir}/4_multiqc_fastqc_${day}/script
mkdir ${dir}/4_multiqc_fastqc_${day}/script/log.out

${multiqc_path}/multiqc ${dir}/4_multiqc_fastqc_${day}/data/before

${multiqc_path}/multiqc ${dir}/4_multiqc_fastqc_${day}/data/after

