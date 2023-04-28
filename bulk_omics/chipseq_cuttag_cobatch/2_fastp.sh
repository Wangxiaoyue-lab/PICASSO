#parameter:
#$1/$dir: the path of project (only have the child folder called 'rawdata')
##fq files should be renamed to fastq and not uncompressed
##fq files should be stored in pairs

day=$(date +'%Y-%m')
dir=$1 
fastp_dir=/public/home/caojun/software
shell_name=$0

mkdir ${dir}/2_fastp_${day}
mkdir ${dir}/2_fastp_${day}/data
ln -s ${dir}/rawdata/ ${dir}/2_fastp_${day}/data
mkdir ${dir}/2_fastp_${day}/output
mkdir ${dir}/2_fastp_${day}/script
mkdir ${dir}/2_fastp_${day}/script/log.out


ls ${dir}/2_fastp_${day}/data | while read sample;
do
  mkdir ${dir}/2_fastp_${day}/output/${sample} 
  ${fastp_dir}/fastp -i ${dir}/2_fastp_${day}/data/${sample}/*R1* -I ${dir}/2_fastp_${day}/data/${sample}/*R2* -o ${dir}/2_fastp_${day}/output/${sample}/${sample}.R1.clean.fastq.gz -O ${dir}/2_fastp_${day}/output/${sample}/${sample}.R2.clean.fastq.gz 1> ${dir}/2_fastp_${day}/script/log.out/${shell_name}_${sample}.out 2> ${dir}/2_fastp_${day}/script/log.out/${shell_name}_${sample}.err
done
