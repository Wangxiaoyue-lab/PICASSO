###############################################
#>File name:5_bwa_para.sh
#>Author:Cao Jun
#>Contact:caojundudu@qq.com
#>Reference:No
#>Notes:No
###############################################
   
#The script is designed to implemet the mapping in parallel
#for chipseq/atacseq/wes/wgs

#e.g.
# map_para.sh -d ~/project/chipseq -t bwa_para -r ~/reference/mouse/GRCm39.genome.fa -m bwa -W MODULE -c 2_fastp_2023-03 -n 2

###############################################

# 1 Parameters
#-d: The path of project (only have the child folder called 'rawdata')
    ##fq files should be renamed to fastq and not uncompressed
    ##fq files should be stored in pairs instead of in chaos
#-t: The name of the task file. I suggest '5_map_para'
#-r/$reference: FATSA file
#-m: The methods to map. 
    ##bwa or bowtie2
#-c: The path of clean data(trimmed by fastp or Trimmomatic). 
#-n: The number of group 
#-W: The way to load the method

while getopts ":d:t:r:m:c:n:W:" opt
do
case $opt in
d)
dir=$OPTARG
;;
t)
task=$OPTARG
;;
r)
reference=$OPTARG
;;
m)
methods=$OPTARG
;;
n)
para_core=$OPTARG
;;
c)
clean_fq=$OPTARG
;;
W)
WAY=$OPTARG
;;
esac
done

# 2 Set variables
day=$(date +'%Y-%m') 
#------------#
if [ ${WAY} == 'MODULE' ];then
module load apps/bwa/0.7.17-gnu485
else
#bwa_dir
export PATH=$PATH:${WAY}
fi
#------------#
shell_name=$0


# 3 Prepare files
if [ ! -d "${dir}/${task}" ];then
mkdir ${dir}/${task}_${day}
else
mv ${dir}/${task} ${dir}/${task}_${day}
fi

mkdir ${dir}/${task}_${day}/data
ln -s ${dir}/${clean_fq}/output ${dir}/${task}_${day}/data/cleandata
ln -s ${reference} ${dir}/${task}_${day}/data/reference
mkdir ${dir}/${task}_${day}/data/reference_i
mkdir ${dir}/${task}_${day}/output
mkdir ${dir}/${task}_${day}/script
mkdir ${dir}/${task}_${day}/script/parallel
mkdir ${dir}/${task}_${day}/script/log.out


# 3 Run the main program
number_samples=$( ls | wc -l )
number_per_group=$(( number_samples / para_core ))
mod=$(( number_samples % para_core ))
if [ mod != 0 ];then
cores=$(( para_core + 1 ))
else
cores=$para_core
fi

#Start_time=$(date +%s)
bwa index ${dir}/${task}_${day}/data/reference

for i in $(seq $cores);do
up=$((i \* number_per_group ))
j=$((i - 1))
down=$((j \* number_per_group ))
if [ $up > $number_samples ];then
up=$number_samples
fi
echo -e "#!/bin/bash\n#\n#$ -cwd\n\n#$ -S /bin/bash\n#$ -l p=15" > ${dir}/${task}_${day}/script/parallel/map_parallel_$i.sh
#echo -e "#!bin/bash\n#SBATCH -J jobname\n#SBATCH -N 1\n#SBATCH -p normal\n#SBATCH -n 15\n#SBATCH -o \n#SBATCH -e //" > ${dir}/${task}_${day}/script/parallel/map_parallel_$i.sh

echo -e 'Start_time=$(date +%s)\necho ${Start_time}' > ${dir}/${task}_${day}/script/parallel/map_parallel_$i.sh
ls ${dir}/${task}_${day}/data/cleandata | awk -v d=$down -v u=$up 'NR==d,NR==u{print $0}' | xargs -I {} echo "export PATH=$PATH:${WAY}\nmkdir ${dir}/${task}_${day}/data/reference_i/refer_$i\nln -s ${reference} ${dir}/${task}_${day}/data/reference_i/refer_$i/ref\nbwa index ${dir}/${task}_${day}/data/reference_i/refer_$i/ref\nr1=${dir}/${task}_${day}/data/cleandata/{}/*R1*\nr2=${dir}/${task}_${day}/data/cleandata/{}/*R2*\nbwa mem ${r1} ${r2} > ${dir}/${task}_${day}/output/{}.sam" >> ${dir}/${task}_${day}/script/parallel/map_parallel_$i.sh

echo -e 'End_time=$(date +%s)\nrunning_time=$End_time-$Start_time\necho ${End_time}\necho ${running_time}' >> ${dir}/${task}_${day}/script/parallel/map_parallel_$i.sh

done


#---------------------------------------#
#The samples is few
#ls ${dir}/${task}_${day}/data/cleandata | while read sample;
#do
#r1=${dir}/${task}_${day}/data/cleandata/${sample}/*R1*
#r2=${dir}/${task}_${day}/data/cleandata/${sample}/*R2*
#bwa mem ${r1} ${r2} > ${dir}/${task}_${day}/output/{}.sam
#done
#---------------------------------------#

#End_time=$(date +%s)
#running_time=$End_time-$Start_time

# 4 Record

# version
# time
# equipment
