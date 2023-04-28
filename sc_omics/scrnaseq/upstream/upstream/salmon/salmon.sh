#用salmon进行计量

module purge
module load apps/salmon/1.10.0-x86_64

r1=
r2=
index=#index文件夹路径

output=#输出结果会在一个指定文件夹中

salmon quant -i ${index} -l A -p 6 --validateMappings \
         --gcBias --numGibbsSamples 20 -o ${output} \
         -1 ${r1} -2 ${r2}