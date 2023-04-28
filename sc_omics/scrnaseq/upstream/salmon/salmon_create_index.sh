#create salmon index
module purge
module load apps/salmon/1.10.0-x86_64

genome="GRCm38.primary_assembly.genome.fa.gz"
transcript="gencode.vM23.transcripts.fa.gz"
decoy="decoys.txt"
gentrome="gentrome.fa.gz"
salmon_index=#最后构建的index是一个文件夹

#Preparing
#create decoys file
grep "^>" <(gunzip -c ${genome}) | cut -d " " -f 1 > ${decoy}
sed -i.bak -e 's/>//g' ${decoy}
#create concatenated transcriptome and genome reference file
cat ${transcript} ${genome} > ${gentrome}

#Salmon Indexing
salmon index -t ${gentrome} -d ${decoy} -p 12 -i ${salmon_index} --gencode

