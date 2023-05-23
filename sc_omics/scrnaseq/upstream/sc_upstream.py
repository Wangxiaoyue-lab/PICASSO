import os
import re 

# 1 fastq 命名规范
##[Sample Name]_S1_L00[Lane Number]_[Read Type]_001.fastq
#--------------------------------------------------------------
#[Sample Name]：自己对样本的命名
#[Lane Number]：同一个样本不同测序泳道的序号
#[Read Type]：四种类型，如下
#I1: Sample index read (这个文件可要可不要)
#R1: Read 1
#R2: barcode（这个文件必须要有）


def check_sc_fqs(directory):
    pattern = re.compile(r'^.+_S1_L00\d+_[^_]+_001\.fastq(\.gz)?$')
    for filename in os.listdir(directory):
        if not pattern.match(filename):
            print(f'文件 {filename} 不符合命名规则')


