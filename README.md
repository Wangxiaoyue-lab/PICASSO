# PICASSO: Pipeline, Code and Analysis for Scalable Single-cell Omics

**Temporarily designed for single-cell learning**
![PICASSO.png](https://raw.githubusercontent.com/Moloch0/PICASSO/main/picture/PICASSO.png)

# Usage

Firstly, clone this repository.

```shell
mkdir ONE_DRIECTORY
cd ONE_DRIECTORY
git clone https://github.com/Moloch0/PICASSO.git
pwd # Gets the absolute path to the folder
DRIECTORY_OF_PICASSO
```

Then, if you want to use R scripts:

```r
source(paste0("DRIECTORY_OF_PICASSO","/PICASSO/picasso.R"), chdir=TRUE)
```

## file management and scripts usage

The file management system and the scripts usage system are synchronized. Here, the second-level header of the document is the pipeline, and the third-level header is the module. `choose_pipeline("scrnaseq")` loads all functions of that pipeline.

```r
> choose_pipeline()
#choose_pipeline(<pipelie>,<module>)

#---- base ----#
--> plot 

#---- bulk_omics ----#
--> bulkrnaseq 
--> chipseq 

#---- other_bioinfo ----#
--> annotation 
--> functional_analysis 
--> sequence_analysis 

#---- sc_omics ----#
--> scrnaseq 
--> smartseq 
--> vdj 

#---- stastics ----#
--> clinical_analysis 
--> deep_learning 
--> experiment_analysis 
--> machine_learning 
Please specify a pipeline and press 'Enter'
scrnaseq # type the pipeline
This job starts at: 2023-04-16 23:35:52 
This script is in: /public/home/luoliheng/PICASSO 
Succeed to load script: scrnaseq
# Or
choose_pipeline("scrnaseq")
```
 


# Acknowledgements

- 2023/4/16
  ***[Cao Jun](caojundudu@qq.com)*** first proposed building a pipeline for better scientific training. ***[Luo Liheng](1351570198@qq.com)*** created the github repository. Along with ***[Zhang Xu](1351570198@qq.com)*** who joined this repository, they were the original contributors.

***Wish you all the best in bioinfomatics research!***
