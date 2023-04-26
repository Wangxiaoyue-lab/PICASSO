# Methods

| value | tool             | time | platform | link | description |  |
| ----- | ---------------- | ---- | -------- | ---- | ----------- | - |
|       | DoubletDetection | 2017 | python   |      |             |  |
|       |                  |      |          |      |             |  |
| **    | scrublet         | 2018 | Python   |      |             |  |
|       | DobuletDecon     | 2018 | R        |      |             |  |
| ****  | DoubletFinder    | 2019 | R        |      |             |  |
| ***   | cxds             | 2019 | R        |      |             |  |
|       | chord            | 2021 | R        |      |             |  |

# Benchmark

Benchmarking Computational Doublet-DetectionMethods for Single-Cell RNA Sequencing Data

https://www.ncbi.nlm.nih.gov/pmc/articles/PMC7897250/

DoubletFinder最准

cxds最快

 **doubletCells** : The method was executed by following the instruction at [https://bioconductor.statistik.tu-dortmund.de/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-6-doublet.html](https://bioconductor.statistik.tu-dortmund.de/packages/3.8/workflows/vignettes/simpleSingleCell/inst/doc/work-6-doublet.html). Doublet scores were obtained from the *dblCells* function in R package *scran* (v 1.16.0) with parameters set to default.

 **Scrublet** : R package *reticulate* (v 1.16) was used to execute the python module *scrublet* (v 0.2.1). The parameters were set by following the instruction at [https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb](https://github.com/AllonKleinLab/scrublet/blob/master/examples/scrublet_basics.ipynb). Doublet scores were obtained from the function  *Scrublet.scrub_doublets* .

 **cxds** , **bcds** and  **hybrid** : These three methods were executed by following the instructions at [https://github.com/kostkalab/scds](https://github.com/kostkalab/scds). Doublet scores were obtained from the functions *cxds, bcds* and *cxds_bcds_hybrid* in R package *scds* (v 1.2.0) with parameters set to default.

 **DoubletDetection** : R package *reticulate* (v 1.16) was used to execute the python module  *doubletdetection* . The parameters were set by following the instruction at [https://nbviewer.jupyter.org/github/JonathanShor/DoubletDetection/blob/master/tests/notebooks/PBMC_8k_vignette.ipynb](https://nbviewer.jupyter.org/github/JonathanShor/DoubletDetection/blob/master/tests/notebooks/PBMC_8k_vignette.ipynb). The parameter *n_iters* was set to 5, as larger values were found to increase the running time significantly, but with little improvement in performance. Doublet scores were obtained from the function  *doubletdetection.BoostClassifier.fit* .

 **DoubletFinder** : The method was executed by following the instruction at [https://github.com/chris-mcginnis-ucsf/DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder). Doublet scores were obtained from the function *doubletFinder_v3* in R package *DoubletFinder* (v 2.0.3) with parameters set to default.

 **DoubletDecon** : The method was executed by following the instruction at [https://github.com/EDePasquale/DoubletDecon](https://github.com/EDePasquale/DoubletDecon). Doublet predictions were obtained from the function *Main_Doublet_Decon* in R package *DoubletDecon* (v 1.1.5) with parameters set to default.

 **Solo** : The method was executed by following the instruction at the GitHub repository [https://github.com/calico/Solo](https://github.com/calico/Solo). Every scRNA-seq count matrix was transformed into the *loom* format as required by the method. The parameters were set the same as those in the file  *Solo_params_example.json* , which was downloaded from the GitHub repository. Doublet scores were obtained from the file  *softmax_scores.npy* .

# Details

## DoubletFinder


## cxds



# Comments
