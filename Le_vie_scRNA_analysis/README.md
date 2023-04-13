# Instructions of PICASSO\scRNA_analysis

**Compiled: 2023-04-04**

**By: Luo Liheng**

## Before all the functions

1. Define `work_dir` and `project_name`
   (Define it first, otherwise you cannot use the `save_file `function)
   like

```R
work_dir <- "/public/home/yujia/analysis/gz+jz/re_analysis_jz/"
project_name <- "MSC"
```

2. Source the script and use all scripts in this PICASSO\scRNA_analysis

```R
source("/public/home/SINGLE/scRNA_analysis/utils.R", chdir = TRUE)
# The output will help you to manage your project
[1] "This scripts is in: /public/home/SINGLE/scRNA_analysis"
[1] "Your work directory is in: /public/home/yujia/analysis/gz+jz/re_analysis_jz/"
[1] "Your project name is: MSC"
```

The script will load `Seurat` and `tidyverse` automatically

## Use the functions provided in utils.R

### `save_file`

1. `save_file `generates a file name and writes data to that file. It is used when saving data, plots, and other outputs.
2. It can be used when saving a single object:
   `save_file(file = "test.rds", data = data, fun = saveRDS)`
3. Or when saving a large data frame to a csv file.
   `save_file(file = "test.csv", data = data, fun = write.csv, row.names = FALSE, col.names = TRUE)`
4. Or when saving a pdf file:
   `save_file(file = "test.pdf", fun = pdf, width = 8, height = 6)`
5. If no file name is provided, the function pastes `work_dir`, `project_name`, and `name_string`, which I often use.

```R
# Arguments
save_file(file = NULL, data = NULL, fun = NULL, name_string = NULL, ...) 
# only saveRDS, write.csv, pdf are supported for now

# Example 1
save_file(fun = pdf, name_string = "markers")
# /public/home/yujia/analysis/gz+jz/re_analysis_jz/MSC_markers.pdf

# Example 2
save_file(data = all_data, fun = saveRDS, name_string = "qc")
# /public/home/yujia/analysis/gz+jz/re_analysis_jz/MSC_qc.RDS

# Example 3
readRDS(save_file(fun = saveRDS, name_string = "qc"))
# Will read the files above
```

### `find_assay`

This function prints the current assay used for the analysis and all assays in the object. The function takes the object as a parameter, which is all_data by default.

```R
# Arguments
find_assay(all_data = all_data) 
```

### `read_refdata`

 This function is a part of `qc_check` which can read files in PICASSO/scRNA_analysis/Refgenome

```R
# Arguments
read_refdata(species, file_type)
# support cell_cycle_markers and annotations in mm and hs for now
```

### `in_data_markers`

This function checks if the genes(vector) passed as argument are present in the dataset.If not, it returns a warning message and fliters them. If all the genes are present, it returns the marker genes.

```R
# Arguments
in_data_markers(genes, dataset)
```

### `time_it`

The function prints the execution time of the wrapped function. The function returns the result of the wrapped function.

```R
# Arguments
time_it(f)
```

### `check_expression`

This function is used to check the expression of certain genes in the dataset.  It takes the dataset, a vector of genes or a dataframe of genes, and returns the average expression of the genes. If a dataframe of genes is provided, the function will return the average expression of each gene and each idents. Also the sum and variance between idents. Other arguments will be passed into `AverageExpression`.

```R
# Arguments
check_expression(all_data = all_data, feats = NULL, ...)

# Example 1
require(readxl)
marker <- read_excel("/public/home/UMAP_markers.xlsx", sheet = 1)
# read markers
all_data<-all_data[,all_data$orig.ident=="a"]
# Find these cells
Idents(all_data)<-all_data$type
# Class them by type
a<-check_expression(all_data,marker)
# check_expression
save_file(data=a,fun=write.csv,name_string="exp_in_types")
```

## Use the functions provided in sc_process.R

### `qc_check`

This function checks the quality of the data. Takes in the all_data Seurat object and a species name, calculates the percent of Mito/Ribo/hemoglobin gene expression, scores the cells for cell cycle phases, normalizes scales the data, finds the variable features, runs PCA and UMAP, finds the doublets by `Find_doublet`  in the data.

```R
# Arguments
qc_check(all_data = all_data, species = "hs",Find_doublet = TRUE)
# Example 1
all_data <- qc_check(all_data, species = "hs", Find_doublet = TRUE)
save_file(data = all_data, fun = saveRDS, name_string = "before_qc")

```

### `qc_check_plot`

The function plots the elbow plot, also the `feats` plot with `VlnPlot` which is grouped by `vln_group`,`DimPlot` group by `dim_group`, `FeaturePlot_scCustom` group by `feats`, colors in FeaturePlot is `colors`

```R
# Arguments
qc_check_plot(all_data, feature_scatter = TRUE, vln_group, dim_group, feats = NULL, colors = pal)

# Example 1
# This function plots QC metrics for each sample, including the number of reads per cell, the number of genes detected per cell,  the percentage of mitochondrial reads per cell, and the UMI vs. gene counts plot. 
pal <- viridis::viridis(n = 10)
save_file(fun = pdf, name_string = "before_qc")
qc_check_plot(all_data,
  feature_scatter = TRUE,
  vln_group = "type",
  dim_group = c("orig.ident", "type", "doublet_info", "Phase"),
  feats = c("nCount_RNA", "nFeature_RNA", "percent_mito", "percent_ribo", "percent_hb"),
  colors = pal
)
dev.off()
```

### `qc_process`

This function runs the QC process on the Seurat object. It takes as input the Seurat object, the number of dimensions to use, the resolutions to use, a boolean indicating whether to run Harmony, and a boolean indicating whether to run SCTransform. It returns the Seurat object after running the QC process.

```R
# Arguments
qc_process(all_data,
    dim_use = 20,
    resolutions = c(0.1, 0.2, 0.3, 0.5),
    run_harmony = TRUE,
    run_sctransform = TRUE,
    group_in_harmony = "orig.ident",
    vars_to_regress = c("percent_mito", "S.Score", "G2M.Score")
)


# Example 1
all_data <- subset(all_data,
        subset = percent_mito < 20 &
            percent_hb < 5 &
            doublet_info == "Singlet" &
            nFeature_RNA > 500 &
            nFeature_RNA < 7500
    ) %>% qc_process()
# Note, if you choose to run_sctransform, the default assay will be "SCT", and the feature will be nFeature_SCT
```

### `qc_process_plot`

This code is used to create a QC report for the data. `DimPlot` group by `dim_group`, `FeaturePlot_scCustom` group by `feats`, colors in FeaturePlot is `colors`. It also includes a plot of the clusters for each of the four resolutions.

```R
# Arguments
qc_process_plot(
    all_data,
    dim_group = c("orig.ident", "type", "Phase"),
    feats = c("nFeature_RNA", "nCount_RNA", "percent_mito", "percent_hb"),
    colors = pal,
    resolutions = c(0.1, 0.2, 0.3, 0.5))
# Note the default ident will be the smallest resolution

# Example 1
save_file(fun = pdf, name_string = "after_qc")
pal <- viridis::viridis(n = 10)
qc_process_plot(all_data)
dev.off()
```

### `find_markers`

The function is to find marker genes. Determine whether `all `markers are found. Choose `ident `to use. And `loop_var `could be used to find markers for each cluster.  Additional arguments will be passed to `FindMarkers` or `FindAllMarkers` .

```R
# Arguments
find_markers(all_data, all = TRUE, ident = NULL, loop_var = NULL, species,...)
# Note that the output will be different from the original result
# Percentage difference, type of marker gene(up or down),gene description will be added

# Example 1
marker_genes <- find_markers(all_data,
  all = FALSE, ident = NULL, loop_var = 0:3, species="mm",
  ident.1 = "ident_A", group.by = "type", only.pos = FALSE
)
save_file(data = marker_genes, fun = write.csv, name_string = "markers", row.names = FALSE, col.names = TRUE)
# all = FALSE. Use FindMarkers
# ident = NULL. Otherwise, this function will be executed
# all_data <- SetIdent(all_data, value = ident)
# loop_var = 0:3. Check your defult ident. 
# Here I want to loop from  cluster 0 to cluster 3
# species will be used to add gene description by gene name.
# group.by="type". In a cluster, I want to compare columns with group name "type"
# ident.1 = "ident_A". In this group, I want to compare "ident_A" with other "type"s
```

### `plot_makers`

1. This equation takes as input a marker vector or a data frame of markers with column names of the specified name. Also specify that the name will be the upper-left identifier for all graphics.
2. Optionally, output dotplot or/and featureplot. `max_markers `identifies the maximum number of genes drawn in a dotplot, and markers for more than one dotplot.
   The cluster accepts numbers as clustering numbers; other forms do not output `Clustered_DotPlot`(note that this graph uses data from scale.data).
3. `col_dot `and `col_fea `at the end specify the colors for dotplot and featureplot. The remaining arguments are passed to the first function of the plot. Only dotplot accepts arguments if dotplot and featureplot are plotted together. Featureplot has only 95 percentiles by default.
4. This function uses `in_data_markers` to test the entered gene.

```R
# Arguments
plot_makers(unfilterd_markers,
  all_data, dot_plot = TRUE, max_markers = 40, cluster = 10,
  feature_plot = TRUE, col_dot = viridis_plasma_dark_high, col_fea = pal, ...
)

# Example 1
require(readxl)
marker <- read_excel("/public/home/files/markers.xlsx", sheet = 1)
pal <- viridis::viridis(n = 10)

save_file(fun = pdf, name_string = "markers")
plot_makers(marker, all_data,
  dot_plot = TRUE, max_markers = 40, cluster = "no",
  feature_plot = TRUE, col_dot = viridis_plasma_dark_high, col_fea = pal
)
dev.off()
```

### `time_it`

 The function prints the execution time of the wrapped function. The function returns the result of the wrapped function.

```R
# Arguments


```

### `time_it`

 The function prints the execution time of the wrapped function. The function returns the result of the wrapped function.

```R
# Arguments



```

## Please try to install them, follow their tutorials and their citation rules

- [Seurat](https://satijalab.org/seurat/)
- [tidyverse](https://www.tidyverse.org/)
- [scCustomize](https://samuel-marsh.github.io/scCustomize/)
- [DoubletFinder](https://github.com/chris-mcginnis-ucsf/DoubletFinder)
- [clustree](https://lazappi.github.io/clustree/)
