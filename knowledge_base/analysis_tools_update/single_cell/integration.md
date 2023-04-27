# Methods

| value | tool | time | platform | link | description |  |
| ----- | ---- | ---- | -------- | ---- | ----------- | - |
|       |      |      |          |      |             |  |
|       |      |      |          |      |             |  |
|       |      |      |          |      |             |  |
|       |      |      |          |      |             |  |
|       |      |      |          |      |             |  |
|       |      |      |          |      |             |  |
|       |      |      |          |      |             |  |

# Benchmark

# Benchmark

Benchmark论文的一点看法，Benchmark有很多方面，不存在一个客观的评估，只存在一个某方面有一定参考价值的评估。
[**2020 Jan** / *Genome Biol* / A benchmark of batch-effect correction methods for single-cell RNA sequencing data](https://pubmed.ncbi.nlm.nih.gov/31948481/) # 单模态

[**2022 Jun** / *Nat Methods* / Benchmarking spatial and single-cell transcriptomics integration methods for transcript distribution prediction and cell type deconvolution](https://www.nature.com/articles/s41592-022-01480-9) # 多模态
[https://github.com/QuKunLab/SpatialBenchmarking](https://github.com/QuKunLab/SpatialBenchmarking)

[**2022 Jan** / *Nat Methods* /Benchmarking atlas-level data integration in single-cell genomics](https://www.nature.com/articles/s41592-021-01336-8) # 单模态但针对百万级细胞，Seurat V5没有参与评估但有可能效果不错
[https://github.com/theislab/scib-reproducibility](https://github.com/theislab/scib-reproducibility)
[https://theislab.github.io/scib-reproducibility/](https://theislab.github.io/scib-reproducibility/)
[https://github.com/theislab/scib-pipeline](https://github.com/theislab/scib-pipeline)

# Details

## GLUE

[**2022 Oct** / *Nat Biotechnol* /Multi-omics single-cell data integration and regulatory inference with graph-linked embedding](https://www.nature.com/articles/s41587-022-01284-4)
[https://github.com/gao-lab/GLUE](https://github.com/gao-lab/GLUE) # 多模态整合结合机器学习的先进技术

# Comments

`230426.LLH` : 本文主要基于SCALEX论文： [2022 Oct / Nat Commun / Online single-cell data integration through projecting heterogeneous datasets into a common cell-embedding space](https://www.nature.com/articles/s41467-022-33758-z)

## 整合分类

整合有许多方面：

- 批次：
  同一样本的生物学/技术重复
  同一样本在不同时间点处理
  不同建库策略，10X平台，Drop/SMART2
  不同测序平台，BGI/Illumina
  不同样本
  同一样本的不同因素处理
- 技术：
  RNA ATAC CHIP ...
- 模态：
  CITE SHARE ASAP ...

主要应该分为多模态（数据结构差异巨大）和单模态（数据结构单一）

## 整合目的

将不同来源的数据汇总，我们就可以更好地分辨其同质性和异质性

- 计算方便（多个组或者不同来源数据）
- 可视化方便

## 整合思想

一个数据样本或者多个数据样本的，其构成其内部异质性的因素有两种，*生物学因素*，即我们关注的方面（基因表达，处理，模态）；或者*非生物学因素*，我们想要舍弃的方面（批次，物种，平台）。
舍弃非生物学因素，这不就是**质控**吗？
那么我们有任何先验知识来确定到底哪些是生物学因素或者非生物学因素吗？
***没有。***
**使用技术就要受到技术的影响**，我们可以认定某些因素的表达是有规律的，如基因表达的零膨胀负二项分布以及dropout出现的泊松分布，凡是不符合这个分布的，我们认为其存在一定问题。但是还是要根据具体情况判断。

## 方法分类

两类方法并不是完全独立的，以及这种分类并不一定科学。Andrews T S. et al. Nat Protoc. 2021.

- 第一类：以细胞标签为中心（部分）

  - 寻找基因表达模式相近的细胞（锚点），并归类为同种细胞
  - 利于细胞注释
  - MNN, Seurat, Harmony, BBKNN, Scanorama
- 第二类：跨数据集整合（整体）

  - 应用一定的计算方法混合样本
  - 具体结果和使用的算法相关
  - online iNMF, LIGER, SCALEX

### 基于KNN（MNN，Seurat，BBKNN, Scanorama, ...）

- MNN假设：Haghverdi. et al. Nat Biotechnol. 2018.

  1. 两个批次中至少存在同一个细胞群
  2. 批次效应几乎与生物效应正交
  3. 批次效应变化远小于不同细胞类型之间的生物效应变化
- MNN流程：

  1. 整体归一化
  2. 计算细胞间的欧氏距离
  3. 对于单个细胞，寻找另一批次的k个近邻细胞
  4. 重复此操作于每个细胞
  5. 相互近邻细胞被标记为同一类型细胞
- Harmony流程：Korsunsky I. et al. Nat Methods. 2019.

  1. 采用PCA方法并进行模糊聚类（将细胞尽可能分配给多个类群）-> 保留整体拓扑结构
  2. 用每群细胞的质心计算矫正因子 -> 质心对应潜在细胞类型和细胞状态
  3. 校正细胞并重新分配聚类，循环直至收敛

### 第二类：LIGER，iNMD

- PCA：线性降维，用一个低维矩阵表示高维矩阵的大部分信息
- 非负矩阵分解（NMD）：非线性降维，要求分量为非负值 Welch JD. et al. Cell. 2019.
- iNMD分解矩阵的因子，相当于PCA降维之后的PC轴，但iNMD的因子比 PCA的PC轴可解释性更强
- 根据最大因子loading为每个细胞分配一个标签，使用KNN替代原有标签，最后连接具有相似表达模式的细胞
- SCALEX, GLUE（没有用过SCALEX，文章影响因子也不是很高，了解一下即可）
- SCALEX的基本设计概念：

  1. 实现一个广义投影函数，分离批次处理相关组件和批次处理无关组件
  2. 将批次处理无关组件投射到公共嵌入空间中
- 设计过程：
- 核心元件：variational autoencoder (VAE)变分自编码器

  1. 从输入单细胞数据中提取生物相关的潜在特征
  2. 合并批次处理信息
  3. 重构原始数据
  4. Domain-Specific Batch Normalization(DSBN):在单细胞数据重构过程中处理批次特异性特征
     **两种特征的分离是在线整合能够实现的基础**

## 方法总结

- **第一类：寻找基因表达模式相近的细胞，配对后进行后续分析（每加入一个新数据集就需要重新匹配锚点）**
- **第二类：使用算法分离细胞差异，使用主要差异进行后续分析（部分方法可以运算后直接将细胞投入到公共空间）**

文中也使用了一些评估方法和公式

$$
ARI = \frac{\sum_{ij}\binom{n_{ij}}{2} - [\sum_i\binom{a_i}{2}\sum_j\binom{b_j}{2}]/\binom{n}{2}}{\frac{1}{2}[\sum_i\binom{a_i}{2}+\sum_j\binom{b_j}{2}] - [\sum_i\binom{a_i}{2}\sum_j\binom{b_j}{2}]/\binom{n}{2}}
$$

Adjusted Rand Index (ARI)：聚类分布之间的相似度得分（基于随机分布）
其中，nij表示第i个聚类中属于第j个类别的元素数量，ai表示第i个聚类中的元素数量，bj表示第j个类别中的元素数量，n表示元素总数。

$$
NMI(U, V) = \frac{2I(U, V)}{H(U) + H(V)}
$$

Normalized Mutual Information (NMI)：预测和真实聚类的分布（基于信息熵）
其中，U和V分别表示两个聚类结果，I(U,V)表示U和V之间的互信息，H(U)和H(V)分别表示U和V的熵。

## 整合总结

***整合注意事项：***

- 整合不能改变细胞类型，只是为了方便可视化以及后续分析
- 整合没有统一的评价标准，直观的评价是UMAP
- 整合需要考虑整个分析流程的兼容性

***整合方法根据整合目的来确定：***

- 简单批次整合 -> Seurat自带的批次矫正函数**SCT/Harmony**
- 复杂批次整合(不同来源/平台) -> Seurat V3 CCA+MNN/Harmony
- 多模态数据整合 -> Seurat V4/相应文章的流程分析方法
- 大型样本整合(250K) -> GLUE/SCALEX

***值得注意的是，Seurat V5可能吊打以上所有。详情值得单独一章。***
