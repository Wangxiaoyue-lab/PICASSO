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

# Details

# Comments
`230426.LLH`
## 基本问题
[10x注释建议](https://www.10xgenomics.com/resources/analysis-guides/web-resources-for-cell-type-annotation)

- 不同细胞类型之间的区别是什么？由什么决定？
    细胞间形态，结构，功能的区别
**scRNA-seq取得某一实验条件或时空环境中细胞的部分转录组，希望由此推测细胞类型**

- 需要跨越的生物学障碍
  1. 第一，基因表达的差异在不同细胞类型间是模糊的
  2. 第二，基因表达的差异并不总是转化为细胞功能的差异

- 因此，我们面临的问题是：如何关联细胞表达谱与其类型？有三个层面可以回答它
  - marker(s)
  - 功能相关通路
  - 基因表达模式


## 注释指南： [2021 May / Nat Protoc / Tutorial: guidelines for annotating single-cell transcriptomic maps using automated and manual methods](https://www.nature.com/articles/s41596-021-00534-0)
*首先声明，建议采用手动方法，但是自动方法也可以看一看*

- 推荐一个三步工作流，包括：
    自动注释(只要可能)
    手动注释
    验证

1. 自动注释
使用“marker”(即已知细胞类型中特异表达的基因)或参考数据(即专业注释的单细胞数据)，将它们的基因表达模式(特征)与已知细胞类型的基因表达模式相匹配来识别单个细胞或细胞群


2. 手动注释
研究特定基因及其功能；
验证自动细胞注释并识别新的细胞类型和状态

3. 验证
使用独立的方法确认所选细胞类型的身份和功能

### 自动注释

**注释方法从数据来源上有两类**

1. 基于marker
   1. 从数据库或文献中获得marker与细胞类型之间的关系
   2. 标记特异表达marker基因或基因集的细胞
   3. marker需要获得于特征良好的生物样本

2. 基于参考数据
   1. 将要注释的单细胞数据(查询数据集)与经过专业注释的数据集(参考数据集)进行比较
   2. 将标签从参考细胞转移到查询数据集中相似的细胞
   3. 只有在高质量和相关注释数据集可用的情况下才有可能
3. 特殊方法：整合
   1. 将不同来源的数据集直接合并，消除非生物学因素后降维可视化
   2. 同样依赖高质量数据集，且存在过拟合/批次效应的问题


| 工具 | 类型 | 语言 | 分辨率 | 方法 | 注意 |
| --- | --- | --- | --- | --- | --- |
| singleCell | Net | R | 单细胞 | 相对表达基因对+随机森林 | 比其他方法慢10 - 100倍;精度高 |
| scmap-cluster | 参考数据 | R | 单细胞 | 相关一致性 | 最快的方法;平衡假阳性和假阴性;包括与大型预构建参考数据或自定义参考数据一起使用的web界面 |
| scmap-cell | 参考数据 | R | 单细胞 | 近邻法 | 将单个细胞分配给参考数据中最近的相邻细胞;允许绘制细胞轨迹;快速且可扩展 |
| singleR | 参考数据 | R | 单细胞 | 层次聚类和斯皮尔曼相关性 | 需要大的marker 参考数据;不能扩展到≥10,000个细胞的数据集;包括网页界面与marker数据库 |
| Scikit-learn | 参考数据 | Python | 单细胞 | k近邻，支持向量机，随机森林，最近平均分类器和线性判别分析 | 正确设计分类器，适当训练分类器，避免过度训练 |
| AUCell/ SCENIC | Marker | R | 单细胞 | 曲线下面积估计marker基因集富集程度 | 由于单细胞水平的低检出率，每一种细胞类型都需要许多marker |
| SCINA | Marker | R | 单细胞  |期望最大化，高斯混合模型 | 同时注释群和细胞;robust;包含不正确的marker基因 |
| GSEA/GSVA  |Marker | R/Java | 细胞群 | 富集分析 | marker基因列表必须重新格式化为GMT格式; marker必须在聚类中以相同的方向差异表达 |
| Harmony  |整合 | R  |单细胞  |迭代聚类与调整  |仅集成数据的低维投影;无缝集成到Seurat管道;可能会过拟合 |
| Seurat V4  |整合 | R |单细胞  |MNN锚点+CCA典型相关分析 | 准确性取决于MNN锚点的准确性，这些锚点是跨数据集自动识别的对应细胞 |
| mnnCorrect  |整合  |R  |单细胞  |MNN对+奇异值分解  |准确性取决于MNN对(数据集之间匹配的细胞)的准确性 |
| LIGER  |整合 | R  |单细胞  |非负矩阵分解  |允许解释数据集特定的和共有的变异因素 |

### 手动标注

**当自动方法导致较低的置信度，冲突或缺失的细胞标签时，需要手动注释**
*进行手动注释时，首先应当进行降维可视化（t-SNE/UMAP/PCA）简单查看marker*



- Marker理想来源是来自相关生物体、器官和疾病背景的单细胞图谱

  - 同一组织来源的Bulk RNA-seq
  - 鉴于蛋白质表达可能与mRNA表达
  - 组织染色(免疫组织化学或免疫荧光)
  - 流式细胞术
western blot

- Marker注释的困难：
主宰细胞命运的转录因子通常是比细胞表面蛋白更好的marker
集合来自不同来源的多个marker是困难的
    - 例如，PanglaoDB包含220个B细胞marker，CellMarker包含1426个，但只有66个是一致的

- 细胞群体与marker关系
    - 不表达任何已知细胞类型的marker
    - 表达一种以上细胞类型的marker(双胞体)


- 整体上：差异分析
    - 通过计算群和其他群之间的差异来识别潜在marker(一般方法假阳性率高)
    - 通路/本体富集分析同时对多个功能相关基因进行评分，比单个基因更敏感

- *注意*
注释细胞状态和梯度
稳定的细胞类型具有同质的基因表达，在2D投影图中紧凑
表达梯度表明细胞群中存在的连续差异
    - 细胞周期
    免疫激活
    空间模式
    发育阶段

*必须小心区分生物学上有意义的细胞状态和批次效应，它们会以类似的方式表现出来*

- 亚群
连续差异也许以多峰分布的形式表达，考虑细胞亚群的存在

### 验证
- 其它层次
相似细胞类型之间的区别可能在转录上不可见
染色质状态
DNA甲基化
蛋白

- 实验
mRNA水平只能部分定义细胞类型和功能
新细胞类型必须经过其它手段的补充验证

- 肿瘤样本
SNP
CNV

### 注释总结
|分析阶段|分析角度|注意|建议|
|----:|----:|----:|----:|
|自动|自动方法|快，但对特征差的细胞无效|手动注释|
| |分群注释|错过细胞之间的重要差异|手动细化标签；使用多种基于聚类的方法并比较结果|
| |注释单个细胞|理想状态，但有技术限制|分群注释|
| |基于marker|并非所有类型的细胞都能获取marker；marker冲突或缺失|需要专业知识来获得详细的marker列表|
| |基于参考数据|参考数据不完整或不匹配|使用基于marker的方法|
| | |需要批次校正，降低结果的准确性|分析强生物意义的参考数据；好的实验方案可以抵消批次效应|
| | |参考数据中的错误会延续到结果中|在使用前仔细分析参考数据|
| |比较结果|结果可能不一致|比较各自标签的置信度得分，并考虑标签的一致性(多数原则);使用手动注释解决冲突|
| | | |考虑细胞亚型的可能性，新的细胞类型或连续的细胞状态|
|手动|手动方法|缓慢,劳动密集|只要可能，就从自动注释开始确定细胞标签|
| | |主观|与专家合作;考虑多种细胞类型|
| |基于marker|单个marker无法区分细胞类型|为每种细胞类型使用多个marker|
| | |已知的marker不能区分细胞类型|从文献、额外的实验或专家建议中准备更大的marker列表|
| | |来源之间冲突的marker基因集|选择最能代表数据中生物信号的marker基因集；例如，寻找亚型，使用广泛的基因集|


## Marker 数据库

<https://www.gsea-msigdb.org/gsea/msigdb>
<https://www.haemosphere.org/>
血细胞<https://toppgene.cchmc.org/>
clusterProfiler/Bioconductor


bdbiosciences.com/cdmarkers
手动数据库：singleCellBase<http://cloud.capitalbiotech.com/SingleCellBase/index.jsp>

## 单细胞数据库

<https://www.humancellatlas.org/>
Human Cell Atlas (HCA)
创建所有人类细胞的综合参考图
<https://hubmapconsortium.org/>
Human Biomolecular Atlas Project (HuBMAP)
生成基本的3D组织图谱，构建人体细胞之间的功能和关系图谱

<https://lifetime-fetflagship.eu/欧洲HCA>
<https://www.ebi.ac.uk/gxa/sc/home> 单细胞表达图谱 EMBL-EBI
<https://www.ncbi.nlm.nih.gov/geo/>

<https://www.lungmap.net/肺>
<https://kpmp.org/about-kpmp/肾>
<https://www.gudmap.org/泌尿生殖系统>
<https://humantumoratlas.org/癌症>
<https://tabula-muris.ds.czbiohub.org/小鼠>
...
等等和研究数据集实验方法相似，细胞类型相似的文献数据及各种数据库的数据

## 方法介绍：用于单细胞RNA测序的自动细胞类型识别方法:[2021 / Comput Struct Biotechnol J/Automatic cell type identification methods for single-cell RNA sequencing](https://www.sciencedirect.com/science/article/pii/S2001037021004499)

- 讨论和评估32种已发表的用于scRNA-seq数据分析的自动方法


- 手动方法
输入测试数据集
采用无监督方法聚类，检测差异基因
细胞类型由差异基因中的marker指定

- 自动方法输入：
  - 积极学习和消极学习
训练数据集和测试数据集
  - marker学习
每种细胞类型的marker和测试数据集

- 自动方法实现：
  - 积极学习
分类器
  - 消极学习
近邻细胞
  - marker学习
评分函数


- 积极学习 eager learning
算法判断前先利用训练集数据通过训练得到一个目标函数
在需要进行判断时利用已经训练好的函数进行决策

- 消极学习 lazy learning
不会根据已有的样本创建目标函数
对新进入的样本进行判断的时候才开始分析

- 应用到细胞注释中：
积极学习
收集细胞类型对训练数据集进行分组，然后将测试细胞映射到最近的预注释组
消极学习
基于训练数据集投射细胞以识别最近的邻居细胞，然后根据近邻细胞确定细胞类型
marker学习
利用在给定细胞类型中高度表达的典型细胞marker，使用数学模型分配测试数据集

- 积极学习和消极学习方法
**相同：都需要训练数据集来训练分类器**
**不同：积极学习的训练数据集已经根据细胞类型分组**
**此外：marker学习相较前两种只利用marker信息**

软件列表:

- 积极学习
ACTINN, CaSTLe, CHETAH, clustifyr, Garnett, MarkerCount, MARS, scCapsNet, scClassifR, SciBet, scID, scLearn, scmap-cluster, scMatch, scHPL, scPred, scPretrain, scVI, Seurat V4, SingleCellNet, SingleR，Superscan
- 消极学习
CellAtlasSearch, CELLBLAST, CellFishing.jl, scmap-cell
- marker学习
CellAssign, DigitalCellSorter, markcount, scCATCH, SCINA, SCSA，sctype

文章中并非最新，还有
Cell-ID
TOSICA
Seurat V5等

- 步骤
  - 特征选择（除了scMatch和深度学习方法）
使模型更加准确，缩短训练时间，避免维度诅咒，增强泛化能力
    - 总结
scmap-cell：
高drop-out基因的表现优于高变异基因或随机基因
SciBet：
具有高熵差的基因具有最高的分类精度
  - 预测模型
    - 比较查询数据集细胞与训练数据集的相似性

### 自动方法总结
**积极学习与消极学习优于marker学习的原因？**

基于数据集的学习比基于marker的学习利用更多信息(全基因的表达模式vs几个marker特异表达)
特征选择后的基因比marker实际上更有代表性(高熵差的基因具有最高的分类精度)
预测模型性能不是最关键的因素，数据集质量才是

**使用建议**

当预先标注的数据集可用时，从准确性和速度上看，积极学习和消极学习更可取
不可用时，建议使用基于marker(marker学习)的方法

无监督聚类的工作流假设最小的先验知识，非常适合细胞类型发现
**基于文献的手动方法仍然是金标准**


### 分类单位以及日常术语共识

- 细胞类型 type:一个特定物种的个体在细胞谱系树内和跨细胞谱系树的发育起源和潜能的重复模式，通常反映在共享的分子特性上。

- 细胞状态 state:一种细胞类型内不影响其发育潜力的分子表型的变化(例如，细胞周期，随机波动)。

- 细胞特征/身份 identity:在某一特定时刻仅以其分子表型为特征的单个细胞。

- 细胞系 lineage:一个有机体的细胞间的关系，仅由从一个合子开始的一系列细胞分裂来定义。

- 细胞轨迹 trajectory:从分子表型的相似性推断细胞发育关系的顺序，可能概括发育细胞谱系关系。
