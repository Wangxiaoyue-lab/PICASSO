# 0 创建对象
## 创建坐标矩阵
library(spdep)
coords_matrix <- matrix(c(
    1, 1,
    2, 1,
    1, 2,
    2, 2
), nrow = 4, ncol = 2, byrow = TRUE)
rownames(coords_matrix) <- c("A", "B", "C", "D")
coords_matrix

## 创建邻接矩阵
### 获得listw对象(基于list)
listw <- mat2listw(neighbours_matrix) # matrix→listw
listw <- nb2listw(nb) # nb→listw

### 创建nb对象(基于list)
# 根据坐标矩阵和距离阈值创建邻接矩阵
nb <- spdep::dnearneigh(coords_matrix, 0, 1.5) # 第2、3个参数给出了判定为相邻的距离区间

nb <- mat2listw(neighbours_matrix)$neighbours # listw→nb


## 坐标矩阵转化为邻接矩阵
nb2listw(nb)$neighbours


## 创建空间权重矩阵
listw <- nb2listw(your_neighbours, style = "W")
