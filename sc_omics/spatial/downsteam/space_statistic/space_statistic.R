# 0 创建对象
## 创建坐标矩阵
library(spdep)
coords <- matrix(c(
    1, 1,
    2, 1,
    1, 2,
    2, 2
), nrow = 4, ncol = 2, byrow = TRUE)
rownames(coords) <- c("A", "B", "C", "D")
coords

## 创建邻接矩阵
### 创建nb对象


### 创建listw对象



## 坐标矩阵转化为邻接矩阵
# 根据坐标矩阵和距离阈值创建邻接矩阵
nb <- spdep::dnearneigh(coords, 0, 1.5) # 第2、3个参数给出了判定为相邻的距离区间
nb2listw(nb)$neighbours


## 创建空间权重矩阵
listw <- nb2listw(your_neighbours, style = "W")
