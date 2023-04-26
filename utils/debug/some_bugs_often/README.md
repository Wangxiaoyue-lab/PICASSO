# install

## monocle

monocle 安装容易失败，因为缺少软件icu里的库文件

解决方案是安装icu，将里面的lib文件夹添加到LD_LIBRARY_PATH里去


## Rcistarget

依赖的arrow包很难安，可能出现报错显示undefined symbol error

最好不安CRAN版本

按这个

https://github.com/apache/arrow/blob/main/r/R/install-arrow.R

```R
 install_arrow(nightly=T)
```
