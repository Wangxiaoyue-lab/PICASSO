# download

Firstly, please download the software by the url

```shell
cd ~/software
software=X
mkdir ${software} && cd ${software}

```

# dompression


```shell
tar -zxf *
```


# generate makefile 

两种方式都可以生成makefile文件，不同软件选择的方式可能不一样

```shell
./configure --prefix=<>
```

来自于autoconf



```shell
module load compiler/cmake/3.20.1
cd ${software}
mkdir bulid && cd build
cmake ..
```

比configue新


# compile and install

```shell
make && make install
```


并行加速

```
make -j4
```
