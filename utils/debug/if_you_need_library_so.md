# 1 you should check the directionary of two local library firstly

```
/lib64
/usr/lib64
```

# 2 if you don't find it, you can search it

maybe you can insult gpt like new bing

# 3 choose the suitable version and download

sometimes the lib.so belongs to a complex software, so you should isntall the software which should be specified the version.
And you can achieve the lib.so file in the `bin` or `lib` folder.
What's more, maybe you can copy it from other's account.

# 4 install it in the software path

```
~/software
```

# 5 `cp` and `ln -s` the key `.so` file to the lib path

my lib path is

```
~/R/lib_add
```

# 6 modificate the ~/.bashrc file

add the code

```
export LD_LIBRARY_PATH=~/R/lib_add:$LD_LIBRARY_PATH
```

please notice that the code can search the sub directionary

```shell
source ~/.bashrc
```

# 7 some important so

monole →libicuuc.so.58→ icu

SCopeLoomR→libhdf5_hl.so.100→hdf5_1.10.6

rgdal→proj → https://proj.org/download.html

rgdal→gdal
