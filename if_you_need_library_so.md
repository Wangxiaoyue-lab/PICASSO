# 1 you should check the directionary of two local library firstly
```
/lib64
/usr/lib64
```
# 2 if you don't find it, you can search it
maybe you can insult gpt like new bing

# 3 choose the suitable version and download

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
export LD_LIBRARY_PATH=/public/home/caojun/R/lib_add:$LD_LIBRARY_PATH
```

 