#!/bin/bash

# 1 填写yaml必填项
## 包括name、aim和discussion
## 由于result部分的尽量自动化。所以结果是尽量标号、对结果的讨论则位于discussion部分

# 2 yaml生成html预览

# 3 html加入input部分description.yaml的数据描述作为data部分

# 4 html根据log.out里的日志总结使用的包或者软件、提取报错部分、提取统计部分

# 5 html加入将script输给gpt模型的总结返回

# 6 html加入output下文件情况与sha256值

# 7 picture下如果存在pdf图片则转为像素图并添加水印

# 8 html加入output下水印版pdf图片并且编号A、B、C、D

# 9 html转pdf，如若每页过长则换页
