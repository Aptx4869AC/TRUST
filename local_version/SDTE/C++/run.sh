#!/bin/bash

# 判断 build 目录是否存在
if [ -d "build" ]; then
  # 如果存在，删除整个目录
  rm -rf build
fi

mkdir build
cd build

# 运行 cmake 命令生成 Makefile
cmake ..

# 编译项目
cmake --build .

./Cipher_Test


