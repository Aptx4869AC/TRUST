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


#./Plain_Test_Student
#./Plain_Test_Fish


# 运行程序并将输出重定向到临时文件
./Plain_Test_Student > output.log
#./Plain_Test_Fish > output.log
#./Plain_Test_Fish_double > output.log

# 使用 sed 提取 Loss 信息
sed -n 's/.*Loss: \([0-9]*\.[0-9]*\).*/\1/p' output.log > loss_log.txt

# 输出结果
echo "Losses have been logged to build/loss_log.txt"


