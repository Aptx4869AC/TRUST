#!/bin/bash

# 编译 CP.cpp
g++ -o CP CP.cpp -lgmp -fopenmp

# 检查编译是否成功
if [ $? -ne 0 ]; then
    echo "编译 CP.cpp 失败"
    exit 1
fi

echo "CP.cpp 编译成功"

# 编译 CSP.cpp
g++ -o CSP CSP.cpp -lgmp -fopenmp

# 检查编译是否成功
if [ $? -ne 0 ]; then
    echo "编译 CSP.cpp 失败"
    exit 1
fi

echo "CSP.cpp 编译成功"

# 编译 Client.cpp
g++ -o Client Client.cpp -lgmp -fopenmp

# 检查编译是否成功
if [ $? -ne 0 ]; then
    echo "编译 Client.cpp 失败"
    exit 1
fi

echo "Client.cpp 编译成功"

echo "所有程序编译完成"
