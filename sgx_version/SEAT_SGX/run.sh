#!/bin/bash
# Check if ./app exists and then clean up
if [ -f "./app" ]; then
    make clean
fi

# Build the project
make

#./app

# 定义不同的 epoch 数值
epochs=(200)

# 遍历每个 epoch 值
for epoch in "${epochs[@]}"
do
    # 更新 C++ 程序中的 epoch 参数
    # 假设你的 C++ 程序接受 epoch 作为命令行参数
    ./app $epoch > output_$epoch.log  # 输出到不同的文件

    # 使用 sed 提取 Loss 信息
    sed -n 's/.*Loss: \([0-9]*\.[0-9]*\|[0-9]*\).*/\1/p' output_$epoch.log > loss_log_$epoch.txt

    # 输出结果
    echo "Epoch $epoch: Losses have been logged to loss_log_$epoch.txt"
done


# Check if ./app exists and then clean up
if [ -f "./app" ]; then
    make clean
fi


