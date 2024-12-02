cmake_minimum_required(VERSION 3.0)
project(MyProject)

# 添加可执行文件并指定源文件
add_executable(main Test.cpp)

# 设置链接的库
target_link_libraries(main -lgmp -lpthread -fopenmp)
