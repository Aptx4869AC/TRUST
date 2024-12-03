cmake_minimum_required(VERSION 3.0)
project(MyProject)

# 添加可执行文件并指定源文件
add_executable(Plain_Test_Fish Plain_Test_Fish.cpp)
add_executable(Plain_Test_Student Plain_Test_Student.cpp)

# 设置链接的库
target_link_libraries(Plain_Test_Fish -lgmp -lpthread -fopenmp)
target_link_libraries(Plain_Test_Student -lgmp -lpthread -fopenmp)
