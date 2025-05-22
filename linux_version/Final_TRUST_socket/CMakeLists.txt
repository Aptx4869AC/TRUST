cmake_minimum_required(VERSION 3.16.3)
project(TRUST)  
 

# 包含头文件目录
include_directories("${PROJECT_SOURCE_DIR}/include")

# 添加三个可执行文件
add_executable(CP CP.cpp) 
add_executable(CSP CSP.cpp) 
add_executable(Client Client.cpp) 
 
# 为每个目标链接OpenMP和GMP
foreach(target CP CSP Client)
    target_link_libraries(${target}
            -lgmp
            -lpthread
            -fopenmp
    )
endforeach()
