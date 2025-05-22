cmake_minimum_required(VERSION 3.27.8) 
project(TRUST)  
 

# 包含头文件目录
include_directories("${PROJECT_SOURCE_DIR}/include")
# 链接库目录（如果libgmp.a或libgmp.dll 在自定义路径中）
link_directories("${PROJECT_SOURCE_DIR}/lib")
 
# 查找OpenMP
find_package(OpenMP REQUIRED)

# 查找并链接pthread
find_package(Threads REQUIRED) 
 
# 添加三个可执行文件
add_executable(CP CP.cpp) 
add_executable(CSP CSP.cpp) 
add_executable(Client Client.cpp) 
 
# 为每个目标链接OpenMP和GMP
foreach(target CP CSP Client)
    target_link_libraries(${target} 
        OpenMP::OpenMP_CXX 
        gmp
        ws2_32
        Threads::Threads
    )
endforeach()
 

# mingw32-make