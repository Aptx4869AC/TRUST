# TRUST: A Toolkit for TEE-Assisted Secure Outsourced Computation over Integers

基于门限 Paillier 加密算法（paililerTD, FastpaiTD），涵盖了论文 [POCF](https://ieeexplore.ieee.org/abstract/document/7500106)、[SOCI](https://ieeexplore.ieee.org/abstract/document/9908577)、[SOCI+](https://ieeexplore.ieee.org/abstract/document/10531248) 和 [TRUST](https://arxiv.org/abs/2412.01073) 等协议的实现。


## :memo: 文件概览

- #### **local_version**  
  
  - `Baseline`: mini-batch 梯度下降线性回归的 C++ 实现，底层训练，明文版本
  - `SDTE_C++`: mini-batch 梯度下降线性回归基于 `SDTE` 的 C++ 实现，支持 AES 加密，密文版本 
  - `SDTE_SGX`: mini-batch 梯度下降线性回归基于 `SDTE` 的 SGX 实现，结合 AES 加密与 TEE，密文版本
  - `SEAT_C++`: mini-batch 梯度下降线性回归基于 `SEAT` 的 C++ 实现，支持 FastPai 加密，密文版本
  - `POCF_C++`: `POCF` 协议的本地 C++ 可运行版本
  - `SOCI_C++`: `SOCI` 协议的本地 C++ 可运行版本
  - `TRUST_C++`: `SOCI+` 和 `TRUST` 协议的本地 C++ 可运行版本，相当于  `SOCI+_C++` 和`TRUST_C++` 的实现，涵盖对比方案
  
- #### **socket_version**  
  
  - `POCF_socket` :  `POCF` 协议的 socket 可运行版本，需要开多端口
  - `SOCI_socket` :  `SOCI` 协议的 socket 可运行版本，需要开多端口
  - `SOCI+_socket` : 包含 `SOCI+` 和 `TRUST` 协议的 socket 可运行版本，相当于  `SOCI+_socket` 和`TRUST_socket` 的实现，涵盖对比方案，需要开多端口

- #### **sgx_version**

  - `TRUST_SGX` : `TRUST` 的 SGX 实现
  - `SEAT_SGX` : `SEAT` 的 SGX 实现

- #### **diy_version**

  - 做任何你想做的操作（e.g., 浮点数，负数）

- #### **simple_version**

  - `TRUST_C++`: `TRUST` 协议的本地 C++ 可运行版本，无对比方案

- #### **linux_version**

  - `Final_TRUST_socket` : 主动权在`CP.cpp`

- #### **window_version**

  - `Final_TRUST_socket` : 主动权在`CP.cpp`
  - `Final_TRUST_socket_v2` : 主动权在 `Client.cpp`

  
## :pencil2: Remark

Windows和Linux平台的Socket实现存在显著性能差异。相同代码（`Final_TRUST_socket`）在Windows上效率接近串行版本，而在Linux上则较慢。
这主要源于：

  - 内核网络栈实现差异
  - 进程调度机制不同
  - TCP/IP协议栈在不同系统实现细节不同
  - 内存管理和缓冲区优化程度不同



### Window

```shell
[epoch = 80] --- m_1 = 363, m_2 = 19, m_3 = 217
[epoch = 81] --- m_1 = -613, m_2 = -297, m_3 = 254
[epoch = 82] --- m_1 = 1000, m_2 = 383, m_3 = 203
[epoch = 83] --- m_1 = -999, m_2 = -455, m_3 = 252
[epoch = 84] --- m_1 = 292, m_2 = 317, m_3 = 204
[epoch = 85] --- m_1 = 416, m_2 = -862, m_3 = 201
[epoch = 86] --- m_1 = 50, m_2 = 556, m_3 = 180
[epoch = 87] --- m_1 = -595, m_2 = 524, m_3 = 158
[epoch = 88] --- m_1 = 974, m_2 = 612, m_3 = 142
[epoch = 89] --- m_1 = -21, m_2 = 100, m_3 = 217
[epoch = 90] --- m_1 = -255, m_2 = 692, m_3 = 243
[epoch = 91] --- m_1 = -725, m_2 = 746, m_3 = 135
[epoch = 92] --- m_1 = 547, m_2 = -392, m_3 = 205
[epoch = 93] --- m_1 = 504, m_2 = -277, m_3 = 165
[epoch = 94] --- m_1 = 187, m_2 = -350, m_3 = 150
[epoch = 95] --- m_1 = -596, m_2 = 696, m_3 = 182
[epoch = 96] --- m_1 = -246, m_2 = 498, m_3 = 202
[epoch = 97] --- m_1 = -542, m_2 = -165, m_3 = 153
[epoch = 98] --- m_1 = 555, m_2 = 378, m_3 = 225
[epoch = 99] --- m_1 = 900, m_2 = 798, m_3 = 195
[epoch = 100] --- m_1 = -10, m_2 = -106, m_3 = 140
Trust_FMUL (100 average) time is  ------  14.620004 ms
Trust_FCMP (100 average) time is  ------  14.590008 ms      
Trust_FEQL (100 average) time is  ------  29.439998 ms      
Trust_FABS (100 average) time is  ------  14.040010 ms      
Trust_FTRN (100 average) time is  ------  27.969992 ms
```

### Linux

```shell
[epoch = 80] --- m_1 = 363, m_2 = 19, m_3 = 217
[epoch = 81] --- m_1 = -613, m_2 = -297, m_3 = 254
[epoch = 82] --- m_1 = 1000, m_2 = 383, m_3 = 203
[epoch = 83] --- m_1 = -999, m_2 = -455, m_3 = 252
[epoch = 84] --- m_1 = 292, m_2 = 317, m_3 = 204
[epoch = 85] --- m_1 = 416, m_2 = -862, m_3 = 201
[epoch = 86] --- m_1 = 50, m_2 = 556, m_3 = 180
[epoch = 87] --- m_1 = -595, m_2 = 524, m_3 = 158
[epoch = 88] --- m_1 = 974, m_2 = 612, m_3 = 142
[epoch = 89] --- m_1 = -21, m_2 = 100, m_3 = 217
[epoch = 90] --- m_1 = -255, m_2 = 692, m_3 = 243
[epoch = 91] --- m_1 = -725, m_2 = 746, m_3 = 135
[epoch = 92] --- m_1 = 547, m_2 = -392, m_3 = 205
[epoch = 93] --- m_1 = 504, m_2 = -277, m_3 = 165
[epoch = 94] --- m_1 = 187, m_2 = -350, m_3 = 150
[epoch = 95] --- m_1 = -596, m_2 = 696, m_3 = 182
[epoch = 96] --- m_1 = -246, m_2 = 498, m_3 = 202
[epoch = 97] --- m_1 = -542, m_2 = -165, m_3 = 153
[epoch = 98] --- m_1 = 555, m_2 = 378, m_3 = 225
[epoch = 99] --- m_1 = 900, m_2 = 798, m_3 = 195
[epoch = 100] --- m_1 = -10, m_2 = -106, m_3 = 140
Trust_FMUL (100 average) time is  ------  25.115513 ms
Trust_FCMP (100 average) time is  ------  25.743598 ms
Trust_FEQL (100 average) time is  ------  38.334907 ms
Trust_FABS (100 average) time is  ------  24.601104 ms
Trust_FTRN (100 average) time is  ------  36.389111 ms
```



## License

如果学到了，**不要忘记点个Star** :sparkling_heart:

Copyright :copyright:2024 [Aptx4869AC](https://github.com/Aptx4869AC)
