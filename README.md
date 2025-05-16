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
  
- #### **simple_version**

  - `TRUST_C++`: `TRUST` 协议的本地 C++ 可运行版本，无对比方案

- #### **window_version**

  - `Final_TRUST_socket` : 主动权在`CP.cpp`
  - `Final_TRUST_socket_v2` : 主动权在 `Client.cpp`

## :alarm_clock: TODO

  - SGX 版本的 `TRUST_SGX`，待论文录用上传
  - SGX 版本的 `SEAT_SGX`，待论文录用上传



## License

如果学到了，**不要忘记点个Star** :sparkling_heart:
