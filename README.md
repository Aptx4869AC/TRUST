# TRUST

基于门限 Paillier 加密算法（paililerTD, FastpaiTD），涵盖了 POCF、SOCI、SOCI+ 和 TRUST 等协议的实现。

论文链接：[TRUST: A Toolkit for TEE-Assisted Secure Outsourced Computation over Integers](https://arxiv.org/abs/2412.01073)



## :memo: 文件概览

- #### **local_version**  
  
  - `Baseline`: mini-batch 梯度下降线性回归的 C++ 实现，底层训练，明文版本
  - `SDTE_C++`: mini-batch 梯度下降线性回归基于 `SDTE` 的 C++ 实现，支持 AES 加密，密文版本 
  - `SDTE_SGX`: mini-batch 梯度下降线性回归基于 `SDTE` 的 SGX 实现，结合 AES 加密与 TEE，密文版本
  - `POCF_C++`: `POCF` 协议的本地 C++ 可运行版本
  - `SOCI_C++`: `SOCI` 协议的本地 C++ 可运行版本
  - `TRUST_C++`: `SOCI+` 和 `TRUST` 协议的本地 C++ 可运行版本，相当于  `SOCI+_C++` 和`TRUST_C++` 的实现 
  
- #### **socket_version**  
  
  - `POCF_socket` :  `POCF` 协议的 socket 可运行版本，需要开多端口
  - `SOCI_socket` :  `SOCI` 协议的 socket 可运行版本，需要开多端口
  - `SOCI+_socket` : 包含 `SOCI+` 和 `TRUST` 协议的 socket 可运行版本，相当于  `SOCI+_socket` 和`TRUST_socket` 的实现，需要开多端口



## :alarm_clock: 日志更新

<div align="left">
  <img src="https://img.shields.io/badge/-TODO-critical" alt="TODO" width="80"/>
</div>

  - SGX 版本的 `TRUST`，上传时间待定
  - SGX 版本的 `SEAT`，上传时间待定


<div align="left">
  <img src="https://img.shields.io/badge/V3-(2024.12.04)-blue" alt="v3" width="200"/>
</div>

- **Changes**:
  - 提交了`SOCI_C++`
  - 提交了`POCF_socket`
  - 提交了`SOCI_socket`

<div align="left">
  <img src="https://img.shields.io/badge/V2-(2024.12.03)-blue" alt="v2" width="200"/>
</div>

- **Changes**:
  - 提交了 `POCF_C++`
  - 提交了 `Baseline_C++`、`SDTE_C++` 与 `SDTE_SGX`，其中包含了 mini-batch 梯度下降线性回归的底层实现，明文版本和密文版本（AES-GCM、TEE）
  - 提交了 `SOCI+_socket`，其中包含常用协议。该文件也可以看作是 `TRUST_socket` 与 `SOCI+_socket` 的结合

<div align="left">
  <img src="https://img.shields.io/badge/V1-(2024.12.02)-blue" alt="v1" width="200"/>
</div>

- **Changes**:  
  - 提交了 `TRUST_C++`，其中包含常用协议。该文件也可以看作是 `TRUST_C++` 与 `SOCI+_C++` 的结合，调试方便



## License

如果学到了，**不要忘记点个Star** :sparkling_heart:
