# TRUST

基于门限 Paillier 加密算法（paililerTD, FastpaiTD），涵盖了 POCF、SOCI、SOCI+ 和 TRUST 等协议的实现。

论文链接：[TRUST: A Toolkit for TEE-Assisted Secure Outsourced Computation over Integers](https://arxiv.org/abs/2412.01073)



## :memo: 文件概览

- #### **local_version**  
  
  - `Baseline`: 小批次梯度下降线性回归的 C++ 实现，明文版本。  
  - `SDTE_C++`: SDTE 协议的 C++ 实现，支持 AES 加密。  
  - `SDTE_SGX`: SDTE 协议的 SGX 版本，结合 AES 加密与 TEE。  
  - `TRUST_C++`: 包含 `SOCI+` 和 `TRUST` 协议的本地可运行版本，相当于  `SOCI+_C++` 和`TRUST_C++` 的实现。  
  
- #### **socket_version**  
  
  - `soci_plus_socket` 



## :alarm_clock: 日志更新

<div align="left">
  <img src="https://img.shields.io/badge/V2-(2024.12.02)-blue" alt="v2" width="200"/>
</div>

- **Operator**: Aptx4869AC
- **Changes**:
  - 提交了 `Baseline_C++`、`SDTE_C++` 与 `SDTE_SGX`，其中包含了 Minibatch Gradient Descent for Linear Regression 的底层实现，明文版本和密文版本（AES-GCM）。
  - 提交 `soci_plus_socket`，其中包含常用协议。该文件也可以看作是socket版 `TRUST_socket` 与 `SOCI+_socket` 的结合。
- **TODO**:
  - 本地版本的 `POCF_C++` 待上传
  - 本地版本的 `SOCI_C++` 待上传
  - socket 版本的 `POCF_C++` 待上传
  - socket 版本的 `SOCI_C++` 待上传
  - SGX 版本的 TRUST，待定
  - SGX 版本的 SEAT，待定

<div align="left">
  <img src="https://img.shields.io/badge/V1-(2024.12.01)-blue" alt="v1" width="200"/>
</div>

- **Operator**: Aptx4869AC  
- **Changes**:  
  - 提交了一个可直接运行的单程序 `TRUST_C++`，其中包含常用协议，方便调试。该文件也可以看作是本地版 `TRUST_C++` 与 `SOCI+_C++` 的结合。



## License

如果学到了，**不要忘记点个Star** :sparkling_heart:
