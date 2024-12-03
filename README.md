# Changelog
本库基于门限 Paillier 加密算法（paililerTD, FastpaiTD），涵盖了 POCF、SOCI、SOCI+ 和 TRUST 等协议的实现。


## [v2] - 2024-12-02
- **Operator**: Aptx4869AC
- **Changes**:
  - 提交了 `Baseline_c++`、`SDTE_c++` 与 `SDTE_SGX` ，其中包含了 Minibatch Gradient Descent for Linear Regression 的底层实现，明文版本和密文版本（AES-GCM）。
- **TODO**:
  - 本地版本的 `POCF_c++` 待上传
  - 本地版本的  `SOCI_c++` 待上传
  - socket 版本的 `POCF_c++` 待上传
  - socket 版本的 `SOCI_c++` 待上传
  - socket 版本的 `SOCI+_c++` 待上传
  - SGX 版本的 TRUST，待定
  - SGX 版本的 SEAT，待定


## [v1] - 2024-12-01
- **Operator**: Aptx4869AC  
- **Changes**:  
  - 提交了一个可直接运行的单程序 `TRUST_c++`，其中包含常用协议，方便调试。该文件也可以看作是本地版 `TRUST_c++` 与 `SOCI_c++` 的结合。
- **TODO**:  
  - 本地版本的 `POCF_c++` 待上传
  - 本地版本的  `SOCI_c++` 待上传  
  - socket 版本的 `POCF_c++` 待上传  
  - socket 版本的 `SOCI_c++` 待上传 
  - socket 版本的 `SOCI+_c++` 待上传
  - SGX 版本的 TRUST 与 SEAT，等论文公开再上传
