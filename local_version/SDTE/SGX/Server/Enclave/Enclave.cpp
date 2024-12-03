/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-02
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include <stdio.h>
#include <stdarg.h>
#include <sgx_trts.h>
#include "Enclave.h"
#include "Enclave_t.h"
#include "tSgxSSL_api.h"
#include <vector>
#include <cstdio>
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <algorithm>
#include <cmath>
#include <string>

using namespace std;

/**
 * OCall打印
 * @param fmt
 * @param ...
 * @return
 */
int printf(const char *fmt, ...)
{
    char buf[BUFSIZ] = {'\0'};
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, BUFSIZ, fmt, ap);
    va_end(ap);
    ocall_print_string(buf);
    return (int) strnlen(buf, BUFSIZ - 1) + 1;
}

/**
 * AES加密
 * @param plaintext
 * @param plaintext_len
 * @param aad
 * @param aad_len
 * @param key
 * @param iv
 * @param iv_len
 * @param ciphertext
 * @param tag
 * @return
 */
int gcm_encrypt(unsigned char *plaintext, int plaintext_len,
                unsigned char *aad, int aad_len,
                unsigned char *key,
                unsigned char *iv, int iv_len,
                unsigned char *ciphertext,
                unsigned char *tag)
{
    EVP_CIPHER_CTX *ctx;

    int len;
    int ciphertext_len;

    /* Create and initialise the context */
    if (!(ctx = EVP_CIPHER_CTX_new()))
        printf("errors");

    /* Initialise the encryption operation. */
    if (1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL))
        printf("errors");

    /*
     * Set IV length if default 12 bytes (96 bits) is not appropriate
     */
    if (1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, iv_len, NULL))
        printf("errors");

    /* Initialise key and IV */
    if (1 != EVP_EncryptInit_ex(ctx, NULL, NULL, key, iv))
        printf("errors");

    /*
     * Provide any AAD data. This can be called zero or more times as
     * required
     */
    if (1 != EVP_EncryptUpdate(ctx, NULL, &len, aad, aad_len))
        printf("errors");

    /*
     * Provide the message to be encrypted, and obtain the encrypted output.
     * EVP_EncryptUpdate can be called multiple times if necessary
     */
    if (1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintext_len))
        printf("errors");
    ciphertext_len = len;

    /*
     * Finalise the encryption. Normally ciphertext bytes may be written at
     * this stage, but this does not occur in GCM mode
     */
    if (1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len))
        printf("errors");
    ciphertext_len += len;

    /* Get the tag */
    if (1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, tag))
        printf("errors");

    /* Clean up */
    EVP_CIPHER_CTX_free(ctx);

    return ciphertext_len;
}

/**
 * AES解密
 * @param ciphertext
 * @param ciphertext_len
 * @param aad
 * @param aad_len
 * @param tag
 * @param key
 * @param iv
 * @param iv_len
 * @param plaintext
 * @return
 */
int gcm_decrypt(unsigned char *ciphertext, int ciphertext_len,
                unsigned char *aad, int aad_len,
                unsigned char *tag,
                unsigned char *key,
                unsigned char *iv, int iv_len,
                unsigned char *plaintext)
{
    EVP_CIPHER_CTX *ctx;
    int len;
    int plaintext_len;
    int ret;

    /* Create and initialise the context */
    if (!(ctx = EVP_CIPHER_CTX_new()))
        printf("errors");

    /* Initialise the decryption operation. */
    if (!EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL))
        printf("errors");

    /* Set IV length. Not necessary if this is 12 bytes (96 bits) */
    if (!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, iv_len, NULL))
        printf("errors");

    /* Initialise key and IV */
    if (!EVP_DecryptInit_ex(ctx, NULL, NULL, key, iv))
        printf("errors");

    /*
     * Provide any AAD data. This can be called zero or more times as
     * required
     */
    if (!EVP_DecryptUpdate(ctx, NULL, &len, aad, aad_len))
        printf("errors");

    /*
     * Provide the message to be decrypted, and obtain the plaintext output.
     * EVP_DecryptUpdate can be called multiple times if necessary
     */
    if (!EVP_DecryptUpdate(ctx, plaintext, &len, ciphertext, ciphertext_len))
        printf("errors");
    plaintext_len = len;

    /* Set expected tag value. Works in OpenSSL 1.0.1d and later */
    if (!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, 16, tag))
        printf("errors");

    /*
     * Finalise the decryption. A positive return value indicates success,
     * anything else is a failure - the plaintext is not trustworthy.
     */
    ret = EVP_DecryptFinal_ex(ctx, plaintext + len, &len);

    /* Clean up */
    EVP_CIPHER_CTX_free(ctx);

    if (ret > 0)
    {
        /* Success */
        plaintext_len += len;
        return plaintext_len;
    } else
    {
        /* Verify failed */
        return -1;
    }
}


/**
 * 将 vector<vector<int>> 转换为字节数组
 * @param vec
 * @param out_len
 * @return
 */
unsigned char *vectorToBytes_two_dimensional(const vector<vector<int>> &vec, int &out_len)
{
    // 计算总字节大小：所有元素的数据 + 每个子向量的大小信息（子向量长度* sizeof(int)）
    int total_size = sizeof(int) * vec.size();  // 存储每个子向量长度的大小
    for (const auto &sub_vec: vec)
    {
        total_size += sub_vec.size() * sizeof(int);
    }

    // 分配字节数组
    unsigned char *bytes = new unsigned char[total_size];
    unsigned char *ptr = bytes;  // 遍历指针

    // 逐个序列化每个子向量
    for (const auto &sub_vec: vec)
    {
        int sub_size = sub_vec.size();

        // 将当前子向量的大小写入字节数组
        memcpy(ptr, &sub_size, sizeof(int));
        ptr += sizeof(int);

        // 将子向量的数据写入字节数组
        memcpy(ptr, sub_vec.data(), sub_size * sizeof(int));
        ptr += sub_size * sizeof(int);
    }

    out_len = total_size;
    return bytes;
}

/**
 * 将字节数组转换回 vector<vector<int>>
 * @param bytes
 * @param byte_len
 * @return
 */
vector<vector<int>> bytesToVector_two_dimensional(const unsigned char *bytes, int byte_len)
{
    vector<vector<int>> vec;
    const unsigned char *ptr = bytes;  // 遍历指针

    while (ptr < bytes + byte_len)
    {
        // 读取子向量的大小
        int sub_size;
        memcpy(&sub_size, ptr, sizeof(int));
        ptr += sizeof(int);

        // 读取子向量的数据
        vector<int> sub_vec(sub_size);
        memcpy(sub_vec.data(), ptr, sub_size * sizeof(int));
        ptr += sub_size * sizeof(int);

        // 将子向量添加到主向量中
        vec.push_back(sub_vec);
    }

    return vec;
}

/**
 * 将 vector<double> 转换为字节数组
 * @param vec
 * @param out_len
 * @return
 */
unsigned char *DoubleVectorToBytes(const vector<double> &vec, int &out_len)
{
    out_len = vec.size() * sizeof(double);
    unsigned char *bytes = new unsigned char[out_len];
    memcpy(bytes, vec.data(), out_len);
    return bytes;
}

/**
 * 将字节数组转换回 vector<double>
 * @param bytes
 * @param byte_len
 * @return
 */
vector<double> bytesToDoubleVector(const unsigned char *bytes, int byte_len)
{
    vector<double> vec(byte_len / sizeof(double));
    memcpy(vec.data(), bytes, byte_len);
    return vec;
}


/**
 * 将 vector<int> 转换为字节数组
 * @param vec
 * @param out_len
 * @return
 */
unsigned char *vectorToBytes_one_dimensional(const vector<int> &vec, int &out_len)
{
    out_len = vec.size() * sizeof(int);
    unsigned char *bytes = new unsigned char[out_len];
    memcpy(bytes, vec.data(), out_len);
    return bytes;
}

/**
 * 将字节数组转换回 vector<int>
 * @param bytes
 * @param byte_len
 * @return
 */
vector<int> bytesToVector_one_dimensional(const unsigned char *bytes, int byte_len)
{
    vector<int> vec(byte_len / sizeof(int));
    memcpy(vec.data(), bytes, byte_len);
    return vec;
}


/**
 * 将 double 转换为字节数组
 * @param value
 * @return
 */
unsigned char *doubleToBytes(double value)
{
    unsigned char *bytes = new unsigned char[sizeof(double)];
    memcpy(bytes, &value, sizeof(double));
    return bytes;
}

/**
 * 将字节数组转换回 double
 * @param bytes
 * @return
 */
double bytesToDouble(const unsigned char *bytes)
{
    double value;
    memcpy(&value, bytes, sizeof(double));
    return value;
}

/**
 * 二维整型加密
 * @param vtr
 * @param key
 * @param iv
 * @param iv_len
 * @param tag
 * @param ciphertext
 * @param ciphertext_len
 */
void encryptVector(
        const vector<vector<int>> &vtr,                // 二维数据
        unsigned char *key,            // AES 密钥
        unsigned char *iv, int iv_len, // IV 和 IV 长度
        unsigned char **tag,                  // 输出 Tag
        unsigned char **ciphertext,          // 输出密文（通过指针）
        int &ciphertext_len                  // 输出密文长度（通过引用）
)
{
    // 1. 将 vector<int> 转换为字节数组
    int plaintext_len;
    unsigned char *plaintext = vectorToBytes_two_dimensional(vtr, plaintext_len);
    printf("[Enclave] int类型 二维数组明文长度 = %d\n", plaintext_len);

    // 2. 分配密文空间
    *ciphertext = new unsigned char[plaintext_len];

    // 3. 调用 AES-GCM 加密函数
    ciphertext_len = gcm_encrypt(
            plaintext, plaintext_len,
            nullptr, 0,   // 无 AAD
            key, iv, iv_len,
            *ciphertext, *tag
    );

    // 4. 释放明文内存
    delete[] plaintext;
}

/**
 * 一维整型加密
 * @param vtr
 * @param key
 * @param iv
 * @param iv_len
 * @param tag
 * @param ciphertext
 * @param ciphertext_len
 */
void encryptVector(
        const vector<int> &vtr,                // 一维数据
        unsigned char *key,            // AES 密钥
        unsigned char *iv, int iv_len, // IV 和 IV 长度
        unsigned char **tag,                  // 输出 Tag
        unsigned char **ciphertext,          // 输出密文（通过指针）
        int &ciphertext_len                  // 输出密文长度（通过引用）
)
{
    // 1. 将 vector<int> 转换为字节数组
    int plaintext_len;
    unsigned char *plaintext = vectorToBytes_one_dimensional(vtr, plaintext_len);
    printf("[Enclave] int类型 一维数组明文长度 = %d\n", plaintext_len);

    // 2. 分配密文空间
    *ciphertext = new unsigned char[plaintext_len];

    // 3. 调用 AES-GCM 加密函数
    ciphertext_len = gcm_encrypt(
            plaintext, plaintext_len,
            nullptr, 0,   // 无 AAD
            key, iv, iv_len,
            *ciphertext, *tag
    );

    // 4. 释放明文内存
    delete[] plaintext;
}


/**
 * 一维浮点型加密
 * @param vtr
 * @param key
 * @param iv
 * @param iv_len
 * @param tag
 * @param ciphertext
 * @param ciphertext_len
 */
void encryptVector(
        const vector<double> &vtr,                // 一维double数据
        unsigned char *key,            // AES 密钥
        unsigned char *iv, int iv_len, // IV 和 IV 长度
        unsigned char **tag,                  // 输出 Tag
        unsigned char **ciphertext,          // 输出密文（通过指针）
        int &ciphertext_len                  // 输出密文长度（通过引用）
)
{
    // 1. 将 vector<double> 转换为字节数组
    int plaintext_len;
    unsigned char *plaintext = DoubleVectorToBytes(vtr, plaintext_len);
    printf("[Enclave] double类型 一维数组明文长度 = %d\n", plaintext_len);

    // 2. 分配密文空间
    *ciphertext = new unsigned char[plaintext_len];

    // 3. 调用 AES-GCM 加密函数
    ciphertext_len = gcm_encrypt(
            plaintext, plaintext_len,
            nullptr, 0,   // 无 AAD
            key, iv, iv_len,
            *ciphertext, *tag
    );

    // 4. 释放明文内存
    delete[] plaintext;
}

/**
 * 二维整型解密
 * @param ciphertext
 * @param ciphertext_len
 * @param key
 * @param iv
 * @param iv_len
 * @param tag
 * @param vtr
 */
void decryptVector(
        unsigned char *ciphertext, int ciphertext_len, // 输入密文及其长度
        unsigned char *key, unsigned char *iv, int iv_len, // AES 密钥和 IV
        unsigned char *tag, // GCM Tag
        vector<vector<int>> &vtr // 输出解密后的二维数据
)
{
    // 1. 分配解密输出缓冲区
    unsigned char *decryptedtext = new unsigned char[120000];  // 根据数据大小动态调整

    // 2. 调用 AES-GCM 解密函数
    int decrypted_len = gcm_decrypt(
            ciphertext, ciphertext_len,
            nullptr, 0,   // 无 AAD
            tag,
            key, iv, iv_len,
            decryptedtext
    );
    // 3. 检查解密结果
    if (decrypted_len >= 0)
    {
        vtr = bytesToVector_two_dimensional(decryptedtext, decrypted_len);
    } else
    {
        printf("[Enclave] ERROR 解密失败！\n");
        throw;
    }

    // 4. 释放解密输出缓冲区内存
    delete[] decryptedtext;
}

/**
 * 一维整型解密
 * @param ciphertext
 * @param ciphertext_len
 * @param key
 * @param iv
 * @param iv_len
 * @param tag
 * @param vtr
 */
void decryptVector(
        unsigned char *ciphertext, int ciphertext_len, // 输入密文及其长度
        unsigned char *key, unsigned char *iv, int iv_len, // AES 密钥和 IV
        unsigned char *tag, // GCM Tag
        vector<int> &vtr // 输出解密后的一维数据
)
{
    // 1. 分配解密输出缓冲区
    unsigned char *decryptedtext = new unsigned char[20000];  // 根据数据大小动态调整

    // 2. 调用 AES-GCM 解密函数
    int decrypted_len = gcm_decrypt(
            ciphertext, ciphertext_len,
            nullptr, 0,   // 无 AAD
            tag,
            key, iv, iv_len,
            decryptedtext
    );

    // 3. 检查解密结果
    if (decrypted_len >= 0)
    {
        vtr = bytesToVector_one_dimensional(decryptedtext, decrypted_len);
    } else
    {
        printf("[Enclave] ERROR 解密失败！\n");
        throw;
    }

    // 4. 释放解密输出缓冲区内存
    delete[] decryptedtext;
}

/**
 * 一维浮点型解密
 * @param ciphertext
 * @param ciphertext_len
 * @param key
 * @param iv
 * @param iv_len
 * @param tag
 * @param vtr
 */
void decryptVector(
        unsigned char *ciphertext, int ciphertext_len, // 输入密文及其长度
        unsigned char *key, unsigned char *iv, int iv_len, // AES 密钥和 IV
        unsigned char *tag, // GCM Tag
        vector<double> &vtr // 输出解密后的一维数据
)
{
    // 1. 分配解密输出缓冲区
    unsigned char *decryptedtext = new unsigned char[40];  // 根据数据大小动态调整

    // 2. 调用 AES-GCM 解密函数
    int decrypted_len = gcm_decrypt(
            ciphertext, ciphertext_len,
            nullptr, 0,   // 无 AAD
            tag,
            key, iv, iv_len,
            decryptedtext
    );

    // 3. 检查解密结果
    if (decrypted_len >= 0)
    {
        vtr = bytesToDoubleVector(decryptedtext, decrypted_len);
    } else
    {
        printf("[Enclave] ERROR 解密失败！\n");
        throw;
    }

    // 4. 释放解密输出缓冲区内存
    delete[] decryptedtext;
}

/**
 * 将十六进制字符串转换为字节数组
 * @param hex
 * @param buffer
 * @param bufferSize
 * @return
 */
bool hexStringToBytes(const string &hex, unsigned char *buffer, size_t bufferSize)
{
    if (hex.length() != bufferSize * 2)
    {
        printf("[Enclave] Hex string length doesn't match buffer size.\n");
        return false;
    }
    for (size_t i = 0; i < bufferSize; ++i)
    {
        string byteString = hex.substr(i * 2, 2);  // 取两个字符
        buffer[i] = static_cast<unsigned char>( stoul(byteString, nullptr, 16));
    }
    return true;
}


/**
 * 多元线性回归的预测函数：y_hat = w1 * X[i][0] + w2 * X[i][1] + ... + wn * X[i][n-1] + b
 * @param x
 * @param w
 * @param b
 * @return
 */
double predict(const vector<int> &x, const vector<double> &w, double b)
{
    double y_hat = b;
    for (size_t j = 0; j < x.size(); ++j)
    {
        y_hat += w[j] * x[j];
    }
    return y_hat;
}

/**
 * 计算误差损失函数 (Mean Squared Error)
 * @param X
 * @param Y
 * @param w
 * @param b
 * @return
 */
double computeLoss(const vector<vector<int>> &X, const vector<int> &Y, const vector<double> &w, double b)
{
    double loss = 0.0;
    int n = X.size();
    for (int i = 0; i < n; ++i)
    {
        double y_hat = predict(X[i], w, b);
        loss += pow(y_hat - Y[i], 2);
    }
    return loss / n;
}

/**
 * 打印十六进制
 * @param data
 * @param len
 */
void printHex(const unsigned char *data, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        printf("%02x", data[i]);  // %02x：两位十六进制，不足补0
    }
    printf("\n");
}

/**
 * Ecall 梯度下降
 * @param ciphertext_X_train
 * @param ciphertext_len_X_train
 * @param tag_X_train
 * @param ciphertext_Y_train
 * @param ciphertext_len_Y_train
 * @param tag_Y_train
 * @param ciphertext_w
 * @param ciphertext_len_w
 * @param tag_w
 * @param ciphertext_b
 * @param ciphertext_len_b
 * @param tag_b
 * @param ciphertext_learningRate
 * @param ciphertext_len_learningRate
 * @param tag_learningRate
 * @param epochs
 * @param batchSize
 * @param ciphertext_w_final
 * @param ciphertext_len_w_final
 * @param tag_w_final
 * @param ciphertext_b_final
 * @param ciphertext_len_b_final
 * @param tag_b_final
 */
void Enclave_gradientDescent(unsigned char **ciphertext_X_train, size_t ciphertext_len_X_train, unsigned char **tag_X_train,
                     unsigned char **ciphertext_Y_train, size_t ciphertext_len_Y_train, unsigned char **tag_Y_train,
                     unsigned char **ciphertext_w, size_t ciphertext_len_w, unsigned char **tag_w,
                     unsigned char **ciphertext_b, size_t ciphertext_len_b, unsigned char **tag_b,
                     unsigned char **ciphertext_learningRate, size_t ciphertext_len_learningRate, unsigned char **tag_learningRate,
                     size_t epochs, size_t batchSize,
                     unsigned char **ciphertext_w_final, size_t *ciphertext_len_w_final, unsigned char **tag_w_final,
                     unsigned char **ciphertext_b_final, size_t *ciphertext_len_b_final, unsigned char **tag_b_final)
{
    // 密钥和 IV 的十六进制字符串
    string keyHex = "af46723b2704a46fc22415e4093f47d8ed19fd52ca857c2929b0669c9e5c026c";
    string ivHex = "869b58e4d0d3ed27c419cbad0cd760d9";

    unsigned char key[32];  // 256 位密钥
    unsigned char iv[16];   // 128 位 IV

    // 转换并赋值
    hexStringToBytes(keyHex, key, sizeof(key));
    hexStringToBytes(ivHex, iv, sizeof(iv));
    // 输出 key 和 iv 的值
    printf("[Enclave] Key: ");
    printHex(key, sizeof(key));

    printf("[Enclave] IV: ");
    printHex(iv, sizeof(iv));

    // 解密得到明文的数据集X与Y
    vector<vector<int>> X;
    decryptVector(*ciphertext_X_train, ciphertext_len_X_train, key, iv, sizeof(iv), *tag_X_train, X);
    printf("[Enclave]  SUCCESS X_train 解密成功！\n");

    vector<int> Y;
    decryptVector(*ciphertext_Y_train, ciphertext_len_Y_train, key, iv, sizeof(iv), *tag_Y_train, Y);
    printf("[Enclave]  SUCCESS Y_train 解密成功！\n");

    // 解密得到明文的模型权重w
    vector<double> w;
    decryptVector(*ciphertext_w, ciphertext_len_w, key, iv, sizeof(iv), *tag_w, w);
    printf("[Enclave]  SUCCESS w 解密成功！\n");

    // 解密偏置项b
    unsigned char *decryptedtext_b = new unsigned char[sizeof(double)];
    int decrypted_len_b;
    decrypted_len_b = gcm_decrypt(
            *ciphertext_b, ciphertext_len_b,
            nullptr, 0,   // 无 AAD
            *tag_b,
            key, iv, sizeof(iv),
            decryptedtext_b  // 存储解密后的数据
    );
    double b = bytesToDouble(decryptedtext_b);
    printf("[Enclave]  SUCCESS b 解密成功！\n");

    // 解密学习率learningRate
    unsigned char *decryptedtext_learningRate = new unsigned char[sizeof(double)];
    int decrypted_len_learningRate;
    decrypted_len_learningRate = gcm_decrypt(
            *ciphertext_learningRate, ciphertext_len_learningRate,
            nullptr, 0,   // 无 AAD
            *tag_learningRate,
            key, iv, sizeof(iv),
            decryptedtext_learningRate  // 存储解密后的数据
    );
    double learningRate = bytesToDouble(decryptedtext_learningRate);
    printf("[Enclave]  SUCCESS learningRate 解密成功！\n");

    int n = X.size();       // 样本数
    int features = X[0].size(); // 特征数

    printf("=== （明文训练）线性回归模型配置 ===\n");
    printf("[Enclave] 样本数量: %d\n", n);
    printf("[Enclave] 特征数量: %d\n", features);
    printf("[Enclave] 学习率: %f\n", learningRate);
    printf("[Enclave] 训练轮数: %d\n", epochs);
    printf("[Enclave] 批次大小: %d\n", batchSize);

    printf("[Enclave] 初始权重: ");
    for (const auto &weight: w)
    {
        printf("%f ", weight);
    }
    printf("\n");

    printf("[Enclave] 初始偏置项: %f\n", b);
    printf("=========================\n");

    for (int epoch = 1; epoch <= epochs; epoch++)
    {
        // 遍历所有数据，按批次大小进行训练
        for (int i = 0; i < n; i += batchSize)
        {
            vector<double> dw(features, 0.0);
            double db = 0.0;

            // 处理当前批次的数据
            size_t currentBatchSize = min(batchSize, static_cast<size_t>(n - i)); // 防止越界
            for (int j = i; j < i + currentBatchSize; ++j)
            {
                double y_hat = predict(X[j], w, b);
                double error = y_hat - Y[j];
                for (int k = 0; k < features; ++k)
                {
                    dw[k] += error * X[j][k];
                }
                db += error;
            }



            // 平均梯度
            for (int k = 0; k < features; ++k)
            {
                dw[k] *= 2;
                dw[k] /= currentBatchSize;
            }
            db *= 2;
            db /= currentBatchSize;

            // 更新参数
            for (int k = 0; k < features; ++k)
            {
                w[k] -= learningRate * dw[k];
            }
            b -= learningRate * db;


        }

        double loss = computeLoss(X, Y, w, b);
        printf("[Enclave] Epoch %d, Loss: %f, b: %f\n", epoch, loss, b);
        printf("[Enclave] Weights: ");
        for (const auto &weight: w)
        {
            printf("%f ", weight);
        }
        printf("\n");
    }

    printf("\n[Enclave] 训练结束。最终权重和偏置：\n");
    printf("[Enclave] 权重: ");
    for (const auto &weight: w)
    {
        printf("%f ", weight);
    }
    printf("\n");
    printf("[Enclave] 偏置项: %f\n", b);


    // 加密最终模型权重w

    int plaintext_len;
    unsigned char *plaintext = DoubleVectorToBytes(w, plaintext_len);

    // 3. 调用 AES-GCM 加密函数
    *ciphertext_len_w_final = gcm_encrypt(
            plaintext, plaintext_len,
            nullptr, 0,   // 无 AAD
            key, iv, sizeof(iv),
            *ciphertext_w_final, *tag_w_final
    );
    delete[] plaintext;

    // 加密最终模型偏置项b
    unsigned char *plaintext_b = doubleToBytes(b);
    *ciphertext_len_b_final = gcm_encrypt(
            plaintext_b, sizeof(double),
            nullptr, 0,   // 无 AAD
            key, iv, sizeof(iv),
            *ciphertext_b_final, *tag_b_final
    );

//    printHex(*ciphertext_w_final, sizeof(*ciphertext_w_final) / sizeof(*ciphertext_w_final[0]));
//    printHex(*ciphertext_b_final, sizeof(*ciphertext_b_final) / sizeof(*ciphertext_b_final[0]));
//    printHex(*tag_w_final, sizeof(*tag_w_final) / sizeof(*tag_w_final[0]));
//    printHex(*tag_b_final, sizeof(*tag_b_final) / sizeof(*tag_b_final[0]));
//    printf("%d\n", *ciphertext_len_w_final);
//    printf("%d\n", *ciphertext_len_b_final);
}





