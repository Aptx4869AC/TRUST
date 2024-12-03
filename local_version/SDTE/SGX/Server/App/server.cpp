/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-02
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include <iostream>
#include <stdio.h>
#include <string.h>
#include <assert.h>

#include <unistd.h>
#include <pwd.h>
#include <libgen.h>
#include <stdlib.h>
#include <pthread.h>
#include <cstdio>
#include <chrono>
#include <thread>
#include <omp.h>
#include <iomanip>
#include <vector>
#include <algorithm>
#include <cmath>
#include <numeric>
#include <fstream>
#include <sstream>
#include <openssl/conf.h>
#include <openssl/evp.h>
#include <openssl/err.h>
#include <openssl/rsa.h>
#include <openssl/pem.h>
#include <openssl/rand.h>

#include <sgx_urts.h>
#include "server.h"
#include "Enclave_u.h"

#define MAX_PATH FILENAME_MAX
using namespace std;

sgx_enclave_id_t global_eid = 0;

typedef struct _sgx_errlist_t
{
    sgx_status_t err;
    const char *msg;
    const char *sug; /* Suggestion */
} sgx_errlist_t;

/**
 * Error code returned by sgx_create_enclave
 */
static sgx_errlist_t sgx_errlist[] =
        {
                {
                        SGX_ERROR_UNEXPECTED, "Unexpected error occurred.", NULL
                }, {
                SGX_ERROR_INVALID_PARAMETER, "Invalid parameter.", NULL
        }, {
                SGX_ERROR_OUT_OF_MEMORY, "Out of memory.", NULL
        }, {
                SGX_ERROR_ENCLAVE_LOST, "Power transition occurred.", "Please refer to the sample \"PowerTransition\" for details."
        }, {
                SGX_ERROR_INVALID_ENCLAVE, "Invalid enclave image.", NULL
        }, {
                SGX_ERROR_INVALID_ENCLAVE_ID, "Invalid enclave identification.", NULL
        }, {
                SGX_ERROR_INVALID_SIGNATURE, "Invalid enclave signature.", NULL
        }, {
                SGX_ERROR_OUT_OF_EPC, "Out of EPC memory.", NULL
        }, {
                SGX_ERROR_NO_DEVICE, "Invalid Intel® Software Guard Extensions device.", "Please make sure Intel® Software Guard Extensions module is enabled in the BIOS, and install Intel® Software Guard Extensions driver afterwards."
        }, {
                SGX_ERROR_MEMORY_MAP_CONFLICT, "Memory map conflicted.", NULL
        }, {
                SGX_ERROR_INVALID_METADATA, "Invalid enclave metadata.", NULL
        }, {
                SGX_ERROR_DEVICE_BUSY, "Intel® Software Guard Extensions device was busy.", NULL
        }, {
                SGX_ERROR_INVALID_VERSION, "Enclave version was invalid.", NULL
        }, {
                SGX_ERROR_INVALID_ATTRIBUTE, "Enclave was not authorized.", NULL
        }, {
                SGX_ERROR_ENCLAVE_FILE_ACCESS, "Can't open enclave file.", NULL
        },};

/**
 * Check error conditions for loading enclave
 * @param ret
 */
void print_error_message(sgx_status_t ret)
{
    size_t idx = 0;
    size_t ttl = sizeof sgx_errlist / sizeof sgx_errlist[0];

    for (idx = 0; idx < ttl; idx++)
    {
        if (ret == sgx_errlist[idx].err)
        {
            if (NULL != sgx_errlist[idx].sug)
                printf("Info: %s\n", sgx_errlist[idx].sug);
            printf("Error: %s\n", sgx_errlist[idx].msg);
            break;
        }
    }

    if (idx == ttl)
        printf("Error: Unexpected error occurred [0x%x].\n", ret);
}


/**
 * Initialize the enclave:
 * Step 1: retrive the launch token saved by last transaction
 * Step 2: call sgx_create_enclave to initialize an enclave instance
 * Step 3: save the launch token if it is updated
 * @return
 */
int initialize_enclave(void)
{
    char token_path[MAX_PATH] = {'\0'};
    sgx_launch_token_t token = {0};
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;
    int updated = 0;
    /* Step 1: retrive the launch token saved by last transaction */

    /* try to get the token saved in $HOME */
    char cwd[1024];
    const char *home_dir = getcwd(cwd, sizeof(cwd));
    if (home_dir != NULL &&
        (strlen(home_dir) + strlen("/") + sizeof(TOKEN_FILENAME) + 1) <= MAX_PATH)
    {
        /* compose the token path */
        strncpy(token_path, home_dir, strlen(home_dir));
        strncat(token_path, "/", strlen("/"));
        strncat(token_path, TOKEN_FILENAME, sizeof(TOKEN_FILENAME) + 1);
    } else
    {
        /* if token path is too long or $HOME is NULL */
        strncpy(token_path, TOKEN_FILENAME, sizeof(token_path));
    }

    FILE *fp = fopen(token_path, "rb");
    if (fp == NULL && (fp = fopen(token_path, "wb")) == NULL)
    {
        printf("Warning: Failed to create/open the launch token file \"%s\".\n", token_path);
    }
    printf("token_path: %s\n", token_path);
    if (fp != NULL)
    {
        /* read the token from saved file */
        size_t read_num = fread(token, 1, sizeof(sgx_launch_token_t), fp);
        if (read_num != 0 && read_num != sizeof(sgx_launch_token_t))
        {
            /* if token is invalid, clear the buffer */
            memset(&token, 0x0, sizeof(sgx_launch_token_t));
            printf("Warning: Invalid launch token read from \"%s\".\n", token_path);
        }
    }

    /* Step 2: call sgx_create_enclave to initialize an enclave instance */
    /* Debug Support: set 2nd parameter to 1 */

    ret = sgx_create_enclave(ENCLAVE_FILENAME, SGX_DEBUG_FLAG, &token, &updated, &global_eid, NULL);

    if (ret != SGX_SUCCESS)
    {
        print_error_message(ret);
        if (fp != NULL) fclose(fp);

        return -1;
    }

    /* Step 3: save the launch token if it is updated */

    if (updated == FALSE || fp == NULL)
    {
        /* if the token is not updated, or file handler is invalid, do not perform saving */
        if (fp != NULL) fclose(fp);
        return 0;
    }

    /* reopen the file with write capablity */
    fp = freopen(token_path, "wb", fp);
    if (fp == NULL) return 0;
    size_t write_num = fwrite(token, 1, sizeof(sgx_launch_token_t), fp);
    if (write_num != sizeof(sgx_launch_token_t))
        printf("Warning: Failed to save launch token to \"%s\".\n", token_path);
    fclose(fp);

    return 0;
}

/**
 * OCall functions
 * @param str
 */
void ocall_print_string(const char *str)
{
    printf("%s", str);
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
 * 计算 R^2 (决定系数)
 * @param X
 * @param Y
 * @param w
 * @param b
 * @return
 */
double computeR2(const vector<vector<int>> &X, const vector<int> &Y, const vector<double> &w, double b)
{
    double ss_tot = 0.0, ss_res = 0.0;
    double mean_y = accumulate(Y.begin(), Y.end(), 0.0) / Y.size();

    for (size_t i = 0; i < X.size(); ++i)
    {
        double y_hat = predict(X[i], w, b);
        ss_res += pow(Y[i] - y_hat, 2);
        ss_tot += pow(Y[i] - mean_y, 2);
    }

    return 1 - (ss_res / ss_tot);
}

/**
 * 计算 MAE (Mean Absolute Error)
 * @param X
 * @param Y
 * @param w
 * @param b
 * @return
 */
double computeMAE(const vector<vector<int>> &X, const vector<int> &Y, const vector<double> &w, double b)
{
    double mae = 0.0;
    for (size_t i = 0; i < X.size(); ++i)
    {
        double y_hat = predict(X[i], w, b);
        mae += abs(y_hat - Y[i]);
    }
    return mae / X.size();
}

/**
 * 评价函数
 * @param X
 * @param Y
 * @param w
 * @param b
 */
void evaluate(vector<vector<int>> &X, vector<int> &Y, vector<double> &w, double &b)
{
    // 计算评价指标
    double mse = computeLoss(X, Y, w, b);
    double r2 = computeR2(X, Y, w, b);
    double mae = computeMAE(X, Y, w, b);

    cout << "------------------------------------------------------\n";
    cout << "[APP] MSE: " << mse << '\n';
    cout << "[APP] R²: " << r2 << '\n';
    cout << "[APP] MAE: " << mae << '\n';
    cout << "------------------------------------------------------\n";

}


void handleErrors(void)
{
    ERR_print_errors_fp(stderr);
    abort();
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
        handleErrors();

    /* Initialise the encryption operation. */
    if (1 != EVP_EncryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL))
        handleErrors();

    /*
     * Set IV length if default 12 bytes (96 bits) is not appropriate
     */
    if (1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, iv_len, NULL))
        handleErrors();

    /* Initialise key and IV */
    if (1 != EVP_EncryptInit_ex(ctx, NULL, NULL, key, iv))
        handleErrors();

    /*
     * Provide any AAD data. This can be called zero or more times as
     * required
     */
    if (1 != EVP_EncryptUpdate(ctx, NULL, &len, aad, aad_len))
        handleErrors();

    /*
     * Provide the message to be encrypted, and obtain the encrypted output.
     * EVP_EncryptUpdate can be called multiple times if necessary
     */
    if (1 != EVP_EncryptUpdate(ctx, ciphertext, &len, plaintext, plaintext_len))
        handleErrors();
    ciphertext_len = len;

    /*
     * Finalise the encryption. Normally ciphertext bytes may be written at
     * this stage, but this does not occur in GCM mode
     */
    if (1 != EVP_EncryptFinal_ex(ctx, ciphertext + len, &len))
        handleErrors();
    ciphertext_len += len;

    /* Get the tag */
    if (1 != EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_GET_TAG, 16, tag))
        handleErrors();

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
        handleErrors();

    /* Initialise the decryption operation. */
    if (!EVP_DecryptInit_ex(ctx, EVP_aes_256_gcm(), NULL, NULL, NULL))
        handleErrors();

    /* Set IV length. Not necessary if this is 12 bytes (96 bits) */
    if (!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_IVLEN, iv_len, NULL))
        handleErrors();

    /* Initialise key and IV */
    if (!EVP_DecryptInit_ex(ctx, NULL, NULL, key, iv))
        handleErrors();

    /*
     * Provide any AAD data. This can be called zero or more times as
     * required
     */
    if (!EVP_DecryptUpdate(ctx, NULL, &len, aad, aad_len))
        handleErrors();

    /*
     * Provide the message to be decrypted, and obtain the plaintext output.
     * EVP_DecryptUpdate can be called multiple times if necessary
     */
    if (!EVP_DecryptUpdate(ctx, plaintext, &len, ciphertext, ciphertext_len))
        handleErrors();
    plaintext_len = len;

    /* Set expected tag value. Works in OpenSSL 1.0.1d and later */
    if (!EVP_CIPHER_CTX_ctrl(ctx, EVP_CTRL_GCM_SET_TAG, 16, tag))
        handleErrors();

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
 * 打印16进制
 * @param data
 * @param len
 */
void printHex(const unsigned char *data, size_t len)
{
    for (size_t i = 0; i < len; ++i)
    {
        cout << hex << setw(2) << setfill('0') << static_cast<int>(data[i]);
    }
    cout << '\n';
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
        cout << hex.length() << '\n';
        cout << bufferSize << '\n';
        cerr << "Hex string length doesn't match buffer size." << endl;
        throw;
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
 * 读取CSV文件并将内容存储到 X 和 Y 中
 * @param filePath
 * @param X
 * @param Y
 * @return
 */
bool loadCSV(const string &filePath, vector<vector<int>> &X, vector<int> &Y)
{
    ifstream file(filePath);
    if (!file)
    {
        cerr << "文件无法打开: " << strerror(errno) << '\n';
        return false;  // 打开文件失败
    }

    string line, value;
    getline(file, line);  // 跳过表头

    // 逐行读取CSV文件的内容
    while (getline(file, line))
    {
        stringstream ss(line);
        vector<int> row;  // 存储当前行的输入变量

        // 逐列读取，除了最后一列都存入 row
        while (getline(ss, value, ','))
        {
            row.push_back(stoi(value));
        }

        // 最后一列是预测值，存入 Y
        Y.push_back(row.back());
        row.pop_back();  // 删除最后一列

        // 剩下的列存入 X
        X.push_back(row);
    }

    file.close();  // 关闭文件
    return true;   // 文件读取成功
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


void printData(vector<vector<int>> &X, vector<int> &Y)
{
    cout << "读取的数据：" << '\n';
    for (size_t i = 0; i < X.size(); ++i)
    {
        cout << "样本 " << i + 1 << ": ";
        for (size_t j = 0; j < X[i].size(); ++j)
        {
            cout << X[i][j] << " ";
        }
        cout << "| 真实值: " << Y[i] << '\n';
    }
}

/**
 * 加密整型二维向量
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
    cout << "[APP] int类型 二维数组明文长度 = " << plaintext_len << '\n';

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
 * 加密整型一维向量
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
    cout << "[APP] int类型 一维数组明文长度 = " << plaintext_len << '\n';

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
 * 加密浮点型一维向量
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
    cout << "[APP] double类型 一维数组明文长度 = " << plaintext_len << '\n';

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
 * 解密整型二维向量
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
        cout << " ERROR 解密失败！" << '\n';
        throw;
    }

    // 4. 释放解密输出缓冲区内存
    delete[] decryptedtext;
}

/**
 * 解密整型一维向量
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
        cout << " ERROR 解密失败！" << '\n';
        throw;
    }

    // 4. 释放解密输出缓冲区内存
    delete[] decryptedtext;
}

/**
 * 解密浮点型一维向量
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
        cout << " ERROR 解密失败！" << '\n';
        throw;
    }

    // 4. 释放解密输出缓冲区内存
    delete[] decryptedtext;
}


/* Application entry */
int main(int argc, char **argv)
{
    (void) (argc);
    (void) (argv);

    /* Changing dir to where the executable is.*/
    char absolutePath[MAX_PATH];
    char *ptr = NULL;

    ptr = realpath(dirname(argv[0]), absolutePath);

    if (ptr == NULL || chdir(absolutePath) != 0)
        return 1;

    /* Initialize the enclave */
    if (initialize_enclave() < 0)
        return 1;

    double start_time, end_time;

    vector<vector<int>> X_train;  // 多维输入变量
    vector<int> Y_train;           // 预测值

    string train_filePath = "./Data/train_dataset.csv";
    if (loadCSV(train_filePath, X_train, Y_train))
    {
        cout << "[APP] 训练 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "[APP] 训练 CSV 文件读取失败！" << '\n';
    }

    vector<vector<int>> X_test;  // 多维输入变量
    vector<int> Y_test;           // 预测值

    string test_filePath = "./Data/test_dataset.csv";
    if (loadCSV(test_filePath, X_test, Y_test))
    {
        cout << "[APP] 测试 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "[APP] 测试 CSV 文件读取失败！" << '\n';
    }

    // 密钥和 IV 的十六进制字符串
    string keyHex = "af46723b2704a46fc22415e4093f47d8ed19fd52ca857c2929b0669c9e5c026c";
    string ivHex = "869b58e4d0d3ed27c419cbad0cd760d9";

    unsigned char key[32];  // 256 位密钥
    unsigned char iv[16];   // 128 位 IV

    // 转换并赋值
    hexStringToBytes(keyHex, key, sizeof(key));
    hexStringToBytes(ivHex, iv, sizeof(iv));

//    unsigned char *key = new unsigned char[32];
//    unsigned char *iv = new unsigned char[16];
//    // 初始化 key、iv
//    RAND_bytes(key, sizeof(key));
//    RAND_bytes(iv, sizeof(iv));

    // 输出 key 和 iv 的值
    cout << "[APP] Key: ";
    printHex(key, sizeof(key));

    cout << "[APP] IV: ";
    printHex(iv, sizeof(iv));


    // 声明密文指针和长度
    unsigned char *ciphertext_X_train = nullptr;
    int ciphertext_len_X_train = 0;
    unsigned char *tag_X_train = new unsigned char[16];  // 128 位 GCM Tag
    unsigned char *ciphertext_Y_train = nullptr;
    int ciphertext_len_Y_train = 0;
    unsigned char *tag_Y_train = new unsigned char[16];  // 128 位 GCM Tag

    // 数据集加密处理
    start_time = omp_get_wtime();
    encryptVector(X_train, key, iv, sizeof(iv), &tag_X_train, &ciphertext_X_train, ciphertext_len_X_train); // 加密 X_train
    encryptVector(Y_train, key, iv, sizeof(iv), &tag_Y_train, &ciphertext_Y_train, ciphertext_len_Y_train); // 加密 Y_train
    end_time = omp_get_wtime();
    printf("[APP] 加密 X_train 与 Y_train time is  ------  %f ms\n", (end_time - start_time) * 1000);


    int features = X_train[0].size(); // 获取特征数量
    // 初始化线性回归参数
    vector<double> w(features, 0.0); // 初始化权重为 0
    double b = 0.0;  // 偏置
    double learningRate = 0.0001;  // 学习率
    int epochs = 800;  // 训练轮数
    int batchSize = 64; // 批次大小

    // 加密w
    start_time = omp_get_wtime();
    unsigned char *ciphertext_w = nullptr;
    int ciphertext_len_w = 0;
    unsigned char *tag_w = new unsigned char[16];
    encryptVector(w, key, iv, sizeof(iv), &tag_w, &ciphertext_w, ciphertext_len_w);

    // 加密 b
    unsigned char *plaintext_b = doubleToBytes(b);
    unsigned char *ciphertext_b = new unsigned char[sizeof(double)];
    int ciphertext_len_b;
    unsigned char *tag_b = new unsigned char[16];
    ciphertext_len_b = gcm_encrypt(
            plaintext_b, sizeof(double),
            nullptr, 0,   // 无 AAD
            key, iv, sizeof(iv),
            ciphertext_b, tag_b
    );

    // 加密 learningRate
    unsigned char *plaintext_learningRate = doubleToBytes(learningRate);
    unsigned char *ciphertext_learningRate = new unsigned char[sizeof(double)];
    int ciphertext_len_learningRate;
    unsigned char *tag_learningRate = new unsigned char[16];
    ciphertext_len_learningRate = gcm_encrypt(
            plaintext_learningRate, sizeof(double),
            nullptr, 0,   // 无 AAD
            key, iv, sizeof(iv),
            ciphertext_learningRate, tag_learningRate
    );
    end_time = omp_get_wtime();
    printf("[APP] 加密 w、b、learningRate time is  ------  %f ms\n", (end_time - start_time) * 1000);

    // 记录返回的w与b
    unsigned char *ciphertext_w_final = new unsigned char[w.size() * sizeof(double)];
    size_t ciphertext_len_w_final = 0;
    unsigned char *tag_w_final = new unsigned char[16];
    unsigned char *ciphertext_b_final = new unsigned char[sizeof(double)];
    size_t ciphertext_len_b_final = 0;
    unsigned char *tag_b_final = new unsigned char[16];


    start_time = omp_get_wtime();
    // 使用梯度下降进行训练
    Enclave_gradientDescent(global_eid, &ciphertext_X_train, ciphertext_len_X_train, &tag_X_train, // 输入，密文
                            &ciphertext_Y_train, ciphertext_len_Y_train, &tag_Y_train, // 输入，密文
                            &ciphertext_w, ciphertext_len_w, &tag_w,// 输入，密文
                            &ciphertext_b, ciphertext_len_b, &tag_b, // 输入，密文
                            &ciphertext_learningRate, ciphertext_len_learningRate, &tag_learningRate, // 输入，密文
                            epochs, batchSize, // 输入，明文
                            &ciphertext_w_final, &ciphertext_len_w_final, &tag_w_final, // 输出，密文
                            &ciphertext_b_final, &ciphertext_len_b_final, &tag_b_final); // 输出，密文

    end_time = omp_get_wtime();
    cout << "------------------------------------------------------\n";
    printf("[Enclave] Ecall AES-LinR time is  ------  %f ms\n", (end_time - start_time) * 1000);
    cout << "------------------------------------------------------\n";

//    printHex(ciphertext_w_final, sizeof(ciphertext_w_final) / sizeof(ciphertext_w_final[0]));
//    printHex(ciphertext_b_final, sizeof(ciphertext_b_final) / sizeof(ciphertext_b_final[0]));
//    printHex(tag_w_final, sizeof(tag_w_final) / sizeof(tag_w_final[0]));
//    printHex(tag_b_final, sizeof(tag_b_final) / sizeof(tag_b_final[0]));
//    cout << ciphertext_len_w_final << '\n';
//    cout << ciphertext_len_b_final << '\n';


    // 解密 w
    vector<double> final_w;
    decryptVector(ciphertext_w_final, ciphertext_len_w_final, key, iv, sizeof(iv), tag_w_final, final_w);
    cout << "[APP] final_w 权重: ";
    for (double weight: final_w)
    {
        cout << weight << " ";
    }
    cout << '\n';


    // 分配解密后的缓冲区
    unsigned char *decryptedtext_b = new unsigned char[sizeof(double)];
    int decrypted_len_b;

    // 解密 b
    decrypted_len_b = gcm_decrypt(
            ciphertext_b_final, ciphertext_len_b_final,
            nullptr, 0,   // 无 AAD
            tag_b_final,
            key, iv, sizeof(iv),
            decryptedtext_b  // 存储解密后的数据
    );

    // 将字节数组转换回 double
    double final_b = bytesToDouble(decryptedtext_b);
    cout << "[APP] final_b = " << final_b << '\n';

    // 评估
    evaluate(X_test, Y_test, final_w, final_b);

    // destroy Encalve
    sgx_destroy_enclave(global_eid);
    printf("Info: successfully returned.\n");
    return 0;

}

