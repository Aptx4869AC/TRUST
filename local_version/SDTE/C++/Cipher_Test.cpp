/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-02
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include <omp.h>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <cstring>
#include <random>
#include <fstream>
#include <sstream>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include "./include/AES.h"
using namespace std;

unsigned char key[32];  // 256 位密钥
unsigned char iv[16];   // 128 位 IV


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
    cout << "MSE: " << mse << '\n';
    cout << "R²: " << r2 << '\n';
    cout << "MAE: " << mae << '\n';
    cout << "------------------------------------------------------\n";

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
 * 整型二维加密
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
    cout << "[int] 二维数组明文长度 = " << plaintext_len << '\n';

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
 * 浮点型一维加密
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
    cout << "[double] 一维数组明文长度 = " << plaintext_len << '\n';

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
 * 整型一维加密
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
    cout << "[int] 一维数组明文长度 = " << plaintext_len << '\n';

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
 * 整型二维解密
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
 * 浮点型一维解密
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


/**
 * 整型一维解密
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
 * 线性回归小批次梯度下降，更新 w 和 b 参数
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
void gradientDescent(unsigned char *ciphertext_X_train, int ciphertext_len_X_train, unsigned char *tag_X_train,
                     unsigned char *ciphertext_Y_train, int ciphertext_len_Y_train, unsigned char *tag_Y_train,
                     unsigned char *ciphertext_w, int ciphertext_len_w, unsigned char *tag_w,
                     unsigned char *ciphertext_b, int ciphertext_len_b, unsigned char *tag_b,
                     unsigned char *ciphertext_learningRate, int ciphertext_len_learningRate, unsigned char *tag_learningRate,
                     int epochs, int batchSize,
                     unsigned char **ciphertext_w_final, int &ciphertext_len_w_final, unsigned char **tag_w_final,
                     unsigned char **ciphertext_b_final, int &ciphertext_len_b_final, unsigned char **tag_b_final)
{

    // 解密得到明文的数据集X与Y
    vector<vector<int>> X;
    decryptVector(ciphertext_X_train, ciphertext_len_X_train, key, iv, sizeof(iv), tag_X_train, X);
    printf(" SUCCESS X_train 解密成功！\n");
    vector<int> Y;
    decryptVector(ciphertext_Y_train, ciphertext_len_Y_train, key, iv, sizeof(iv), tag_Y_train, Y);
    printf(" SUCCESS Y_train 解密成功！\n");

    // 解密得到明文的模型权重w
    vector<double> w;
    decryptVector(ciphertext_w, ciphertext_len_w, key, iv, sizeof(iv), tag_w, w);
    printf(" SUCCESS w 解密成功！\n");

    // 解密偏置项b
    unsigned char *decryptedtext_b = new unsigned char[sizeof(double)];
    int decrypted_len_b;
    decrypted_len_b = gcm_decrypt(
            ciphertext_b, ciphertext_len_b,
            nullptr, 0,   // 无 AAD
            tag_b,
            key, iv, sizeof(iv),
            decryptedtext_b  // 存储解密后的数据
    );
    double b = bytesToDouble(decryptedtext_b);
    printf(" SUCCESS b 解密成功！\n");

    // 解密学习率learningRate
    unsigned char *decryptedtext_learningRate = new unsigned char[sizeof(double)];
    int decrypted_len_learningRate;
    decrypted_len_learningRate = gcm_decrypt(
            ciphertext_learningRate, ciphertext_len_learningRate,
            nullptr, 0,   // 无 AAD
            tag_learningRate,
            key, iv, sizeof(iv),
            decryptedtext_learningRate  // 存储解密后的数据
    );
    double learningRate = bytesToDouble(decryptedtext_learningRate);
    printf(" SUCCESS learningRate 解密成功！\n");

    int n = X.size();       // 样本数
    int features = X[0].size(); // 特征数

    printf("=== （明文训练）线性回归模型配置 ===\n");
    printf("样本数量: %zu\n", X.size());
    printf("特征数量: %d\n", features);
    printf("学习率: %f\n", learningRate);
    printf("训练轮数: %d\n", epochs);
    printf("批次大小: %d\n", batchSize);

    printf("初始权重: ");
    for (const auto& weight : w) {
        printf("%f ", weight);
    }
    printf("\n初始偏置: %f\n", b);
    printf("=========================\n");

    for (int epoch = 1; epoch <= epochs; epoch++)
    {
        // 遍历所有数据，按批次大小进行训练
        for (int i = 0; i < n; i += batchSize)
        {
            vector<double> dw(features, 0.0);
            double db = 0.0;

            // 处理当前批次的数据
            int currentBatchSize = min(batchSize, n - i); // 防止越界
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
        printf("Epoch %d, Loss: %f, b: %f\n", epoch, loss, b);
        printf("Weights: ");
        for (const auto& weight : w) {
            printf("%f ", weight);
        }
        printf("\n");
    }

    printf("\n训练结束。最终权重和偏置：\n");
    printf("权重: ");
    for (const auto& weight : w) {
        printf("%f ", weight);
    }
    printf("\n偏置: %f\n", b);


    // 加密最终模型权重w
    encryptVector(w, key, iv, sizeof(iv), tag_w_final, ciphertext_w_final, ciphertext_len_w_final);

    // 加密最终模型偏置项b
    unsigned char *plaintext_b = doubleToBytes(b);
    *ciphertext_b_final = new unsigned char[sizeof(double)];
    ciphertext_len_b_final = gcm_encrypt(
            plaintext_b, sizeof(double),
            nullptr, 0,   // 无 AAD
            key, iv, sizeof(iv),
            *ciphertext_b_final, *tag_b_final
    );
}


int main()
{

    double start_time, end_time;

    vector<vector<int>> X_train;  // 多维输入变量
    vector<int> Y_train;           // 预测值

    string train_filePath = "../Data/train_dataset.csv";
    if (loadCSV(train_filePath, X_train, Y_train))
    {
        cout << "训练 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "训练 CSV 文件读取失败！" << '\n';
    }

    vector<vector<int>> X_test;  // 多维输入变量
    vector<int> Y_test;           // 预测值

    string test_filePath = "../Data/test_dataset.csv";
    if (loadCSV(test_filePath, X_test, Y_test))
    {
        cout << "测试 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "测试 CSV 文件读取失败！" << '\n';
    }

    // 初始化 key、iv
    RAND_bytes(key, sizeof(key));
    RAND_bytes(iv, sizeof(iv));

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
    printf("[加密 X_train 与 Y_train ] time is  ------  %f ms\n", (end_time - start_time) * 1000);


    int features = X_train[0].size(); // 获取特征数量
    // 初始化参数
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
    printf("[加密 w、b、learningRate ] time is  ------  %f ms\n", (end_time - start_time) * 1000);

    // 记录返回的w与b
    unsigned char *ciphertext_w_final = nullptr;
    int ciphertext_len_w_final = 0;
    unsigned char *tag_w_final = new unsigned char[16];
    unsigned char *ciphertext_b_final = nullptr;
    int ciphertext_len_b_final = 0;
    unsigned char *tag_b_final = new unsigned char[16];


    start_time = omp_get_wtime();
    // 使用梯度下降进行训练
    gradientDescent(ciphertext_X_train, ciphertext_len_X_train, tag_X_train, // 输入，密文
                    ciphertext_Y_train, ciphertext_len_Y_train, tag_Y_train, // 输入，密文
                    ciphertext_w, ciphertext_len_w, tag_w,// 输入，密文
                    ciphertext_b, ciphertext_len_b, tag_b, // 输入，密文
                    ciphertext_learningRate, ciphertext_len_learningRate, tag_learningRate, // 输入，密文
                    epochs, batchSize, // 输入，密文
                    &ciphertext_w_final, ciphertext_len_w_final, &tag_w_final, // 输出，密文
                    &ciphertext_b_final, ciphertext_len_b_final, &tag_b_final); // 输出，密文
    end_time = omp_get_wtime();
    printf("AES-LinR time is  ------  %f ms\n", (end_time - start_time) * 1000);


    // 解密 w
    vector<double> final_w;
    decryptVector(ciphertext_w_final, ciphertext_len_w_final, key, iv, sizeof(iv), tag_w_final, final_w);
    cout << "final_w 权重: ";
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
    cout << "final_b = " << final_b << '\n';

    // 评估
    evaluate(X_test, Y_test, final_w, final_b);


    printf("Info: successfully returned.\n");
    return 0;

}
