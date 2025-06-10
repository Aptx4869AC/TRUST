/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-5-27
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <pwd.h>
#include <omp.h>
#include <gmp.h>
#include <vector>
#include <random>
#include <ctime>
#include <cstring>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <netinet/tcp.h>
#include <fcntl.h>
#include <sys/uio.h>
#include <errno.h>
#include <fstream>
#include <sstream>
#include "sgx_urts.h"
#include "App.h"
#include "Enclave_u.h"
#include "fastPai.h"
#include "offline.h"

using namespace std;
using namespace PHESPACE;
using namespace OFFLINESPACE;

sgx_enclave_id_t global_eid = 0;
int sigma;
static ssize_t total_bytes;
long long multiExponent; // 控制最多支持的小数位数
int multi;

typedef struct _sgx_errlist_t
{
    sgx_status_t err;
    const char *msg;
    const char *sug; /* Suggestion */
} sgx_errlist_t;

/* Error code returned by sgx_create_enclave */
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
                SGX_ERROR_NO_DEVICE, "Invalid SGX device.", "Please make sure SGX module is enabled in the BIOS, and install SGX driver afterwards."
        }, {
                SGX_ERROR_MEMORY_MAP_CONFLICT, "Memory map conflicted.", NULL
        }, {
                SGX_ERROR_INVALID_METADATA, "Invalid enclave metadata.", NULL
        }, {
                SGX_ERROR_DEVICE_BUSY, "SGX device was busy.", NULL
        }, {
                SGX_ERROR_INVALID_VERSION, "Enclave version was invalid.", NULL
        }, {
                SGX_ERROR_INVALID_ATTRIBUTE, "Enclave was not authorized.", NULL
        }, {
                SGX_ERROR_ENCLAVE_FILE_ACCESS, "Can't open enclave file.", NULL
        }, {
                SGX_ERROR_MEMORY_MAP_FAILURE, "Failed to reserve memory for the enclave.", NULL
        },};


/* Check error conditions for loading enclave */
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
        printf("Error code is 0x%X. Please refer to the \"Intel SGX SDK Developer Reference\" for more details.\n",
               ret);
}

/* Initialize the enclave:
 *   Call sgx_create_enclave to initialize an enclave instance
 */
int initialize_enclave(void)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;

    /* Call sgx_create_enclave to initialize an enclave instance */
    /* Debug Support: set 2nd parameter to 1 */
    ret = sgx_create_enclave(ENCLAVE_FILENAME, SGX_DEBUG_FLAG, NULL, NULL, &global_eid, NULL);
    if (ret != SGX_SUCCESS)
    {
        print_error_message(ret);
        return -1;
    }

    return 0;
}


/* Other functions */

void multiplication(mpz_t &res, mpz_t &value, mpz_t multiple)
{
    mpz_mul(res, value, multiple);
}

void multiplication(mpz_t &res, mpz_t &value, long long multiple)
{
    mpz_mul_si(res, value, multiple);
}


/**
 * 将value扩大multiple倍存入mpz类型的res
 * @param res
 * @param value
 * @param multiple
 */
void multiplication(mpz_t res, double value, long long multiple)
{
    if (value == 0.0)
    {
        mpz_set_si(res, 0);
        return;
    }
    // 获取整数部分和小数部分
    long long integer = static_cast<long long>(value);  // 转换为整数部分
    double fractional = value - integer;  // 小数部分是原数减去整数部分

    mpz_set_si(res, integer);
    mpz_mul_si(res, res, multiple);
    long long ans = fractional * multiple;

    mpz_t si_tmp;
    mpz_init(si_tmp);
    mpz_set_si(si_tmp, ans);
    mpz_add(res, res, si_tmp);
}

/**
 * 将value缩小multiple并恢复成double类型的res
 * @param res
 * @param value
 * @param multiple
 */
void division(double &res, mpz_t &value, long long multiple)
{
    if (multiple == 0)
    {
        cerr << "错误：除数不能为零！" << '\n';
        return;
    }
    mpz_t integer, remainder, fractional;
    mpz_inits(integer, remainder, fractional, NULL);

    if (mpz_cmp_si(value, 0) > 0)
    {
        mpz_div_ui(integer, value, multiple);

        mpz_t mod;
        mpz_init(mod);
        mpz_init_set_si(mod, multiple);
        mpz_mod(fractional, value, mod);

        // 转换为 int
        long long integer_val = mpz_get_si(integer);
        double fractional_val = mpz_get_si(fractional) / (multiple * 1.0);

        res = integer_val + fractional_val;
        mpz_clear(mod);
    } else
    {
        mpz_t neg_value;
        mpz_init(neg_value);
        mpz_neg(neg_value, value);

        mpz_div_ui(integer, neg_value, multiple);

        mpz_t mod;
        mpz_init(mod);
        mpz_init_set_si(mod, multiple);
        mpz_mod(fractional, neg_value, mod);

        // 转换为 int
        long long integer_val = mpz_get_si(integer);
        double fractional_val = mpz_get_si(fractional) / (multiple * 1.0);

        res = integer_val + fractional_val;
        res *= -1;
        mpz_clear(mod);
        mpz_clear(neg_value);
    }
}


/**
 * 序列化
 * @param num
 * @param buffer
 * @param size
 */
void serialize_mpz(mpz_t num, unsigned char **buffer, size_t *size)
{
    *size = (mpz_sizeinbase(num, 2) + 7) / 8; // 计算所需字节数
    *buffer = new unsigned char[*size];
    mpz_export(*buffer, size, 1, 1, 0, 0, num); // 导出为字节数组
}

/**
 * 反序列化
 * @param res
 * @param buffer
 * @param size
 */
void deserialize_mpz(mpz_t res, unsigned char *buffer, size_t size)
{
    mpz_import(res, size, 1, 1, 0, 0, buffer); // 从字节数组导入
}

void test_serialize_deserialize()
{
    mpz_t original_var, deserialized_var;
    mpz_inits(original_var, deserialized_var, NULL);

    // Assign a value to original_var (replace with your value)
    mpz_set_ui(original_var, 123456789);

    // Serialize original_var
    unsigned char *buffer;
    size_t size;
    serialize_mpz(original_var, &buffer, &size);

    // Deserialize into deserialized_var
    deserialize_mpz(deserialized_var, buffer, size);

    // Print original_var and deserialized_var
    gmp_printf("Original variable: %Zd\n", original_var);
    gmp_printf("Deserialized variable: %Zd\n", deserialized_var);

    // Free memory allocated for buffer
    free(buffer);

    // Clear mpz variables
    mpz_clears(original_var, deserialized_var, NULL);
}


/* OCall functions */
void ocall_print_string(const char *str)
{
    printf("%s", str);
}

/* ECall functions */

void REE_Trust_PDec(mpz_t res, mpz_t eres, PaillierThdDec cp)
{

    mpz_t tmp_1, tmp_2;
    mpz_inits(tmp_1, tmp_2, NULL);
    cp.PDec(tmp_1, eres);
    unsigned char *seal_tmp_2;
    size_t tmp_2_size;
    unsigned char *seal_eres;
    size_t eres_size;
    serialize_mpz(eres, &seal_eres, &eres_size);
    tmp_2_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_tmp_2 = new unsigned char[tmp_2_size];

    Enclave_PDec(global_eid, &seal_tmp_2, &tmp_2_size, &seal_eres, &eres_size);
    deserialize_mpz(tmp_2, seal_tmp_2, tmp_2_size);
    cp.TDec(res, tmp_1, tmp_2);

    free(seal_tmp_2);
    free(seal_eres);
}


/**
 * 密文乘
 * @param res
 * @param em_1
 * @param em_2
 * @param cp
 * @param record_bytes_flag
 */
void REE_Trust_FMUL(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

    // Step 1, REE computes
    mpz_t r, neg_r, er;
    mpz_inits(r, neg_r, er, NULL);
    mpz_set(r, tuple_S0.trust_r);
    cp.pai.add(er, tuple_S0.e_trust_r, tuple_S0.e0); //refresh

    mpz_t X, Y, X_PDec_1;
    mpz_inits(X, Y, X_PDec_1, NULL);
    // X = [m_1 + r]
    cp.pai.add(X, em_1, er);

    // Y = [-r * m_2]
    mpz_neg(neg_r, r);
    cp.pai.scl_mul(Y, em_2, neg_r);

    // PDec(X)
    cp.PDec(X_PDec_1, X);

    /* send <[m_1 + r], (PDec)[m_1 + r], [m_2], [-r * m_2]> to TEE */
    unsigned char *seal_res, *seal_X, *seal_X_PDec_1, *seal_Y, *seal_em_2;
    size_t res_size, X_size, X_PDec_1_size, Y_size, em_2_size;
    serialize_mpz(X, &seal_X, &X_size);
    serialize_mpz(X_PDec_1, &seal_X_PDec_1, &X_PDec_1_size);
    serialize_mpz(Y, &seal_Y, &Y_size);
    serialize_mpz(em_2, &seal_em_2, &em_2_size);
    res_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_res = new unsigned char[res_size];
    Enclave_Trust_FMUL(global_eid, &seal_res, &res_size, &seal_X, &X_size, &seal_X_PDec_1, &X_PDec_1_size,
                       &seal_Y, &Y_size, &seal_em_2, &em_2_size);


    deserialize_mpz(res, seal_res, res_size);

    mpz_clears(r, neg_r, er, X, Y, X_PDec_1, NULL);
    delete[] seal_res;
    delete[] seal_X;
    delete[] seal_X_PDec_1;
    delete[] seal_Y;
    delete[] seal_em_2;
}

/**
 * 密文比较
 * @param res
 * @param em_1
 * @param em_2
 * @param cp
 * @param record_bytes_flag
 */
void REE_Trust_FCMP(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;
    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);

    mpz_t ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1;
    mpz_inits(ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);
    mpz_set(r1, tuple_S0.trust_r1);
    mpz_neg(neg_r1, r1);

    // Step 1
    int pi = (int) random() % 2;
    if (pi == 0)
    {
        cp.pai.scl_mul(em_1r1, em_1, r1);
        cp.pai.scl_mul(neg_em_2r1, em_2, neg_r1);
        cp.pai.add(ed, em_1r1, neg_em_2r1);
        cp.pai.add(e_trust_r1_add_r2, tuple_S0.e_trust_r1_add_r2, tuple_S0.e0);  //refresh
        cp.pai.add(ed, ed, e_trust_r1_add_r2);
    } else
    {
        cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
        cp.pai.scl_mul(em_2r1, em_2, r1);
        cp.pai.add(ed, neg_em_1r1, em_2r1);
        cp.pai.add(e_trust_r2, tuple_S0.e_trust_r2, tuple_S0.e0);  //refresh
        cp.pai.add(ed, ed, e_trust_r2);
    }
    cp.PDec(ed_1, ed);

    mpz_t epi;
    mpz_inits(epi, NULL);
    if (pi == 0)
    {
        cp.pai.add(epi, tuple_S0.e0, tuple_S0.e0);   //refresh

    } else
    {
        cp.pai.add(epi, tuple_S0.e1, tuple_S0.e0);   //refresh
    }


    /* send <[d], (PDec)[d], [pai]> to TEE */
    unsigned char *seal_res, *seal_ed, *seal_ed_1, *seal_epi;
    size_t res_size, ed_size, ed_1_size, epi_size;
    serialize_mpz(ed, &seal_ed, &ed_size);
    serialize_mpz(ed_1, &seal_ed_1, &ed_1_size);
    serialize_mpz(epi, &seal_epi, &epi_size);
    res_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_res = new unsigned char[res_size];
    Enclave_Trust_FCMP(global_eid, &seal_res, &res_size, &seal_ed, &ed_size, &seal_ed_1, &ed_1_size,
                       &seal_epi, &epi_size);


    deserialize_mpz(res, seal_res, res_size);


    mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, epi, NULL);
    delete[] seal_res;
    delete[] seal_ed;
    delete[] seal_ed_1;
    delete[] seal_epi;
}


void EncDataset(vector<vector<mpz_t>> &eX, vector<mpz_t> &eY, vector<vector<mpz_t>> &X, vector<mpz_t> &Y, Paillier pai)
{
    for (int i = 0; i < X.size(); i++)
    {

        vector <mpz_t> row(X[i].size());

        for (int j = 0; j < X[i].size(); j++)
        {
            mpz_init(row[j]);
            mpz_set(row[j], X[i][j]);
            pai.encrypt(row[j], row[j]);
        }

        eX.push_back(move(row));
        mpz_set(eY[i], Y[i]);
        pai.encrypt(eY[i], Y[i]);
    }
}

void show_dec(mpz_t eres, PaillierThdDec cp, Paillier pai)
{
    mpz_t neg_one, half_N;
    mpz_inits(neg_one, half_N, NULL);
    mpz_set_si(neg_one, -1);
    mpz_div_ui(half_N, cp.pk.N, 2);
    gmp_printf("① [原始数据] eres = %Zd\n", eres);

    mpz_t tmp_dec;
    mpz_init(tmp_dec);
    pai.decrypt(tmp_dec, eres);
    gmp_printf("解密后 tmp_dec = %Zd\n", tmp_dec);
    mpz_clear(tmp_dec);

    mpz_t tmp_res;
    mpz_init(tmp_res);
    REE_Trust_PDec(tmp_res, eres, cp);
    if (mpz_cmp(tmp_res, half_N) > 0)
    {
        mpz_sub(tmp_res, tmp_res, cp.pk.N);
        gmp_printf("② [解密] res = %Zd\n", tmp_res);
    } else
    {
        gmp_printf("② [解密] res = %Zd\n", tmp_res);
    }

    double value;
    division(value, tmp_res, multiExponent);
    printf("③ [转换成浮点数] res = %.10f\n", value);

    mpz_clear(tmp_res);
}

/**
 * 多元线性回归的预测函数：y_hat = w1 * X[i][0] + w2 * X[i][1] + ... + wn * X[i][n-1] + b
 * @param y_hat
 * @param x
 * @param w
 * @param b
 */
void predict_PlaintextVersion(mpz_t y_hat, vector<mpz_t> &x, vector<double> &w, double b)
{
    multiplication(y_hat, b, multiExponent);
    for (int j = 0; j < x.size(); j++)
    {
        mpz_t value;
        mpz_init(value);
        multiplication(value, w[j], multiExponent);
        mpz_mul(value, x[j], value); // w[j] * x[j]
        mpz_add(y_hat, y_hat, value);
        mpz_clear(value);
    }
}


void predict_CipherVersion(mpz_t ey_hat, vector<mpz_t> &ex, vector<mpz_t> &ew, mpz_t &eb, PaillierThdDec cp)
{

    mpz_set(ey_hat, eb);
    for (int j = 0; j < ex.size(); j++)
    {
        mpz_t value;
        mpz_init(value);
        REE_Trust_FMUL(value, ex[j], ew[j], cp);
        cp.pai.add(ey_hat, ey_hat, value);
        mpz_clear(value);
    }
}

/**
 * 计算误差损失函数 (Mean Squared Error)
 * @param X
 * @param Y
 * @param w
 * @param b
 * @return
 */
double computeLoss(vector<vector<mpz_t>> &X, vector<mpz_t> &Y, vector<double> &w, double b)
{
    double loss = 0.0;
    for (int i = 0; i < X.size(); ++i)
    {
        mpz_t y_hat, y, error;
        mpz_inits(y_hat, y, error, NULL);
        predict_PlaintextVersion(y_hat, X[i], w, b);
        multiplication(y, Y[i], multiExponent);
        mpz_sub(error, y_hat, y);
        double errorFloat;
        division(errorFloat, error, multiExponent);
        loss += pow(errorFloat, 2);
    }
    return loss / X.size();
}

/**
 * 计算 R^2 (决定系数)
 * @param X
 * @param Y
 * @param w
 * @param b
 * @return
 */
double computeR2(vector<vector<mpz_t>> &X, vector<mpz_t> &Y, vector<double> &w, double b)
{
    double ss_tot = 0.0, ss_res = 0.0;

    mpz_t sum;
    mpz_init(sum);
    mpz_set_si(sum, 0);
    for (int i = 0; i < Y.size(); i++)
    {
        mpz_add(sum, sum, Y[i]);
    }
    double mean_y = mpz_get_si(sum) * 1.0;
    mean_y /= Y.size();
    mpz_clear(sum);

    for (int i = 0; i < X.size(); i++)
    {
        mpz_t y_hat, y, error;
        mpz_inits(y_hat, y, error, NULL);
        predict_PlaintextVersion(y_hat, X[i], w, b);
        multiplication(y, Y[i], multiExponent);
        mpz_sub(error, y_hat, y);
        double errorFloat;
        division(errorFloat, error, multiExponent);
        ss_res += pow(errorFloat, 2);
        ss_tot += pow(mpz_get_si(Y[i]) - mean_y, 2);
        mpz_clears(y_hat, y, error, NULL);
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
double computeMAE(vector<vector<mpz_t>> &X, vector<mpz_t> &Y, vector<double> &w, double b)
{
    double mae = 0.0;
    for (int i = 0; i < X.size(); i++)
    {
        mpz_t y_hat, y, error;
        mpz_inits(y_hat, y, error, NULL);
        predict_PlaintextVersion(y_hat, X[i], w, b);
        multiplication(y, Y[i], multiExponent);
        mpz_sub(error, y_hat, y);
        double errorFloat;
        division(errorFloat, error, multiExponent);
        mae += abs(errorFloat);
        mpz_clears(y_hat, y, error, NULL);
    }

    return mae / X.size();
}

void evaluate(vector<vector<mpz_t>> &X, vector<mpz_t> &Y, vector<double> &w, double &b)
{
    // 计算评价指标
    double mse = computeLoss(X, Y, w, b);
    double r2 = computeR2(X, Y, w, b);
    double mae = computeMAE(X, Y, w, b);

    cout << "------------------------------------------------------\n";
    printf("MSE: %.4f\n", mse / (multi * multi));
    printf("R²: %.4f\n", r2);
    printf("MAE: %.4f\n", mae / (multi));
    cout << "------------------------------------------------------\n";

}


void gradientDescent(vector<mpz_t> &final_ew, mpz_t &final_eb, // 输出
                     vector<vector<mpz_t>> &X_test, vector<mpz_t> &Y_test, // 用于测试
                     vector<vector<mpz_t>> &eX, vector<mpz_t> &eY,
                     vector<mpz_t> &ew, mpz_t &eb, mpz_t eLearningRate,
                     int epochs, int batchSize,
                     PaillierThdDec cp, Paillier pai)
{
    int n = eX.size();       // 样本数
    int features = eX[0].size(); // 特征数

    mpz_t neg_one, half_N, multiExponent_mpz;
    mpz_inits(neg_one, half_N, multiExponent_mpz, NULL);
    mpz_set_si(neg_one, -1);
    mpz_div_ui(half_N, cp.pk.N, 2);
    mpz_set_si(multiExponent_mpz, multiExponent);

    for (int epoch = 1; epoch <= epochs; epoch++)
    {

        // 遍历所有数据，按批次大小进行训练
        for (int i = 0; i < n; i += batchSize)
        {

            vector <mpz_t> edw(features);
            for (int index = 0; index < edw.size(); index++)
            {
                mpz_init(edw[index]);
                mpz_set_si(edw[index], 0);
                cp.pai.encrypt(edw[index], edw[index]);
            }
            mpz_t edb;
            mpz_init(edb);
            mpz_set_si(edb, 0);
            cp.pai.encrypt(edb, edb);

            // 处理当前批次的数据
            int currentBatchSize = min(batchSize, n - i); // 防止越界
            for (int j = i; j < i + currentBatchSize; j++)
            {
                mpz_t ey_hat, ey, eError;
                mpz_inits(ey_hat, ey, eError, NULL);
                predict_CipherVersion(ey_hat, eX[j], ew, eb, cp); // 获取密文预测值
                cp.pai.scl_mul(ey, eY[j], multiExponent * (-1)); // 计算密文误差
                cp.pai.add(eError, ey_hat, ey);

                // 计算损失函数（均分误差）对w的偏导
                for (int k = 0; k < features; ++k)
                {
                    mpz_t res, eres;
                    mpz_inits(res, eres, NULL);
                    REE_Trust_FMUL(eres, eError, eX[j][k], cp);
                    cp.pai.add(edw[k], edw[k], eres);
                    mpz_clears(res, eres, NULL);
                }
                // 计算损失函数（均分误差）对b的偏导
                cp.pai.add(edb, edb, eError);
                mpz_clears(ey_hat, ey, eError, NULL);
            }
//            for (int k = 0; k < features; ++k)
//            {
//                show_dec(edw[k], cp, pai);
//            }
//            cout << "------\n";
//            show_dec(edb, cp, pai);
//            cout << "-------------------------------------------------\n";


            // 乘以学习率
            for (int k = 0; k < features; k++)
            {
                REE_Trust_FMUL(edw[k], eLearningRate, edw[k], cp);
            }
            REE_Trust_FMUL(edb, eLearningRate, edb, cp);
//            for (int k = 0; k < features; ++k)
//            {
//                show_dec(edw[k], cp, pai);
//            }
//            cout << "------\n";
//            show_dec(edb, cp, pai);
//            cout << "-------------------------------------------------\n";


            // 归一化处理
            // Step-1：REE 执行
            /* 随机生成r1，用于后续实现密文除法，针对edb */
            mpz_t r1, er1;
            mpz_inits(r1, er1, NULL);
            mpz_urandomb(r1, randstate, sigma);
            multiplication(er1, r1, multiExponent_mpz);
            multiplication(er1, er1, multiExponent_mpz);
            cp.pai.encrypt(er1, er1);

            mpz_t r1_div_n, er1_div_n;
            mpz_inits(r1_div_n, er1_div_n, NULL);
            multiplication(r1_div_n, r1, multiExponent_mpz);
            multiplication(r1_div_n, r1_div_n, multiExponent_mpz);
            mpz_div_ui(r1_div_n, r1_div_n, currentBatchSize);
            mpz_div(r1_div_n, r1_div_n, multiExponent_mpz);
            cp.pai.encrypt(er1_div_n, r1_div_n);
            cp.pai.scl_mul(er1_div_n, er1_div_n, neg_one);

            mpz_t mask_edb;
            mpz_init(mask_edb);
            cp.pai.add(mask_edb, edb, er1);
            mpz_t mask_edb_1, edb_add_r1_div_n;
            mpz_inits(mask_edb_1, edb_add_r1_div_n, NULL);
            cp.PDec(mask_edb_1, mask_edb);

            /* send <db_1, mask_edb, currentBatchSize, multiExponent> to TEE */
            unsigned char *seal_edb_add_r1_div_n, *seal_Mask_edb, *seal_Mask_edb_1;
            size_t edb_add_r1_div_n_size, Mask_edb_size, Mask_edb_1_size;
            serialize_mpz(mask_edb, &seal_Mask_edb, &Mask_edb_size);
            serialize_mpz(mask_edb_1, &seal_Mask_edb_1, &Mask_edb_1_size);
            edb_add_r1_div_n_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
            seal_edb_add_r1_div_n = new unsigned char[edb_add_r1_div_n_size];
            size_t size_currentBatchSize = static_cast<size_t>(currentBatchSize);
            size_t size_multiExponent = static_cast<size_t>(multiExponent);
            Enclave_Trust_Mask(global_eid, &seal_edb_add_r1_div_n, &edb_add_r1_div_n_size, &seal_Mask_edb_1, &Mask_edb_1_size, &seal_Mask_edb, &Mask_edb_size,
                               &size_currentBatchSize, &size_multiExponent);

            deserialize_mpz(edb_add_r1_div_n, seal_edb_add_r1_div_n, edb_add_r1_div_n_size);
            // Step-3：REE 执行
            cp.pai.add(edb, edb_add_r1_div_n, er1_div_n); // 得到密文 [db/n]
            cp.pai.scl_mul(edb, edb, 2);
            mpz_clears(mask_edb, mask_edb_1, edb_add_r1_div_n, NULL);
            delete[] seal_edb_add_r1_div_n;
            delete[] seal_Mask_edb;
            delete[] seal_Mask_edb_1;

            /* 随机生成r2，用于后续实现密文除法，针对edw */
            // Step-1：REE 执行
            mpz_t r2, er2;
            mpz_inits(r2, er2, NULL);
            mpz_urandomb(r2, randstate, sigma);
            multiplication(er2, r2, multiExponent_mpz);
            multiplication(er2, er2, multiExponent_mpz);
            cp.pai.encrypt(er2, er2);

            mpz_t r2_div_n, er2_div_n;
            mpz_inits(r2_div_n, er2_div_n, NULL);
            multiplication(r2_div_n, r2, multiExponent_mpz);
            multiplication(r2_div_n, r2_div_n, multiExponent_mpz);
            mpz_div_ui(r2_div_n, r2_div_n, currentBatchSize);
            mpz_div(r2_div_n, r2_div_n, multiExponent_mpz);
            cp.pai.encrypt(er2_div_n, r2_div_n);
            cp.pai.scl_mul(er2_div_n, er2_div_n, neg_one);
            for (int k = 0; k < features; ++k)
            {
                mpz_t mask_edw;
                mpz_init(mask_edw);
                cp.pai.add(mask_edw, edw[k], er2);
                mpz_t mask_edw_1, edw_add_r2_div_n;
                mpz_inits(mask_edw_1, edw_add_r2_div_n, NULL);
                cp.PDec(mask_edw_1, mask_edw);

                /* send <dw_1, mask_edw, currentBatchSize, multiExponent> to TEE */
                unsigned char *seal_edw_add_r2_div_n, *seal_Mask_edw, *seal_Mask_edw_1;
                size_t edw_add_r2_div_n_size, Mask_edw_size, Mask_edw_1_size;
                serialize_mpz(mask_edw, &seal_Mask_edw, &Mask_edw_size);
                serialize_mpz(mask_edw_1, &seal_Mask_edw_1, &Mask_edw_1_size);
                edw_add_r2_div_n_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
                seal_edw_add_r2_div_n = new unsigned char[edw_add_r2_div_n_size];
                Enclave_Trust_Mask(global_eid, &seal_edw_add_r2_div_n, &edw_add_r2_div_n_size, &seal_Mask_edw_1, &Mask_edw_1_size, &seal_Mask_edw, &Mask_edw_size,
                                   &size_currentBatchSize, &size_multiExponent);

                deserialize_mpz(edw_add_r2_div_n, seal_edw_add_r2_div_n, edw_add_r2_div_n_size);
                // Step-3：REE 执行
                cp.pai.add(edw[k], edw_add_r2_div_n, er2_div_n); // 得到密文 [db/n]
                cp.pai.scl_mul(edw[k], edw[k], 2);
                mpz_clears(mask_edw, mask_edw_1, edw_add_r2_div_n, NULL);

                delete[] seal_edw_add_r2_div_n;
                delete[] seal_Mask_edw;
                delete[] seal_Mask_edw_1;
            }


            mpz_t eValue;
            mpz_init(eValue);
            // 更新 w
            for (int k = 0; k < features; k++)
            {
                cp.pai.scl_mul(eValue, edw[k], neg_one);
                cp.pai.add(ew[k], ew[k], eValue);
            }
            // 更新 b
            cp.pai.scl_mul(eValue, edb, neg_one);
            cp.pai.add(eb, eb, eValue);


            for (int k = 0; k < features; ++k)
            {
                mpz_clear(edw[k]);
            }
            mpz_clear(edb);
            mpz_clear(eValue);
            mpz_clears(r1, er1, r1_div_n, er1_div_n, NULL);
            mpz_clears(r2, er2, r2_div_n, er2_div_n, NULL);

        }


        // 测试输出看看效果
        vector<double> w(features);
        double b;

        for (int k = 0; k < features; ++k)
        {
            // 解密
            mpz_t w_res;
            mpz_init(w_res);
            REE_Trust_PDec(w_res, ew[k], cp);
            if (mpz_cmp(w_res, half_N) > 0)
            {
                mpz_sub(w_res, w_res, cp.pk.N);
            }
            division(w[k], w_res, multiExponent); //缩小multiExponent倍
            mpz_clear(w_res);
        }

        // 解密
        mpz_t b_res;
        mpz_init(b_res);
        REE_Trust_PDec(b_res, eb, cp);
        if (mpz_cmp(b_res, half_N) > 0)
        {
            mpz_sub(b_res, b_res, cp.pk.N);
        }
        division(b, b_res, multiExponent); //缩小multiExponent倍
        mpz_clear(b_res);

        double loss = computeLoss(X_test, Y_test, w, b);
        cout << "Epoch " << epoch << ", Loss: " << loss / (multi * multi) << ", b: " << b << '\n';
        cout << "Weights: ";
        for (double weight: w)
        {
            cout << weight << " ";
        }
        cout << '\n';
        printf("epoch = %d is OK!\n", epoch);
        cout << "------------------------------------------------------------------\n";
    }

    // 训练结束赋值
    for (int k = 0; k < features; ++k)
    {
        mpz_set(final_ew[k], ew[k]);
    }
    mpz_set(final_eb, eb);

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

void printData(vector<vector<mpz_t>> &X, vector<mpz_t> &Y)
{
    cout << "读取的数据：" << '\n';
    size_t sampleIndex = 1;

    for (auto &row: X)  // 遍历每个样本
    {
        cout << "样本 " << sampleIndex++ << ": ";
        for (auto &value: row)  // 遍历样本中的每个特征
        {
            gmp_printf("%Zd ", value);
        }

        cout << "| 真实值: ";
        gmp_printf("%Zd\n", Y[sampleIndex - 2]);  // 打印预测值
    }
}

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

/* Application entry */
int main(int argc, char *argv[])
{


    /* Initialize the enclave */
    if (initialize_enclave() < 0)
    {
        printf("Enter a character before exit ...\n");
        getchar();
        return -1;
    }

    double start_time, end_time;

    int k = 64; //BIT SECURITY, 64、80、104、112、124、128
    sigma = 128;
    switch (k)
    {
        case 64:
            printf("N 的比特长度为 512 bits\n");
            break;
        case 80:
            printf("N 的比特长度为 1024 bits\n");
            break;
        case 104:
            printf("N 的比特长度为 1536 bits\n");
            break;
        case 112:
            printf("N 的比特长度为 2048 bits\n");
            break;
        case 124:
            printf("N 的比特长度为 2560 bits\n");
            break;
        case 128:
            printf("N 的比特长度为 3072 bits\n");
            break;
        case 192:
            printf("N 的比特长度为 7680 bits\n");
            break;
    }
    setrandom(&randstate);
    KeyGen keyGen(k, sigma);

    Paillier pai(keyGen.pk, keyGen.sk, sigma);
    PaillierThirdParty *psk = keyGen.pai_third_parties;
    PaillierThdDec cp = PaillierThdDec(psk[0].N, psk[0].partial_key, keyGen.pk, sigma);
    printf("REE 初始化成功\n");

    // 序列化TEE部分的私钥内容
    PaillierThdDec enclave = PaillierThdDec(psk[1].N, psk[1].partial_key, keyGen.pk, sigma);
    unsigned char *seal_partial_key, *seal_n, *seal_h_N;
    size_t partial_key_size, n_size, h_N_size;
    serialize_mpz(enclave.partial_key, &seal_partial_key, &partial_key_size);
    serialize_mpz(enclave.N, &seal_n, &n_size);
    serialize_mpz(enclave.pk.h_N, &seal_h_N, &h_N_size);
    Enclave_Init(global_eid, &seal_partial_key, &partial_key_size, &seal_n, &n_size, &seal_h_N, &h_N_size,
                 enclave.pk.L_k);

    vector <vector<int>> X_train_int;  // 多维输入变量
    vector<int> Y_train_int;           // 预测值

    string train_filePath = "./Data/train_dataset_student.csv";
    if (loadCSV(train_filePath, X_train_int, Y_train_int))
    {
        cout << "训练 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "训练 CSV 文件读取失败！" << '\n';
    }

    vector <vector<int>> X_test_int;  // 多维输入变量
    vector<int> Y_test_int;           // 预测值

    string test_filePath = "./Data/test_dataset_student.csv";
    if (loadCSV(test_filePath, X_test_int, Y_test_int))
    {
        cout << "测试 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "测试 CSV 文件读取失败！" << '\n';
    }

    // 数据集精度处理
    multiExponent = pow(10, 10);

    // 初始化训练集和测试集
    vector <vector<mpz_t>> X_train, X_test;
    vector <mpz_t> Y_train(Y_train_int.size()), Y_test(Y_test_int.size());

    // X_train
    for (int i = 0; i < X_train_int.size(); i++)
    {
        vector <mpz_t> row(X_train_int[i].size());
        for (int j = 0; j < X_train_int[i].size(); j++)
        {
            mpz_init(row[j]);
            mpz_set_si(row[j], X_train_int[i][j]);
        }
        X_train.push_back(move(row));
    }

    // X_test
    for (int i = 0; i < X_test_int.size(); i++)
    {
        vector <mpz_t> row(X_test_int[i].size());
        for (int j = 0; j < X_test_int[i].size(); j++)
        {
            mpz_init(row[j]);
            mpz_set_si(row[j], X_test_int[i][j]);
        }
        X_test.push_back(move(row));
    }

    // Y_train
    for (int i = 0; i < Y_train_int.size(); i++)
    {
        mpz_init(Y_train[i]);
        mpz_set_si(Y_train[i], Y_train_int[i]);
    }

    // Y_test
    for (int i = 0; i < Y_test_int.size(); i++)
    {
        mpz_init(Y_test[i]);
        mpz_set_si(Y_test[i], Y_test_int[i]);
    }
    vector <vector<mpz_t>> eX_train, eX_test;
    vector <mpz_t> eY_train(Y_train.size()), eY_test(Y_test.size());

    // 加密 X_train 与 Y_train
    start_time = omp_get_wtime();
    EncDataset(eX_train, eY_train, X_train, Y_train, pai);
    end_time = omp_get_wtime();
    printf("[加密 X_train 与 Y_train ] time is  ------  %f ms\n", (end_time - start_time) * 1000);

    // 加密 X_test 与 Y_test
    start_time = omp_get_wtime();
    EncDataset(eX_test, eY_test, X_test, Y_test, pai);
    end_time = omp_get_wtime();
    printf("[加密 X_test 与 Y_test ] time is  ------  %f ms\n", (end_time - start_time) * 1000);

    multi = 1;

    // 初始化参数
    int features = X_train[0].size(); // 获取特征数量
    vector<double> w(features, 0.0); // 初始化权重为 0
    double b = 0.0;  // 偏置
    double learningRate = 0.0001 / (multi * multi);  // 学习率
    int epochs = atoi(argv[1]);
//    int epochs = 100;
    int batchSize = 64; // 批次大小


    printf("=== （密文训练）线性回归模型配置 ===\n");
    printf("样本数量: %ld\n", X_train.size());
    printf("特征数量: %d\n", features);
    printf("学习率: %.10f\n", learningRate);
    printf("训练轮数: %d\n", epochs);
    printf("批次大小: %d\n", batchSize);
    printf("初始权重: ");
    for (const auto &weight: w)
    {
        printf("%.10f ", weight);
    }
    printf("\n初始偏置: %f\n", b);
    printf("=========================\n");


    vector <mpz_t> ew(w.size());
    for (int j = 0; j < ew.size(); j++)
    {
        mpz_init(ew[j]);
        multiplication(ew[j], w[j], multiExponent);
    }
    mpz_t eLearningRate, eb;
    mpz_inits(eLearningRate, eb, NULL);
    multiplication(eLearningRate, learningRate, multiExponent);
    multiplication(eb, b, multiExponent);

    // 加密 learningRate, w, b
    start_time = omp_get_wtime();
    for (int k = 0; k < features; k++)
    {
        pai.encrypt(ew[k], ew[k]);
    }
    pai.encrypt(eLearningRate, eLearningRate);
    pai.encrypt(eb, eb);
    end_time = omp_get_wtime();
    printf("[加密 learningRate, w, b ] time is  ------  %f ms\n", (end_time - start_time) * 1000);

    vector <mpz_t> final_ew(w.size()), final_w(w.size());
    for (int k = 0; k < features; k++)
    {
        mpz_init(final_ew[k]);
        mpz_init(final_w[k]);
    }
    mpz_t final_eb, final_b;
    mpz_inits(final_eb, final_b, NULL);

    // SEAT
    start_time = omp_get_wtime();
    gradientDescent(final_ew, final_eb, X_train, Y_train, eX_train, eY_train, ew, eb, eLearningRate, epochs, batchSize, cp, pai);
    end_time = omp_get_wtime();
    printf("SEAT LinR time is  ------  %f ms\n", (end_time - start_time) * 1000);
    cout << "------------------------------------------------------\n";


    mpz_t half_N;
    mpz_init(half_N);
    mpz_div_ui(half_N, cp.pk.N, 2);

    vector<double> LinR_w(features);
    double LinR_b;

    // 用户解密final_ew与final_eb
    for (int k = 0; k < features; k++)
    {
        pai.decrypt(final_w[k], final_ew[k]);
        if (mpz_cmp(final_w[k], half_N) > 0)
        {
            mpz_sub(final_w[k], final_w[k], cp.pk.N);
        }
        gmp_printf("final_w[%d] = %Zd\n", k, final_w[k]);
        division(LinR_w[k], final_w[k], multiExponent);
        mpz_clear(final_w[k]);

    }
    pai.decrypt(final_b, final_eb);
    if (mpz_cmp(final_b, half_N) > 0)
    {
        mpz_sub(final_b, final_b, cp.pk.N);
    }
    gmp_printf("final_b = %Zd\n", final_b);
    division(LinR_b, final_b, multiExponent);
    mpz_clear(final_b);

    cout << "最终权重: ";
    for (auto &LinR_weight: LinR_w)
    {
        printf("%.10f ", LinR_weight);
    }
    printf("\n最终偏置: %.10f\n", LinR_b);
    evaluate(X_test, Y_test, LinR_w, LinR_b);// 模型评价


    Enclave_Release(global_eid);
    /* Destroy the enclave */
    sgx_destroy_enclave(global_eid);
    printf("Info: successfully returned.\n");
    return 0;
}

