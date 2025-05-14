/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-05
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/



#include "include/protocol.h"
#include "include/fastPai.h"
#include <fstream>
#include <sstream>
#include <gmp.h>
#include <omp.h>
#include <iomanip>
#include <algorithm>
#include <cstdlib>
#include <ctime>

using namespace PROTOCOLSPACE;
using namespace PHESPACE;


long long multiExponent; // 控制最多支持的小数位数
protocol sc;
bool flag;

void multiplication(mpz_t &res, mpz_t &value, mpz_t multiple)
{
    mpz_mul(res, value, multiple);
}

void multiplication(mpz_t &res, mpz_t &value, long long multiple)
{
    mpz_mul_si(res, value, multiple);
}


/**
 * 将value扩大multiple倍存入res
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
        if (flag)
        {
            // 打印 long long int，使用 %lld
            printf("integer_val = %lld\n", integer_val);
            // 打印浮点数（double）
            printf("fractional_val = %.10f\n", fractional_val);
            flag = false;
        }

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
        if (flag)
        {
            // 打印 long long int，使用 %lld
            printf("integer_val = %lld\n", integer_val);
            // 打印浮点数（double）
            printf("fractional_val = %.10f\n", fractional_val);
            flag = false;
        }

        res = integer_val + fractional_val;
        res *= -1;
        mpz_clear(mod);
        mpz_clear(neg_value);
    }
}

void show_dec(mpz_t eres, PaillierThdDec cp, PaillierThdDec csp)
{
    mpz_t neg_one, half_N;
    mpz_inits(neg_one, half_N, NULL);
    mpz_set_si(neg_one, -1);
    mpz_div_ui(half_N, cp.pk.N, 2);

//    gmp_printf("① [原始数据] eres = %Zd\n", eres);
    mpz_t tmp_1, tmp_2, tmp_res;
    mpz_inits(tmp_1, tmp_2, tmp_res, NULL);
    cp.PDec(tmp_1, eres);
    csp.PDec(tmp_2, eres);
    cp.TDec(tmp_res, tmp_1, tmp_2);

    if (mpz_cmp(tmp_res, half_N) > 0)
    {
        mpz_sub(tmp_res, tmp_res, csp.pk.N);
        gmp_printf("② [解密] res = %Zd\n", tmp_res);
    } else
    {
        gmp_printf("② [解密] res = %Zd\n", tmp_res);
    }

    double value;
    division(value, tmp_res, multiExponent);
    printf("④ [转换成浮点数] res = %.10f\n", value);
}


/**
 * 多元线性回归的预测函数：y_hat = w1 * X[i][0] + w2 * X[i][1] + ... + wn * X[i][n-1] + b
 * @param y_hat
 * @param x
 * @param w
 * @param b
 */
void predict(mpz_t y_hat, vector <mpz_t> &x, vector<double> &w, double b)
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


void Cipher_predict(mpz_t ey_hat, vector <mpz_t> &ex, vector <mpz_t> &ew, mpz_t &eb, PaillierThdDec cp, PaillierThdDec csp)
{

    mpz_set(ey_hat, eb);
    for (int j = 0; j < ex.size(); j++)
    {
        mpz_t value;
        mpz_init(value);
        sc.trust_fmul(value, ex[j], ew[j], cp, csp); // 密文乘法
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
double computeLoss(vector <vector<mpz_t>> &X, vector <mpz_t> &Y, vector<double> &w, double b)
{
    double loss = 0.0;
    for (int i = 0; i < X.size(); ++i)
    {
        mpz_t y_hat, y, error;
        mpz_inits(y_hat, y, error, NULL);
        predict(y_hat, X[i], w, b);
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
double computeR2(vector <vector<mpz_t>> &X, vector <mpz_t> &Y, vector<double> &w, double b)
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
        predict(y_hat, X[i], w, b);
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
double computeMAE(vector <vector<mpz_t>> &X, vector <mpz_t> &Y, vector<double> &w, double b)
{
    double mae = 0.0;
    for (int i = 0; i < X.size(); i++)
    {
        mpz_t y_hat, y, error;
        mpz_inits(y_hat, y, error, NULL);

        predict(y_hat, X[i], w, b);
        multiplication(y, Y[i], multiExponent);
        mpz_sub(error, y_hat, y);
        double errorFloat;
        division(errorFloat, error, multiExponent);
        mae += abs(errorFloat);
        mpz_clears(y_hat, y, error, NULL);
    }

    return mae / X.size();
}

void Cipher_gradientDescent(vector <mpz_t> &final_ew, mpz_t &final_eb, // 输出
                            vector <vector<mpz_t>> &X_test, vector <mpz_t> &Y_test, // 用于测试
                            vector <vector<mpz_t>> &eX, vector <mpz_t> &eY,
                            vector <mpz_t> &ew, mpz_t &eb, mpz_t eLearningRate,
                            int epochs, int batchSize,
                            PaillierThdDec cp, PaillierThdDec csp)
{
    int n = eX.size();       // 样本数
    int features = eX[0].size(); // 特征数

    mpz_t neg_one, half_N, multiExponent_mpz, multiExponent_inverse, eMultiExponent_inverse;
    mpz_inits(neg_one, half_N, multiExponent_mpz, multiExponent_inverse, eMultiExponent_inverse, NULL);
    mpz_set_si(multiExponent_mpz, multiExponent);
    mpz_invert(multiExponent_inverse, multiExponent_mpz, csp.pk.N);
//    cp.pai.encrypt(eMultiExponent_inverse, multiExponent_inverse);

    mpz_set_si(neg_one, -1);
    mpz_div_ui(half_N, cp.pk.N, 2);

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
                Cipher_predict(ey_hat, eX[j], ew, eb, cp, csp); // 获取密文预测值
                cp.pai.scl_mul(ey, eY[j], multiExponent * (-1)); // 计算密文误差
                cp.pai.add(eError, ey_hat, ey);
                // 计算损失函数（均分误差）对w的偏导
                for (int k = 0; k < features; ++k)
                {
                    mpz_t res, eres;
                    mpz_inits(res, eres, NULL);

                    sc.trust_fmul(eres, eError, eX[j][k], cp, csp); // 密文乘法
                    cp.pai.add(edw[k], edw[k], eres);
                    mpz_clears(res, eres, NULL);
                }
                // 计算损失函数（均分误差）对b的偏导
                cp.pai.add(edb, edb, eError);
                mpz_clears(ey_hat, ey, eError, NULL);
            }
            
            // 乘以学习率
            for (int k = 0; k < features; k++)
            {
                sc.trust_fmul(edw[k], eLearningRate, edw[k], cp, csp); // 密文乘法
            }
            sc.trust_fmul(edb, eLearningRate, edb, cp, csp); // 密文乘法

            // 归一化处理
            /* 随机生成r1，实现密文除法，针对edb */
            // REE 执行
            mpz_t r1, er1, mask_edb;
            mpz_inits(r1, er1, mask_edb, NULL);
            mpz_urandomb(r1, randstate, sigma);
            multiplication(er1, r1, multiExponent_mpz);
            multiplication(er1, er1, multiExponent_mpz);
            cp.pai.encrypt(er1, er1);
            cp.pai.add(mask_edb, edb, er1);

            mpz_t db_1, db_2, db_add_r1, edb_add_r1_div_n;
            mpz_inits(db_1, db_2, db_add_r1, edb_add_r1_div_n, NULL);
            cp.PDec(db_1, mask_edb);
            // TEE 执行
            csp.PDec(db_2, mask_edb);
            cp.TDec(db_add_r1, db_1, db_2); // 获取明文 db+r
            if (mpz_cmp(db_add_r1, half_N) > 0)
            {
                mpz_sub(db_add_r1, db_add_r1, csp.pk.N);
            }
            mpz_div_ui(edb_add_r1_div_n, db_add_r1, currentBatchSize);
            mpz_div(edb_add_r1_div_n, edb_add_r1_div_n, multiExponent_mpz);
            csp.pai.encrypt(edb_add_r1_div_n, edb_add_r1_div_n);
            // REE 执行
            mpz_t r1_div_n, er1_div_n;
            mpz_inits(r1_div_n, er1_div_n, NULL);
            multiplication(r1_div_n, r1, multiExponent_mpz);
            multiplication(r1_div_n, r1_div_n, multiExponent_mpz);
            mpz_div_ui(r1_div_n, r1_div_n, currentBatchSize);
            mpz_div(r1_div_n, r1_div_n, multiExponent_mpz);
            cp.pai.encrypt(er1_div_n, r1_div_n);
            cp.pai.scl_mul(er1_div_n, er1_div_n, neg_one);
            cp.pai.add(edb, edb_add_r1_div_n, er1_div_n); // 得到密文 [db/n]
            cp.pai.scl_mul(edb, edb, 2);
            mpz_clears(r1, er1, mask_edb, NULL);
            mpz_clears(db_1, db_2, db_add_r1, edb_add_r1_div_n, NULL);
            mpz_clears(r1_div_n, er1_div_n, NULL);

            /* 随机生成r2，实现密文除法，针对edw */
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
                // REE 执行
                mpz_t mask_edw;
                mpz_init(mask_edw);
                cp.pai.add(mask_edw, edw[k], er2);

                mpz_t dw_1, dw_2, dw_add_r2, edw_add_r2_div_n;
                mpz_inits(dw_1, dw_2, dw_add_r2, edw_add_r2_div_n, NULL);
                cp.PDec(dw_1, mask_edw);
                // TEE 执行
                csp.PDec(dw_2, mask_edw);
                cp.TDec(dw_add_r2, dw_1, dw_2); // 获取明文 db+r
                if (mpz_cmp(dw_add_r2, half_N) > 0)
                {
                    mpz_sub(dw_add_r2, dw_add_r2, csp.pk.N);
                }
                mpz_div_ui(edw_add_r2_div_n, dw_add_r2, currentBatchSize);
                mpz_div(edw_add_r2_div_n, edw_add_r2_div_n, multiExponent_mpz);

                csp.pai.encrypt(edw_add_r2_div_n, edw_add_r2_div_n);
                // REE 执行
                cp.pai.add(edw[k], edw_add_r2_div_n, er2_div_n); // 得到密文 [db/n]
                cp.pai.scl_mul(edw[k], edw[k], 2);

                mpz_clear(mask_edw);
                mpz_clears(dw_1, dw_2, dw_add_r2, edw_add_r2_div_n, NULL);

            }
            mpz_clears(r2, er2, r2_div_n, er2_div_n, NULL);

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

        }

        vector<double> w(features);
        double b;

        for (int k = 0; k < features; ++k)
        {
            // 解密
            mpz_t w_1, w_2, w_res;
            mpz_inits(w_1, w_2, w_res, NULL);
            cp.PDec(w_1, ew[k]);
            csp.PDec(w_2, ew[k]);
            cp.TDec(w_res, w_1, w_2);
            if (mpz_cmp(w_res, half_N) > 0)
            {
                mpz_sub(w_res, w_res, csp.pk.N);
            }
            division(w[k], w_res, multiExponent); //缩小multiExponent倍
            mpz_clears(w_1, w_2, w_res, NULL);
        }

        // 解密
        mpz_t b_1, b_2, b_res;
        mpz_inits(b_1, b_2, b_res, NULL);
        cp.PDec(b_1, eb);
        csp.PDec(b_2, eb);
        cp.TDec(b_res, b_1, b_2);
        if (mpz_cmp(b_res, half_N) > 0)
        {
            mpz_sub(b_res, b_res, csp.pk.N);
        }
        division(b, b_res, multiExponent); //缩小multiExponent倍
        mpz_clears(b_1, b_2, b_res, NULL);

        double loss = computeLoss(X_test, Y_test, w, b);
        cout << "Epoch " << epoch << ", Loss: " << loss << ", b: " << b << '\n';
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

void printData(vector <vector<int>> &X, vector<int> &Y)
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

void printData(vector <vector<mpz_t>> &X, vector <mpz_t> &Y)
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


/** 数据集分割函数 **/
void splitDataset(vector <vector<int>> &X, vector<int> &Y,
                  vector <vector<mpz_t>> &X_train, vector <mpz_t> &Y_train,
                  vector <vector<mpz_t>> &X_test, vector <mpz_t> &Y_test, double trainRatio = 0.7)
{
    // 创建索引数组并打乱索引
    vector<int> indices(X.size());
    iota(indices.begin(), indices.end(), 0); // 初始化索引为 0, 1, 2, ..., n-1
    shuffle(indices.begin(), indices.end(), default_random_engine(time(0)));

    // 按比例拆分数据集
    int trainSize = static_cast<size_t>(X.size() * trainRatio);

    for (int i = 0; i < indices.size(); i++)
    {
        // 初始化一行数据
        vector <mpz_t> row(X[indices[i]].size());

        for (int j = 0; j < X[indices[i]].size(); j++)
        {
            mpz_init(row[j]);
            mpz_set_si(row[j], X[indices[i]][j]);
        }

        if (i < trainSize)
        {
            X_train.push_back(move(row));
            mpz_set_si(Y_train[i], Y[indices[i]]);
        } else
        {
            X_test.push_back(move(row));
            mpz_set_si(Y_test[i - trainSize], Y[indices[i]]);
        }
    }
}


void EncDataset(vector <vector<mpz_t>> &eX, vector <mpz_t> &eY, vector <vector<mpz_t>> &X, vector <mpz_t> &Y, Paillier pai)
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

void evaluate(vector <vector<mpz_t>> &X, vector <mpz_t> &Y, vector<double> &w, double &b)
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

// 读取CSV文件并将内容存储到 X 和 Y 中
bool loadCSV(const string &filePath, vector <vector<int>> &X, vector<int> &Y)
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

int main()
{

    double start_time, end_time;

    k = 64;
    sigma = 128;
    setrandom(&randstate);
    KeyGen keyGen(k, sigma);
    Paillier pai(keyGen.pk, keyGen.sk, sigma);

    PaillierThirdParty *psk = keyGen.pai_third_parties;
    PaillierThdDec cp = PaillierThdDec(psk[0].N, psk[0].partial_key, keyGen.pk, sigma);
    PaillierThdDec csp = PaillierThdDec(psk[1].N, psk[1].partial_key, keyGen.pk, sigma);

    vector <vector<int>> X_train_int;  // 多维输入变量
    vector<int> Y_train_int;           // 预测值

    string train_filePath = "../Data/train_dataset.csv";
    if (loadCSV(train_filePath, X_train_int, Y_train_int))
    {
        cout << "训练 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "训练 CSV 文件读取失败！" << '\n';
    }

    vector <vector<int>> X_test_int;  // 多维输入变量
    vector<int> Y_test_int;           // 预测值

    string test_filePath = "../Data/test_dataset.csv";
    if (loadCSV(test_filePath, X_test_int, Y_test_int))
    {
        cout << "测试 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "测试 CSV 文件读取失败！" << '\n';
    }

    // 数据集处理
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


    // 初始化参数
    int features = X_train[0].size(); // 获取特征数量
    vector<double> w(features, 0.0); // 初始化权重为 0
    double b = 0.0;  // 偏置
    double learningRate = 0.0001;  // 学习率
    int epochs = 200;  // 训练轮数
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

    // 密文版本
    start_time = omp_get_wtime();
    Cipher_gradientDescent(final_ew, final_eb, X_test, Y_test, eX_train, eY_train, ew, eb, eLearningRate, epochs, batchSize, cp, csp);
    end_time = omp_get_wtime();
    printf("SEAT LinR time is  ------  %f ms\n", (end_time - start_time) * 1000);
    cout << "------------------------------------------------------\n";


    mpz_t neg_one, half_N;
    mpz_inits(neg_one, half_N, NULL);
    mpz_set_si(neg_one, -1);
    mpz_div_ui(half_N, pai.pk.N, 2);

    vector<double> LinR_w(features);
    double LinR_b;

    // 用户解密final_ew与final_eb
    for (int k = 0; k < features; k++)
    {
        pai.decrypt(final_w[k], final_ew[k]);
        if (mpz_cmp(final_w[k], half_N) > 0)
        {
            mpz_sub(final_w[k], final_w[k], pai.pk.N);
        }
        division(LinR_w[k], final_w[k], multiExponent);
        mpz_clear(final_w[k]);

    }
    pai.decrypt(final_b, final_eb);
    if (mpz_cmp(final_b, half_N) > 0)
    {
        mpz_sub(final_b, final_b, pai.pk.N);
    }
    division(LinR_b, final_b, multiExponent);
    mpz_clear(final_b);

    printf("最终权重: ");
    for (auto &LinR_weight: LinR_w)
    {
        printf("%.10f ", LinR_weight);
    }
    printf("\n最终偏置: %.10f\n", LinR_b);
    evaluate(X_test, Y_test, LinR_w, LinR_b);

    printf("Info: successfully returned.\n");
    return 0;

}
