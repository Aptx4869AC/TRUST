/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-01
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include <iostream>
#include <vector>
#include <fstream>
#include <sstream>
#include <gmp.h>
#include <omp.h>
#include <algorithm>
#include <cstdlib>
#include <ctime>
#include <random>
#include <cstring>
#include <omp.h>

using namespace std;

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
double computeLoss(const vector <vector<int>> &X, const vector<int> &Y, const vector<double> &w, double b)
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
double computeR2(const vector <vector<int>> &X, const vector<int> &Y, const vector<double> &w, double b)
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
double computeMAE(const vector <vector<int>> &X, const vector<int> &Y, const vector<double> &w, double b)
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
 * 梯度下降更新 w 和 b 参数
 * @param X
 * @param Y
 * @param w
 * @param b
 * @param learningRate
 * @param epochs
 * @param batchSize
 */
void gradientDescent(const vector <vector<int>> &X, const vector<int> &Y, vector<double> &w, double &b,
                     double learningRate, int epochs, int batchSize)
{
    int n = X.size();       // 样本数
    int features = X[0].size(); // 特征数

    printf("=== （明文训练）线性回归模型配置 ===\n");
    printf("样本数量: %zu\n", X.size());
    printf("特征数量: %d\n", features);
    printf("学习率: %f\n", learningRate);
    printf("训练轮数: %d\n", epochs);
    printf("批次大小: %d\n", batchSize);

    printf("初始权重: ");
    for (const auto &weight: w)
    {
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
        for (const auto &weight: w)
        {
            printf("%f ", weight);
        }
        printf("\n");
        cout << "------------------------------------------------------------------\n";


    }

    printf("\n训练结束。最终权重和偏置：\n");
    printf("权重: ");
    for (const auto &weight: w)
    {
        printf("%f ", weight);
    }
    printf("\n偏置: %f\n", b);
}


/**
 * 加载CSV文件的函数
 * @param filePath
 * @param X
 * @param Y
 * @return
 */
bool loadCSV(const string &filePath, vector <vector<double>> &X, vector<double> &Y)
{
    ifstream file(filePath);
    if (!file.is_open())
    {
        return false;
    }

    string line;
    while (getline(file, line))
    {
        stringstream ss(line);
        vector<double> row;
        string value;
        double target;

        // 读取每一行中的特征值
        while (getline(ss, value, ','))
        {
            try
            {
                row.push_back(stod(value));
            } catch (const invalid_argument &)
            {
                cerr << "无效的数值: " << value << '\n';
                return false;
            }
        }

        // 将最后一个值作为预测目标值
        target = row.back();
        row.pop_back(); // 移除最后一个值

        X.push_back(row);
        Y.push_back(target);
    }

    file.close();
    return true;
}

/**
 * 打印数据的函数
 * @param X
 * @param Y
 */
void printData(const vector <vector<int>> &X, const vector<int> &Y)
{
    cout << "读取的数据：" << '\n';
    for (size_t i = 0; i < X.size(); ++i)
    {
        cout << "样本 " << i + 1 << ": ";
        for (size_t j = 0; j < X[i].size(); ++j)
        {
            cout << X[i][j] << " ";
        }
        cout << "| 预测值: " << Y[i] << '\n';
    }
}

/**
 * 数据集分割函数
 * @param X
 * @param Y
 * @param X_train
 * @param Y_train
 * @param X_test
 * @param Y_test
 * @param trainRatio
 */
void splitDataset(const vector <vector<int>> &X, const vector<int> &Y,
                  vector <vector<int>> &X_train, vector<int> &Y_train,
                  vector <vector<int>> &X_test, vector<int> &Y_test, double trainRatio = 0.7)
{
    // 创建索引数组并打乱索引
    vector<int> indices(X.size());
    iota(indices.begin(), indices.end(), 0); // 初始化索引为 0, 1, 2, ..., n-1
    shuffle(indices.begin(), indices.end(), default_random_engine(time(0)));

    // 按比例拆分数据集
    size_t trainSize = static_cast<size_t>(X.size() * trainRatio);

    for (size_t i = 0; i < indices.size(); ++i)
    {
        if (i < trainSize)
        {
            X_train.push_back(X[indices[i]]);
            Y_train.push_back(Y[indices[i]]);
        } else
        {
            X_test.push_back(X[indices[i]]);
            Y_test.push_back(Y[indices[i]]);
        }
    }
}

/**
 * 保存数据到CSV文件
 * @param filename
 * @param X
 * @param Y
 */
void saveToCSV(const string &filename, const vector <vector<int>> &X, const vector<int> &Y)
{
    ofstream file(filename, ios::out | ios::trunc);  // 创建空文件或清空旧文件
    if (file.is_open())
    {
        cout << "已创建文件: " << filename << '\n';

        // 写入标题（可选）
        file << "Feature1,Feature2,...,FeatureN,Label\n";  // 根据实际特征数调整

        // 写入数据
        for (size_t i = 0; i < X.size(); ++i)
        {
            for (size_t j = 0; j < X[i].size(); ++j)
            {
                file << X[i][j];
                if (j < X[i].size() - 1) file << ",";
            }
            file << "," << Y[i] << "\n";
        }
        file.close();
        cout << "保存文件成功: " << filename << '\n';
    } else
    {
        cerr << "无法打开文件: " << filename << '\n';
    }
}

/**
 * 读取CSV文件并将内容存储到 X 和 Y 中
 * @param filePath
 * @param X
 * @param Y
 * @return
 */
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

/**
 * 评价函数
 * @param X
 * @param Y
 * @param w
 * @param b
 */
void evaluate(vector <vector<int>> &X, vector<int> &Y, vector<double> &w, double &b)
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


int main()
{

    double start_time, end_time;


    vector <vector<int>> X_train;  // 多维输入变量
    vector<int> Y_train;           // 预测值

    string train_filePath = "../Data/Student_Performance/train_dataset_student.csv";
    if (loadCSV(train_filePath, X_train, Y_train))
    {
        cout << "训练 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "训练 CSV 文件读取失败！" << '\n';
    }

    vector <vector<int>> X_test;  // 多维输入变量
    vector<int> Y_test;           // 预测值

    string test_filePath = "../Data/Student_Performance/test_dataset_student.csv";
    if (loadCSV(test_filePath, X_test, Y_test))
    {
        cout << "测试 CSV 文件读取成功！" << '\n';
    } else
    {
        cerr << "测试 CSV 文件读取失败！" << '\n';
    }


    int features = X_train[0].size(); // 获取特征数量

    // 初始化参数
    vector<double> w(features, 0.0); // 初始化权重为 0
    double b = 0.0;  // 偏置
    double learningRate = 0.0001;  // 学习率
    int epochs = 800;  // 训练轮数
    int batchSize = 64; // 批次大小

    // 明文版本
    start_time = omp_get_wtime();
    gradientDescent(X_train, Y_train, w, b, learningRate, epochs, batchSize);    // 使用梯度下降进行训练
    end_time = omp_get_wtime();
    cout << "------------------------------------------------------\n";
    printf("LinR time is  ------  %f ms\n", (end_time - start_time) * 1000);
    cout << "------------------------------------------------------\n";

    evaluate(X_test, Y_test, w, b);

    printf("Info: successfully returned.\n");
    return 0;

}
