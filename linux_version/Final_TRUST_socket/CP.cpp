/******************************************************************************
 * Author: Aptx4869AC 
 * Created: 2024-12-01 
 * GitHub: https://github.com/Aptx4869AC  
 *****************************************************************************/

#include "include/fastPai.h"
#include "include/cp_protocol.h"
#include "include/network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

int main()
{

    printf("CP run success\n");

    double start_time, end_time;

    double average_time_Trust_FMUL = 0;
    double average_time_Trust_FCMP = 0;
    double average_time_Trust_FEQL = 0;
    double average_time_Trust_FABS = 0;
    double average_time_Trust_FTRN = 0;

    setrandom(&randstate);
    int server_sock, data_owner_sock;
    server_socket_initial(server_sock, data_owner_sock, 8081);

    // 接收数据
    mpz_t sk1, pk_N, pk_h, pk_h_N, pk_L, pk_L_k, pai_alpha;
    mpz_inits(sk1, pk_N, pk_h, pk_h_N, pk_L, pk_L_k, pai_alpha, NULL);
    recv_mpz(data_owner_sock, sk1);
    recv_mpz(data_owner_sock, pk_N);      // N
    recv_mpz(data_owner_sock, pk_h);      // h
    recv_mpz(data_owner_sock, pk_h_N);    // h_N
    recv_mpz(data_owner_sock, pk_L);      // L
    recv_mpz(data_owner_sock, pk_L_k);    // L_k
    recv_mpz(data_owner_sock, pai_alpha); // 用于验证
    server_socket_close(data_owner_sock, server_sock);

    gmp_printf("------ sk1 = %Zd\n", sk1);
    gmp_printf("------ N = %Zd\n", pk_N);
    gmp_printf("------ h = %Zd\n", pk_h);
    gmp_printf("------ h_N = %Zd\n", pk_h_N);
    gmp_printf("------ L = %Zd\n", pk_L);
    gmp_printf("------ L_k = %Zd\n", pk_L_k);

    // 配置自身
    PublicKey pk(pk_N, pk_h, pk_h_N, pk_L, mpz_get_si(pk_L_k));
    PaillierThdDec cp = PaillierThdDec(pk_N, sk1, pk, sigma);

    PrivateKey sk(pai_alpha, pk_N);
    Paillier pai(pk, sk, sigma); // 用于验证
    printf("----------------------------------------------------------\n");

    ifstream file("../Data/training_dataADD.csv"); // 打开CSV文件
    string line, value;

    if (!file.is_open())
    {
        cerr << "Failed to open the file training_dataADD.csv" << '\n';
        return 1;
    }
    // 跳过表头
    getline(file, line);

    int _epoch = 0;

    // Trust 协议
    while (getline(file, line)) // 逐行读取CSV文件的内容
    {
        stringstream ss(line);
        getline(ss, value, ',');
        int value1 = stoi(value); // 转换为整数
        getline(ss, value, ',');
        int value2 = stoi(value); // 转换为整数

        mpz_t m_1, m_2, m_3, em_1, em_2, em_3, eres, res;
        mpz_inits(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
        mpz_set_si(m_1, value1);
        mpz_set_si(m_2, value2);
        mpz_rrandomb(m_3, randstate, 8);
        gmp_printf("[epoch = %d] --- m_1 = %Zd, m_2 = %Zd, m_3 = %Zd\n", _epoch + 1, m_1, m_2, m_3);
        if (mpz_cmp_si(m_1, 0) < 0)
        {
            mpz_add(m_1, m_1, cp.pai.pk.N);
        }
        if (mpz_cmp_si(m_2, 0) < 0)
        {
            mpz_add(m_2, m_2, cp.pai.pk.N);
        }
        cp.pai.encrypt(em_1, m_1);
        cp.pai.encrypt(em_2, m_2);
        cp.pai.encrypt(em_3, m_3);

        // MUL

        start_time = omp_get_wtime();
        CP_TEE_SMUL(eres, em_1, em_2, cp);
        end_time = omp_get_wtime();
        average_time_Trust_FMUL += (end_time - start_time);
        check_smul(_epoch, eres, m_1, m_2, pai); // 验证正确性

        // CMP
        start_time = omp_get_wtime();
        CP_TEE_SCMP_version2(eres, em_1, em_2, cp);
        end_time = omp_get_wtime();
        average_time_Trust_FCMP += (end_time - start_time);
        check_scmp(_epoch, eres, m_1, m_2, pai); // 验证正确性

        // EQL
        start_time = omp_get_wtime();
        CP_Trust_FEQL(eres, em_1, em_2, cp);
        end_time = omp_get_wtime();
        average_time_Trust_FEQL += (end_time - start_time);
        check_feql(_epoch, eres, m_1, m_2, pai); // 验证正确性

        // ABS
        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_Trust_FABS(eres, em_1, cp);
        end_time = omp_get_wtime();
        average_time_Trust_FABS += (end_time - start_time);
        check_fabs(_epoch, eres, m_1, pai); // 验证正确性

        mpz_set_si(m_1, 1);
        pai.encrypt(em_1, m_1);

        // FTRN version1
        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_Trust_FTRN_version1(eres, em_1, em_2, em_3, cp);
        end_time = omp_get_wtime();
        average_time_Trust_FTRN += (end_time - start_time);
        check_ftrn(_epoch, eres, m_1, m_2, m_3, pai); // 验证正确性

        mpz_clears(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);

        _epoch++;
        if (_epoch >= 100) // 只挑100行样本
        {
            break;
        }

        
    }

    printf("Trust_FMUL (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FMUL / _epoch * 1000);
    printf("Trust_FCMP (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FCMP / _epoch * 1000);
    printf("Trust_FEQL (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FEQL / _epoch * 1000);
    printf("Trust_FABS (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FABS / _epoch * 1000);
    printf("Trust_FTRN (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FTRN / _epoch * 1000);
    
    return 0;
}