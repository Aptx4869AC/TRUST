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
    double start_time, end_time;

    double average_time_enc = 0;
    double average_time_dec = 0;
    double average_time_PDec_TDec = 0;
    double average_time_HE_Addition = 0;
    double average_time_HE_ScalarMul = 0;
    double average_time_HE_Subtraction = 0;

    double average_time_TEE_SMUL = 0;
    double average_time_TEE_SCMP_version1 = 0;
    double average_time_TEE_SCMP_version2 = 0;
    double average_time_TEE_SABS = 0;
    double average_time_TEE_SDIV = 0;
    double average_time_SMUL = 0;
    double average_time_SCMP = 0;
    double average_time_SSBA = 0;
    double average_time_SDIV = 0;
    double average_time_Trust_FEQL = 0;
    double average_time_Trust_FABS = 0;
    double average_time_Trust_FTRN_v1 = 0;
    double average_time_Trust_FTRN_v2 = 0;

    double average_bytes_TEE_SMUL = 0;
    double average_bytes_TEE_SCMP_version1 = 0;
    double average_bytes_TEE_SCMP_version2 = 0;
    double average_bytes_TEE_SABS = 0;
    double average_bytes_TEE_SDIV = 0;
    double average_bytes_SMUL = 0;
    double average_bytes_SCMP = 0;
    double average_bytes_SSBA = 0;
    double average_bytes_SDIV = 0;
    double average_bytes_Trust_FEQL = 0;
    double average_bytes_Trust_FABS = 0;
    double average_bytes_Trust_FTRN_v1 = 0;
    double average_bytes_Trust_FTRN_v2 = 0;
    double average_bytes_PK = 0;
    double average_bytes_SK = 0;

    setrandom(&randstate);
    int server_sock, data_owner_sock;
    server_socket_initial(server_sock, data_owner_sock, 8081);

    // 接收数据
    mpz_t sk1, pk_N, pk_h, pk_h_N, pk_L, pk_L_k, ex, ey, x, y, pai_alpha;
    mpz_inits(sk1, pk_N, pk_h, pk_h_N, pk_L, pk_L_k, ex, ey, x, y, pai_alpha, NULL);
    recv_mpz(data_owner_sock, sk1);
    recv_mpz(data_owner_sock, pk_N); // N
    recv_mpz(data_owner_sock, pk_h); // h
    recv_mpz(data_owner_sock, pk_h_N); // h_N
    recv_mpz(data_owner_sock, pk_L); // L
    recv_mpz(data_owner_sock, pk_L_k); //L_k
    recv_mpz(data_owner_sock, ex);
    recv_mpz(data_owner_sock, ey);
    recv_mpz(data_owner_sock, x); //用于验证
    recv_mpz(data_owner_sock, y); //用于验证
    recv_mpz(data_owner_sock, pai_alpha); //用于验证
    server_socket_close(data_owner_sock, server_sock);

    gmp_printf("------ sk1 = %Zd\n", sk1);
    gmp_printf("------ N = %Zd\n", pk_N);
    gmp_printf("------ h = %Zd\n", pk_h);
    gmp_printf("------ h_N = %Zd\n", pk_h_N);
    gmp_printf("------ L = %Zd\n", pk_L);
    gmp_printf("------ L_k = %Zd\n", pk_L_k);
    gmp_printf("------ ex = %Zd\n", ex);
    gmp_printf("------ ey = %Zd\n", ey);



    // 配置自身
    PublicKey pk(pk_N, pk_h, pk_h_N, pk_L, mpz_get_si(pk_L_k));
    PaillierThdDec cp = PaillierThdDec(pk_N, sk1, pk, sigma);

    PrivateKey sk(pai_alpha, pk_N);
    Paillier pai(pk, sk, sigma); // 用于验证
    printf("----------------------------------------------------------\n");

    int epoch = 10;
    bool record_bytes_flag = false; // 控制存储开销的记录按钮

    // 原始SOCI+ 协议
    for (int i = 0; i < epoch; i++)
    {
        // 临时更换新的x，y，ex,ey
        mpz_rrandomb(x, randstate, 8); // 8-bit random number
        mpz_rrandomb(y, randstate, 8); // 8-bit random number
        cp.pai.encrypt(ex, x);
        cp.pai.encrypt(ey, y);

        mpz_t eres, es_x, eu_x, eq, ee, q, r, right_q, right_r;
        mpz_inits(eres, es_x, eu_x, eq, ee, q, r, right_q, right_r, NULL);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_SMUL(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SMUL += end_time - start_time;
        check_smul(i, eres, x, y, pai); // 验证正确性
        average_bytes_SMUL += total_bytes / 1024.0;
//        printf("%d [SMUL] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_SCMP(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SCMP += end_time - start_time;
        check_scmp(i, eres, x, y, pai); // 验证正确性
        average_bytes_SCMP += total_bytes / 1024.0;
//        printf("%d [SCMP] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        SSBA(es_x, eu_x, ex, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SSBA += end_time - start_time;
        check_ssba(i, es_x, eu_x, x, pai); // 验证正确性
        average_bytes_SSBA += total_bytes / 1024.0;
//        printf("%d [SSBA] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        SDIV(eq, ee, ex, ey, 10, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SDIV += (end_time - start_time);
        check_sdiv(i, eq, ee, x, y, pai); // 验证正确性
        average_bytes_SDIV += total_bytes / 1024.0;
//        printf("%d [SDIV] time is  ------  %f ms\n", i + 1, (end_time - start_time) * 1000);

        mpz_clears(eres, es_x, eu_x, eq, ee, q, r, right_q, right_r, NULL);
//        printf("----------------------------------------------------------\n");

    }

    // SOCI-TEE 协议
    for (int i = 0; i < epoch; i++)
    {
        // 临时更换新的x，y，ex,ey
        mpz_rrandomb(x, randstate, 8); // 8-bit random number
        mpz_rrandomb(y, randstate, 8); // 8-bit random number
        cp.pai.encrypt(ex, x);
        cp.pai.encrypt(ey, y);

        mpz_t eres, es_x, eu_x, eq, ee, q, r, right_q, right_r;
        mpz_inits(eres, es_x, eu_x, eq, ee, q, r, right_q, right_r, NULL);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_TEE_SMUL(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_TEE_SMUL += end_time - start_time;
        check_smul(i, eres, x, y, pai); // 验证正确性
        average_bytes_TEE_SMUL += total_bytes / 1024.0;
//        printf("%d [NEW-SMUL] time is %f ms\n", i + 1, (end_time - start_time) * 1000);


        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_TEE_SCMP_version1(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_TEE_SCMP_version1 += end_time - start_time;
        check_scmp(i, eres, x, y, pai); // 验证正确性
        average_bytes_TEE_SCMP_version1 += total_bytes / 1024.0;
//        printf("%d [NEW-SCMP_version1] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_TEE_SCMP_version2(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_TEE_SCMP_version2 += end_time - start_time;
        check_scmp(i, eres, x, y, pai); // 验证正确性
        average_bytes_TEE_SCMP_version2 += total_bytes / 1024.0;
//        printf("%d [NEW-SCMP_version2] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_TEE_SABS(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_TEE_SABS += end_time - start_time;
        check_sabs(i, eres, x, y, pai); // 验证正确性
        average_bytes_TEE_SABS += total_bytes / 1024.0;
//        printf("%d [NEW-SABS] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        TEE_SDIV(eq, ee, ex, ey, 10, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_TEE_SDIV += (end_time - start_time);
        check_sdiv(i, eq, ee, x, y, pai); // 验证正确性
        average_bytes_TEE_SDIV += total_bytes / 1024.0;
//        printf("%d [NEW-SDIV] time is  ------  %f ms\n", i + 1, (end_time - start_time) * 1000);

        mpz_clears(eres, es_x, eu_x, eq, ee, q, r, right_q, right_r, NULL);
//        printf("----------------------------------------------------------\n");
    }

    ifstream file("Data/training_dataADD.csv");  // 打开CSV文件
    string line, value;

    if (!file.is_open())
    {
        cerr << "无法打开文件 training_dataADD.csv" << '\n';
        return 1;
    }
    // 跳过表头
    getline(file, line);

    int _epoch = 1;

    // Trust 协议
    while (getline(file, line))   // 逐行读取CSV文件的内容
    {
        stringstream ss(line);
        getline(ss, value, ',');
        int value1 = stoi(value);  // 转换为整数
        getline(ss, value, ',');
        int value2 = stoi(value);  // 转换为整数

        mpz_t m_1, m_2, m_3, em_1, em_2, em_3, eres, res;
        mpz_inits(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
        mpz_set_si(m_1, value1);
        mpz_set_si(m_2, value2);
        if (mpz_cmp_si(m_1, 0) < 0)
        {
            mpz_add(m_1, m_1, cp.pai.pk.N);
        }
        if (mpz_cmp_si(m_2, 0) < 0)
        {
            mpz_add(m_2, m_2, cp.pai.pk.N);
        }
        mpz_rrandomb(m_3, randstate, 8);
        cp.pai.encrypt(em_1, m_1);
        cp.pai.encrypt(em_2, m_2);
        cp.pai.encrypt(em_3, m_3);

        // EQL
        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_Trust_FEQL(eres, em_1, em_2, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FEQL += (end_time - start_time);
        check_feql(_epoch, eres, m_1, m_2, pai); // 验证正确性
        average_bytes_Trust_FEQL += total_bytes / 1024.0;
//        printf("%d [Trust_FEQL] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        // ABS
        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_Trust_FABS(eres, em_1, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FABS += (end_time - start_time);
        check_fabs(_epoch, eres, m_1, pai); // 验证正确性
        average_bytes_Trust_FABS += total_bytes / 1024.0;
//        printf("%d [Trust_FABS] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        mpz_set_si(m_1, 1);
        pai.encrypt(em_1, m_1);

        // FTRN version1
        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_Trust_FTRN_version1(eres, em_1, em_2, em_3, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FTRN_v1 += (end_time - start_time);
        check_ftrn(_epoch, eres, m_1, m_2, m_3, pai); // 验证正确性
        average_bytes_Trust_FTRN_v1 += total_bytes / 1024.0;
//        printf("%d [Trust_FTRN_v1] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        int pi = (int) random() % 2;
        mpz_set_si(m_1, pi);
        pai.encrypt(em_1, m_1);
        // FTRN version2
        total_bytes = 0.0;
        start_time = omp_get_wtime();
        CP_Trust_FTRN_version2(eres, em_1, em_2, em_3, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FTRN_v2 += (end_time - start_time);
        check_ftrn(_epoch, eres, m_1, m_2, m_3, pai); // 验证正确性
        average_bytes_Trust_FTRN_v2 += total_bytes / 1024.0;
//        printf("%d [Trust_FTRN_v2] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        mpz_clears(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
//        printf("----------------------------------------------------------\n");


        _epoch++;
        if (_epoch > 10)
        {
            break;
        }
    }
    
    // 协议测试
    printf("SMUL (%d average) time is  ------  %f ms\n", epoch, average_time_SMUL / epoch * 1000);
    printf("SCMP (%d average) time is  ------  %f ms\n", epoch, average_time_SCMP / epoch * 1000);
    printf("SSBA (%d average) time is  ------  %f ms\n", epoch, average_time_SSBA / epoch * 1000);
    printf("SDIV (10-bits) (%d average) time is  ------  %f ms\n", epoch, average_time_SDIV / epoch * 1000);
    printf("NEW-SMUL (%d average) time is  ------  %f ms\n", epoch, average_time_TEE_SMUL / epoch * 1000);
    printf("NEW-SCMP_version1 (%d average) time is  ------  %f ms\n", epoch, average_time_TEE_SCMP_version1 / epoch * 1000);
    printf("NEW-SCMP_version2 (%d average) time is  ------  %f ms\n", epoch, average_time_TEE_SCMP_version2 / epoch * 1000);
    printf("NEW-SABS (%d average) time is  ------  %f ms\n", epoch, average_time_TEE_SABS / epoch * 1000);
    printf("NEW-SDIV (%d average) time is  ------  %f ms\n", epoch, average_time_TEE_SDIV / epoch * 1000);
    printf("Trust_FEQL (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FEQL / _epoch * 1000);
    printf("Trust_FABS (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FABS / _epoch * 1000);
    printf("Trust_FTRN_v1 (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FTRN_v1 / _epoch * 1000);
    printf("Trust_FTRN_v2 (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FTRN_v2 / _epoch * 1000);

    if (record_bytes_flag)
    {
        printf("SMUL (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SMUL / epoch);
        printf("SCMP (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SCMP / epoch);
        printf("SSBA (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SSBA / epoch);
        printf("SDIV (10-bits) (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SDIV / epoch);
        printf("NEW-SMUL (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_TEE_SMUL / epoch);
        printf("NEW-SCMP_version1 (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_TEE_SCMP_version1 / epoch);
        printf("NEW-SCMP_version2 (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_TEE_SCMP_version2 / epoch);
        printf("NEW-SABS (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_TEE_SABS / epoch);
        printf("NEW-SDIV (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_TEE_SDIV / epoch);
        printf("Trust_FEQL (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FEQL / _epoch);
        printf("Trust_FABS (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FABS / _epoch);
        printf("Trust_FTRN_v1 (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FTRN_v1 / _epoch);
        printf("Trust_FTRN_v2 (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FTRN_v2 / _epoch);
    }

    return 0;
}