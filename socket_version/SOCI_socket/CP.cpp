/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include "include/paillier.h"
#include "include/cp_protocol.h"
#include "include/network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;


int main()
{
    double start_time, end_time;
    double average_time_SMUL = 0;
    double average_time_SCMP = 0;
    double average_time_SSBA = 0;
    double average_time_SDIV = 0;
    double average_time_Trust_FEQL = 0;
    double average_time_Trust_FABS = 0;
    double average_time_Trust_FTRN_v1 = 0;
    double average_time_Trust_FTRN_v2 = 0;


    double average_bytes_SMUL = 0;
    double average_bytes_SCMP = 0;
    double average_bytes_SSBA = 0;
    double average_bytes_SDIV = 0;
    double average_bytes_Trust_FEQL = 0;
    double average_bytes_Trust_FABS = 0;
    double average_bytes_Trust_FTRN_v1 = 0;
    double average_bytes_Trust_FTRN_v2 = 0;

    setrandom(&randstate);
    int server_sock, data_owner_sock;
    server_socket_initial(server_sock, data_owner_sock, 8081);


    // 接收数据
    mpz_t sk1, n, nsquare, ex, ey;
    mpz_inits(sk1, n, nsquare, ex, ey, NULL);
    recv_mpz(data_owner_sock, sk1);
    recv_mpz(data_owner_sock, n);
    recv_mpz(data_owner_sock, ex);
    recv_mpz(data_owner_sock, ey);

    mpz_t x, y, pai_lambda;
    mpz_inits(x, y, pai_lambda, NULL);
    recv_mpz(data_owner_sock, x); //用于验证
    recv_mpz(data_owner_sock, y); //用于验证
    recv_mpz(data_owner_sock, pai_lambda); //用于验证
    server_socket_close(data_owner_sock, server_sock);

    mpz_mul(nsquare, n, n);
    gmp_printf("------ sk1 = %Zd\n", sk1);
    gmp_printf("------ n = %Zd\n", n);
    gmp_printf("------ nsquare = %Zd\n", nsquare);

    // 配置自身
    PaillierKey pubkey(n);
    PaillierThdPrivateKey share_part1(sk1, n, nsquare);
    PaillierThdDec cp(share_part1, pubkey);

    PaillierPrivateKey prikey(n, pai_lambda); //用于验证
    Paillier pai(pubkey, prikey); //用于验证

    printf("----------------------------------------------------------\n");
    mpz_set_si(neg_one, -1);
    mpz_set_si(one, 1);
    mpz_set_si(zero, 0);

    int epoch = 10;
    bool record_bytes_flag = true;
    for (int i = 0; i < epoch; i++)
    {
        //临时更换新的x，y，ex,ey
        mpz_rrandomb(x, randstate, 8); // 8-bit random number
        mpz_rrandomb(y, randstate, 8); // 8-bit random number
        cp.pai.encrypt(ex, x);
        cp.pai.encrypt(ey, y);

        mpz_t eres, es_x, eu_x, eq, ee, q, r, right_q, right_r;
        mpz_inits(eres, es_x, eu_x, eq, ee, q, r, right_q, right_r, NULL);

        total_bytes = 0;
        start_time = omp_get_wtime();
        CP_SMUL(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SMUL += end_time - start_time;
        check_smul(i, eres, x, y, pai); // 验证正确性
//        printf("%d [SMUL] time is %f ms\n", i + 1, (end_time - start_time) * 1000);
//        printf("total_bytes = %f KB\n", total_bytes / 1024.0);
        average_bytes_SMUL += total_bytes / 1024.0;

        total_bytes = 0;
        start_time = omp_get_wtime();
        CP_SCMP(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SCMP += end_time - start_time;
        check_scmp(i, eres, x, y, pai); // 验证正确性
//        printf("%d [SCMP] time is %f ms\n", i + 1, (end_time - start_time) * 1000);
//        printf("total_bytes = %f KB\n", total_bytes / 1024.0);
        average_bytes_SCMP += total_bytes / 1024.0;

        total_bytes = 0;
        start_time = omp_get_wtime();
        SSBA(es_x, eu_x, ex, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SSBA += end_time - start_time;
        check_ssba(i, es_x, eu_x, x, pai); // 验证正确性
//        printf("%d [SSBA] time is %f ms\n", i + 1, (end_time - start_time) * 1000);
//        printf("total_bytes = %f KB\n", total_bytes / 1024.0);
        average_bytes_SSBA += total_bytes / 1024.0;

        total_bytes = 0;
        start_time = omp_get_wtime();
        SDIV(eq, ee, ex, ey, 10, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SDIV += (end_time - start_time);
        check_sdiv(i, eq, ee, x, y, pai); // 验证正确性
//        printf("%d [SDIV] time is  ------  %f ms\n", i + 1, (end_time - start_time) * 1000);
//        printf("total_bytes = %f KB\n", total_bytes / 1024.0);
        average_bytes_SDIV += total_bytes / 1024.0;
//        printf("----------------------------------------------------------\n");

        mpz_clears(eres, es_x, eu_x, eq, ee, q, r, right_q, right_r, NULL);

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
    int _epoch = 0;
    // Trust 新协议
    while (getline(file, line))   // 逐行读取CSV文件的内容
    {
        if (_epoch > 40)
        {
            break;
        }


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
            mpz_add(m_1, m_1,  cp.pai.pubkey.n);
        }
        if (mpz_cmp_si(m_2, 0) < 0)
        {
            mpz_add(m_2, m_2,  cp.pai.pubkey.n);
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
    }

    printf("SMUL (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SMUL / epoch * 1000);
    printf("SCMP (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SCMP / epoch * 1000);
    printf("SSBA (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SSBA / epoch * 1000);
    printf("SDIV (10-bits) (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SDIV / epoch * 1000);

    /** 以下TRUST供参考用 **/
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

        /** 以下TRUST供参考用 **/
        printf("Trust_FEQL (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FEQL / _epoch);
        printf("Trust_FABS (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FABS / _epoch);
        printf("Trust_FTRN_v1 (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FTRN_v1 / _epoch);
        printf("Trust_FTRN_v2 (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FTRN_v2 / _epoch);
    }


    return 0;
}