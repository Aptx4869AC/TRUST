/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include "include/pcpd.h"
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
    double average_time_SXOR = 0;
    double average_time_SEQ = 0;

    double average_bytes_SMUL = 0;
    double average_bytes_SCMP = 0;
    double average_bytes_SSBA = 0;
    double average_bytes_SDIV = 0;
    double average_bytes_SXOR = 0;
    double average_bytes_SEQ = 0;

    setrandom(&randstate);
    int server_sock, data_owner_sock;
    server_socket_initial(server_sock, data_owner_sock, 8081);

    mpz_inits(zero, one, two, neg_one, N, N_square, N_half, NULL);
    mpz_set_si(zero, 0);
    mpz_set_si(one, 1);
    mpz_set_si(two, 2);
    mpz_set_si(neg_one, -1);

    // 接收数据
    mpz_t lambda, lambda_reverse, ex, ey;
    mpz_inits(lambda, lambda_reverse, ex, ey, NULL);
    recv_mpz(data_owner_sock, lambda);
    recv_mpz(data_owner_sock, lambda_reverse);
    recv_mpz(data_owner_sock, N);
    recv_mpz(data_owner_sock, ex);
    recv_mpz(data_owner_sock, ey);

    mpz_t x, y, pai_lambda, pai_lambda_reverse;
    mpz_inits(x, y, pai_lambda, pai_lambda_reverse, NULL);
    recv_mpz(data_owner_sock, x); //用于验证
    recv_mpz(data_owner_sock, y); //用于验证
    recv_mpz(data_owner_sock, pai_lambda); //用于验证
    recv_mpz(data_owner_sock, pai_lambda_reverse); //用于验证
    server_socket_close(data_owner_sock, server_sock);

    gmp_printf("------ sk1 = %Zd\n", lambda);
    gmp_printf("------ n = %Zd\n", N);

    PublicKey public_key(N, sigma);
    PrivateKey private_key_1(lambda, lambda_reverse, N);
    Paillier_Third_Party cp(public_key, private_key_1);

    PrivateKey pai_private_key(pai_lambda, pai_lambda_reverse, N);
    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);
    mpz_set(N_half, cp.public_key.N_half);
    printf("----------------------------------------------------------\n");

    int epoch = 10;
    bool record_bytes_flag = true;

    // POCF原协议
    for (int i = 0; i < epoch; i++)
    {
        //临时更换新的x，y，ex,ey
        mpz_rrandomb(x, randstate, 8); // 8-bit random number
        mpz_rrandomb(y, randstate, 8); // 8-bit random number
        Enc(ex, cp.public_key, x);
        Enc(ey, cp.public_key, y);

        mpz_t eres, es_x, eu_x, eq, ee, q, r, right_q, right_r;
        mpz_inits(eres, es_x, eu_x, eq, ee, q, r, right_q, right_r, NULL);

        total_bytes = 0;
        start_time = omp_get_wtime();
        CP_SMUL(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SMUL += end_time - start_time;
        check_smul(i, eres, x, y, pai_private_key); // 验证正确性
        average_bytes_SMUL += total_bytes / 1024.0;
//        printf("%d [SMUL] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0;
        start_time = omp_get_wtime();
        CP_SCMP(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SCMP += end_time - start_time;
        check_scmp(i, eres, x, y, pai_private_key); // 验证正确性
        average_bytes_SCMP += total_bytes / 1024.0;
//        printf("%d [SCMP] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0;
        start_time = omp_get_wtime();
        CP_SSBA(es_x, eu_x, ex, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SSBA += end_time - start_time;
        check_ssba(i, es_x, eu_x, x, pai_private_key); // 验证正确性
        average_bytes_SSBA += total_bytes / 1024.0;
//        printf("%d [SSBA] time is %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0;
        start_time = omp_get_wtime();
        SDIV(eq, ee, ex, ey, 10, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SDIV += (end_time - start_time);
        check_sdiv(i, eq, ee, x, y, pai_private_key); // 验证正确性
        average_bytes_SDIV += total_bytes / 1024.0;
//        printf("%d [SDIV] time is  ------  %f ms\n", i + 1, (end_time - start_time) * 1000);

        total_bytes = 0;
        start_time = omp_get_wtime();
        SXOR(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_SXOR += (end_time - start_time);
        average_bytes_SXOR += total_bytes / 1024.0;


        total_bytes = 0;
        start_time = omp_get_wtime();
        SEQ(eres, ex, ey, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        check_seq(i, eres, x, y, pai_private_key);
        average_time_SEQ += (end_time - start_time);
        average_bytes_SEQ += total_bytes / 1024.0;

        mpz_clears(eres, es_x, eu_x, eq, ee, q, r, right_q, right_r, NULL);
//        printf("----------------------------------------------------------\n");
    }

    printf("SMUL (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SMUL / epoch * 1000);
    printf("SCMP (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SCMP / epoch * 1000);
    printf("SSBA (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SSBA / epoch * 1000);
    printf("SDIV (10-bits) (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SDIV / epoch * 1000);
    printf("SXOR (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SXOR / epoch * 1000);
    printf("SEQ  (%d average) Computation Costs is  ---  %f ms\n", epoch, average_time_SEQ / epoch * 1000);
    if (record_bytes_flag)
    {
        printf("SMUL (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SMUL / epoch);
        printf("SCMP (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SCMP / epoch);
        printf("SSBA (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SSBA / epoch);
        printf("SDIV (10-bits) (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SDIV / epoch);
        printf("SXOR (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SXOR / epoch);
        printf("SEQ  (%d average) Communication Costs is  ---  %f KB\n", epoch, average_bytes_SEQ / epoch);
    }
}