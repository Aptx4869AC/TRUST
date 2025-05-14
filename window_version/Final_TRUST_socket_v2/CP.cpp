#include "include/fastPai.h"
#include "include/cp_protocol.h"
#include "include/network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

int main()
{
#ifdef _WIN32
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0)
    {
        std::cerr << "WSAStartup failed: " << WSAGetLastError() << std::endl;
        return EXIT_FAILURE;
    }
#endif

    printf("CP run success\n");

    double start_time, end_time;

    double Time_Trust_FMUL = 0;
    double Time_Trust_FCMP = 0;
    double Time_Trust_FEQL = 0;
    double Time_Trust_FABS = 0;
    double Time_Trust_FTRN = 0;

    setrandom(&randstate);
    int server_sock, data_owner_sock;
    server_socket_initial(server_sock, data_owner_sock, 8081);

    // 接收数据
    mpz_t sk1, pk_N, pk_h, pk_h_N, pk_L, pk_L_k;
    mpz_inits(sk1, pk_N, pk_h, pk_h_N, pk_L, pk_L_k, NULL);
    recv_mpz(data_owner_sock, sk1);    // sk_1
    recv_mpz(data_owner_sock, pk_N);   // N
    recv_mpz(data_owner_sock, pk_h);   // h
    recv_mpz(data_owner_sock, pk_h_N); // h_N
    recv_mpz(data_owner_sock, pk_L);   // L
    recv_mpz(data_owner_sock, pk_L_k); // L_k
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

    printf("CP init is ok!\n");
    printf("----------------------------------------------------------\n");

    server_socket_initial(server_sock, data_owner_sock, 8084);

    // Trust 协议
    while (true)
    {
        mpz_t mpz_type, em_1, em_2, em_3, eres;
        mpz_inits(mpz_type, em_1, em_2, em_3, eres, NULL);
        printf("CP is waiting for input\n");
        recv_mpz(data_owner_sock, mpz_type); // 操作类型

        gmp_printf("------ mpz_type = %Zd\n", mpz_type);

        if (mpz_cmp_si(mpz_type, 1) == 0)    // MUL
        {
            recv_mpz(data_owner_sock, em_1);
            recv_mpz(data_owner_sock, em_2);
            gmp_printf("------ em_1 = %Zd\n", em_1);
            gmp_printf("------ em_2 = %Zd\n", em_2);

            start_time = omp_get_wtime();
            CP_TEE_SMUL(eres, em_1, em_2, cp);
            end_time = omp_get_wtime();
            Time_Trust_FMUL += (end_time - start_time);
            printf("Trust_FMUL time is  ------  %f ms\n", Time_Trust_FMUL * 1000);

            send_mpz(data_owner_sock, eres);
        }
        else if (mpz_cmp_si(mpz_type, 2) == 0) // CMP
        {
            recv_mpz(data_owner_sock, em_1);
            recv_mpz(data_owner_sock, em_2);
            gmp_printf("------ em_1 = %Zd\n", em_1);
            gmp_printf("------ em_2 = %Zd\n", em_2);

            start_time = omp_get_wtime();
            CP_TEE_SCMP_version2(eres, em_1, em_2, cp);
            end_time = omp_get_wtime();
            Time_Trust_FCMP += (end_time - start_time);
            printf("Trust_FCMP time is  ------  %f ms\n", Time_Trust_FCMP * 1000);

            send_mpz(data_owner_sock, eres);
        }
        else if (mpz_cmp_si(mpz_type, 3) == 0) // EQL
        {
            recv_mpz(data_owner_sock, em_1);
            recv_mpz(data_owner_sock, em_2);
            gmp_printf("------ em_1 = %Zd\n", em_1);
            gmp_printf("------ em_2 = %Zd\n", em_2);

            start_time = omp_get_wtime();
            CP_Trust_FEQL(eres, em_1, em_2, cp);
            end_time = omp_get_wtime();
            Time_Trust_FEQL += (end_time - start_time);
            printf("Trust_FEQL time is  ------  %f ms\n", Time_Trust_FEQL * 1000);

            send_mpz(data_owner_sock, eres);
        }
        else if (mpz_cmp_si(mpz_type, 4) == 0) // ABS
        {
            recv_mpz(data_owner_sock, em_1);
            recv_mpz(data_owner_sock, em_2);
            gmp_printf("------ em_1 = %Zd\n", em_1);
            gmp_printf("------ em_2 = %Zd\n", em_2);

            start_time = omp_get_wtime();
            CP_Trust_FABS(eres, em_1, cp);
            end_time = omp_get_wtime();
            Time_Trust_FABS += (end_time - start_time);
            printf("Trust_FABS time is  ------  %f ms\n", Time_Trust_FABS * 1000);

            send_mpz(data_owner_sock, eres);
        }
        else if (mpz_cmp_si(mpz_type, 5) == 0) // FTRN version1
        {
            recv_mpz(data_owner_sock, em_1);
            recv_mpz(data_owner_sock, em_2);
            recv_mpz(data_owner_sock, em_3);
            gmp_printf("------ em_1 = %Zd\n", em_1);
            gmp_printf("------ em_2 = %Zd\n", em_2);
            gmp_printf("------ em_3 = %Zd\n", em_3);

            start_time = omp_get_wtime();
            CP_Trust_FTRN_version1(eres, em_1, em_2, em_3, cp);
            end_time = omp_get_wtime();
            Time_Trust_FTRN += (end_time - start_time);
            printf("Trust_FTRN time is  ------  %f ms\n", Time_Trust_FTRN * 1000);

            send_mpz(data_owner_sock, eres);
        }

        mpz_clears(mpz_type, em_1, em_2, em_3, eres, NULL);
    }

#ifdef _WIN32
    WSACleanup();
#endif

    return 0;
}