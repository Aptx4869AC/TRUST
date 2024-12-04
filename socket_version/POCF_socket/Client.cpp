/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include "include/pcpd.h"
#include "include/network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

int main() {

    double start_time, end_time;
    double average_time_keygen = 0;
    double average_time_enc = 0;
    double average_time_dec = 0;
    double average_time_HE_Addition = 0;
    double average_time_HE_ScalarMul = 0;
    double average_time_HE_Subtraction = 0;
    double average_time_PDec_TDec = 0;

    int epoch = 10;
    int key_len = 1024; // 256, 512, 768, 1024, 1280, 1536
    printf("epoch = %d\n", epoch);
    printf("key_len = %d\n", key_len);
    printf("N 的比特长度为 %d bits\n", 2 * key_len);
    printf("----------------------------------------------------------\n");
    
    Paillier pai;
    pai.KeyGen(key_len, sigma, 0); // k, sigma, yita

    start_time = omp_get_wtime();
    Paillier temp_pai;
    for (int i = 0; i < epoch; i++) {
        temp_pai.KeyGen(key_len, sigma, 0); // k, sigma, yita
    }
    end_time = omp_get_wtime();
    average_time_keygen += (end_time - start_time);

    Paillier_Third_Party cp(pai.public_key, pai.private_key_1);
    Paillier_Third_Party csp(pai.public_key, pai.private_key_2);
    for (int i = 0; i < epoch; i++) {
        mpz_t a, ea;
        mpz_inits(a, ea, NULL);
        mpz_rrandomb(a, randstate, 16);
//        gmp_printf("a = %Zd\n", a);

        start_time = omp_get_wtime();
        Enc(ea, pai.public_key, a);
        end_time = omp_get_wtime();
        average_time_enc += (end_time - start_time);
//        printf("%d [Enc] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        start_time = omp_get_wtime();
        Dec(a, pai.private_key, ea);
        end_time = omp_get_wtime();
        average_time_dec += (end_time - start_time);
//        printf("%d [Dec] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t a1, a2;
        mpz_t a_cp, a_csp;
        mpz_inits(a1, a2, a_cp, a_csp, NULL);
        start_time = omp_get_wtime();
        PDec(a1, cp.private_key, ea);
        PDec(a2, csp.private_key, ea);
        TDec(a_cp, a1, a2, cp.public_key.N_square, cp.public_key.N, cp.public_key.N_half);
        end_time = omp_get_wtime();
        average_time_PDec_TDec += (end_time - start_time);
//        printf("%d [PDec and TDec] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t cz;
        mpz_init(cz);
        start_time = omp_get_wtime();
        mpz_mul(cz, ea, ea);
        mpz_mod(cz, cz, pai.public_key.N_square);
        end_time = omp_get_wtime();
        average_time_HE_Addition += (end_time - start_time);
//        printf("%d [HE_Addition] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t ten;
        mpz_init(ten);
        mpz_set_si(ten, 10);
        start_time = omp_get_wtime();
        mpz_powm(cz, ea, ten, pai.public_key.N_square);
        end_time = omp_get_wtime();
        average_time_HE_ScalarMul += (end_time - start_time);
//        printf("%d [HE_ScalarMul] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t neg_one;
        mpz_init(neg_one);
        mpz_set_si(neg_one, -1);
        start_time = omp_get_wtime();
        mpz_powm(cz, ea, neg_one, pai.public_key.N_square);
        mpz_mul(cz, ea, ea);
        mpz_mod(cz, cz, pai.public_key.N_square);
        end_time = omp_get_wtime();
        average_time_HE_Subtraction += (end_time - start_time);
//        printf("%d [HE_Subtraction] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_clears(a, ea, a1, a2, a_cp, a_csp, NULL);
//        printf("----------------------------------------------------------\n");
    }


    printf("\n\n");
    printf("KeyGen (%d average) time is  ------  %f ms\n", epoch, average_time_keygen / epoch * 1000);
    printf("Enc (%d average) time is  ------  %f ms\n", epoch, average_time_enc / epoch * 1000);
    printf("Dec (%d average) time is  ------  %f ms\n", epoch, average_time_dec / epoch * 1000);
    printf("PDec + TDec (%d average) time is  ------  %f ms\n", epoch, average_time_PDec_TDec / epoch * 1000);
    printf("HE_Addition (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Addition / epoch * 1000);
    printf("HE_ScalarMul (%d average) time is  ------  %f ms\n", epoch, average_time_HE_ScalarMul / epoch * 1000);
    printf("HE_Subtraction (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Subtraction / epoch * 1000);
    printf("----------------------------------------------------------\n");

    mpz_t x, y, ex, ey;
    mpz_inits(x, y, ex, ey, NULL);
    setrandom(&randstate);
    mpz_rrandomb(x, randstate, 8); // 8-bit random number
    mpz_rrandomb(y, randstate, 8); // 8-bit random number
    Enc(ex, pai.public_key, x);
    Enc(ey, pai.public_key, y);

    // 初始化CP
    int sock_to_cp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8081);

    // 发给CP <sk1, pai.prikey.n, ex, ey>
    gmp_printf("Send (PK,SK) to CP, sk1 = %Zd\n", pai.private_key_1.lambda);

    send_mpz(sock_to_cp, pai.private_key_1.lambda);
    send_mpz(sock_to_cp, pai.private_key_1.lambda_reverse);
    send_mpz(sock_to_cp, pai.private_key_1.N);
    send_mpz(sock_to_cp, ex);
    send_mpz(sock_to_cp, ey);
    send_mpz(sock_to_cp, x); //用于验证
    send_mpz(sock_to_cp, y); //用于验证
    send_mpz(sock_to_cp, pai.private_key.lambda); //用于验证
    send_mpz(sock_to_cp, pai.private_key.lambda_reverse); //用于验证

    // 初始化CSP
    int sock_to_csp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8082);

    // 发给CSP <sk2, pai.prikey.n, ex, ey>
    gmp_printf("Send (PK,SK) to CSP, sk2 = %Zd\n", pai.private_key_2.lambda);

    send_mpz(sock_to_csp, pai.private_key_2.lambda);
    send_mpz(sock_to_csp, pai.private_key_2.lambda_reverse);
    send_mpz(sock_to_csp, pai.private_key_2.N);

    gmp_printf("------ n = %Zd\n", pai.public_key.N);
    gmp_printf("------ x = %Zd\n", x);
    gmp_printf("------ y = %Zd\n", y);
    gmp_printf("------ ex = %Zd\n", ex);
    gmp_printf("------ ey = %Zd\n", ey);

    // 收到结束标志
    close(sock_to_cp); // printf("关闭CP的连接\n");
    close(sock_to_csp); // printf("关闭CSP的连接\n");
    return 0;
}