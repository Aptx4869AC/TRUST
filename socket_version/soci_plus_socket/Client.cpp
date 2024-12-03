/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-01
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include "include/fastPai.h"
#include "include/network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

int main() {
    double start_time, end_time;
    double average_time_keygen = 0;
    double average_time_enc = 0;
    double average_time_dec = 0;
    double average_time_PDec_TDec = 0;
    double average_time_HE_Addition = 0;
    double average_time_HE_ScalarMul = 0;
    double average_time_HE_Subtraction = 0;

    int k = 112; //BIT SECURITY, 64、80、104、112、124、128
    KeyGen keyGen(k, sigma);
    Paillier pai(keyGen.pk, keyGen.sk, sigma);

    int epoch = 10;

    start_time = omp_get_wtime();
    for (int i = 0; i < epoch; i++) {
        KeyGen temp_keyGen(k, sigma);
    }
    end_time = omp_get_wtime();
    average_time_keygen += (end_time - start_time);

    // 临时副本，用来测试PDec与TDec
    PaillierThirdParty *psk = keyGen.pai_third_parties;
    PaillierThdDec cp = PaillierThdDec(psk[0].N, psk[0].partial_key, keyGen.pk, sigma);
    PaillierThdDec csp = PaillierThdDec(psk[1].N, psk[1].partial_key, keyGen.pk, sigma);

    for (int i = 0; i < epoch; i++) {
        mpz_t a, ea;
        mpz_inits(a, ea, NULL);
        mpz_rrandomb(a, randstate, 16);
//        gmp_printf("a = %Zd\n", a);

        start_time = omp_get_wtime();
        pai.encrypt(ea, a);
        end_time = omp_get_wtime();
        average_time_enc += (end_time - start_time);
//        printf("%d [Enc] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        start_time = omp_get_wtime();
        pai.decrypt(a, ea);
        end_time = omp_get_wtime();
        average_time_dec += (end_time - start_time);
//        printf("%d [Dec] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t a1, a2;
        mpz_t a_cp, a_csp;
        mpz_inits(a1, a2, a_cp, a_csp, NULL);
        start_time = omp_get_wtime();
        cp.PDec(a1, ea);
        csp.PDec(a2, ea);
        cp.TDec(a_cp, a1, a2);
        end_time = omp_get_wtime();
        average_time_PDec_TDec += (end_time - start_time);
//        printf("%d [PDec and TDec] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t cz;
        mpz_init(cz);
        start_time = omp_get_wtime();
        mpz_mul(cz, ea, ea);
        mpz_mod(cz, cz, pai.pk.N_Square);
        end_time = omp_get_wtime();
        average_time_HE_Addition += (end_time - start_time);
//        printf("%d [HE_Addition] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t ten;
        mpz_init(ten);
        mpz_set_si(ten, 10);
        start_time = omp_get_wtime();
        mpz_powm(cz, ea, ten, pai.pk.N_Square);
        end_time = omp_get_wtime();
        average_time_HE_ScalarMul += (end_time - start_time);
//        printf("%d [HE_ScalarMul] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t neg_one;
        mpz_init(neg_one);
        mpz_set_si(neg_one, -1);
        start_time = omp_get_wtime();
        mpz_powm(cz, ea, neg_one, pai.pk.N_Square);
        mpz_mul(cz, ea, ea);
        mpz_mod(cz, cz, pai.pk.N_Square);
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
    mpz_rrandomb(x, randstate, 8);
    mpz_rrandomb(y, randstate, 8);
    pai.encrypt(ex, x);
    pai.encrypt(ey, y);


    // 初始化CP
    int sock_to_cp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8081);
    mpz_t L_k;
    mpz_init(L_k);
    mpz_set_si(L_k, keyGen.pk.L_k);

    // 发给CP <sk1, pk , ex, ey>
    gmp_printf("Send (PK,SK) to CP, sk1 = %Zd\n", psk[0].partial_key);
    send_mpz(sock_to_cp, psk[0].partial_key);
    send_mpz(sock_to_cp, keyGen.pk.N); // N
    send_mpz(sock_to_cp, keyGen.pk.h); // h
    send_mpz(sock_to_cp, keyGen.pk.h_N); // h_N
    send_mpz(sock_to_cp, keyGen.pk.L); // L
    send_mpz(sock_to_cp, L_k); //L_k
    send_mpz(sock_to_cp, ex);
    send_mpz(sock_to_cp, ey);
    send_mpz(sock_to_cp, x); //用于验证
    send_mpz(sock_to_cp, y); //用于验证
    send_mpz(sock_to_cp, keyGen.sk.alpha); //用于验证

    // 初始化CSP
    int sock_to_csp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8082);

    // 发给CSP <sk2, pk>
    gmp_printf("Send (PK,SK) to CSP, sk2 = %Zd\n", psk[1].partial_key);
    send_mpz(sock_to_csp, psk[1].partial_key);
    send_mpz(sock_to_csp, keyGen.pk.N); // N
    send_mpz(sock_to_csp, keyGen.pk.h); // h
    send_mpz(sock_to_csp, keyGen.pk.h_N); // h_N
    send_mpz(sock_to_csp, keyGen.pk.L); // L
    send_mpz(sock_to_csp, L_k); //L_k

    gmp_printf("------ N = %Zd\n", keyGen.pk.N);
    gmp_printf("------ h = %Zd\n", keyGen.pk.h);
    gmp_printf("------ h_N = %Zd\n", keyGen.pk.h_N);
    gmp_printf("------ L = %Zd\n", keyGen.pk.L);
    gmp_printf("------ L_k = %Zd\n", L_k);
    gmp_printf("------ x = %Zd\n", x);
    gmp_printf("------ y = %Zd\n", y);
    gmp_printf("------ ex = %Zd\n", ex);
    gmp_printf("------ ey = %Zd\n", ey);

    // 收到结束标志
    close(sock_to_cp); // printf("关闭CP的连接\n");
    close(sock_to_csp); // printf("关闭CSP的连接\n");
}