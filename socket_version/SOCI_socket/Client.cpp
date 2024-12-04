/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include "include/paillier.h"
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

    int epoch = 10;
    int key_len = 1024; // 256, 512, 768, 1024, 1280, 1536
    printf("epoch = %d\n", epoch);
    printf("key_len = %d\n", key_len);
    printf("N 的比特长度为 %d bits\n", 2 * key_len);

    Paillier pai;
    pai.keygen(key_len);

    // 密钥拆分
    mpz_t sk1, sk2;
    mpz_inits(sk1, sk2, NULL);
    setrandom(&randstate);
    mpz_rrandomb(sk1, randstate, sigma);    // sk1 is a ranodm number with sigma bits
    mpz_mul(sk2, pai.prikey.lambda, pai.prikey.lmdInv);
    mpz_sub(sk2, sk2, sk1);                // sk2 = lambda ?mu - sk1
    printf("the bits of sk1 = %ld\n", mpz_sizeinbase(sk1, 2));
    printf("the bits of sk2 = %ld\n", mpz_sizeinbase(sk2, 2));

    printf("----------------------------------------------------------\n");
    start_time = omp_get_wtime();
    for (int i = 0; i < epoch; i++) {
        Paillier temp;
        temp.keygen(key_len);
    }
    end_time = omp_get_wtime();
    printf("Keygen (%d average) time is  ------  %f ms\n", epoch, (end_time - start_time) / epoch * 1000);

//    printf("----------------------------------------------------------\n");
//    start_time = omp_get_wtime();
//    for (int i = 0; i < epoch; i++) {
//        PaillierThdPrivateKey *temp_psk = thdkeygen(pai.prikey, sigma);
//        PaillierThdDec temp_cp(temp_psk[0], pai.pubkey);
//        PaillierThdDec temp_csp(temp_psk[1], pai.pubkey);
//    }
//    end_time = omp_get_wtime();
//    printf("Thdkeygen (%d average) time is  ------  %f ms\n", epoch, (end_time - start_time) / epoch * 1000);

    printf("----------------------------------------------------------\n");

    mpz_t x, y, ex, ey, c1, c2;
    mpz_inits(x, y, ex, ey, c1, c2, NULL);
    setrandom(&randstate);
    mpz_rrandomb(x, randstate, 8); // 8-bit random number
    mpz_rrandomb(y, randstate, 8); // 8-bit random number
    pai.encrypt(ex, x);
    pai.encrypt(ey, y);
    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        mpz_rrandomb(cz, randstate, 8); // 8-bit random number
        start_time = omp_get_wtime();
        pai.encrypt(cz, cz);
        end_time = omp_get_wtime();
        average_time_enc += (end_time - start_time);
        mpz_clear(cz);

    }
    printf("Enc (%d average) time is  ------  %f ms\n", epoch, average_time_enc / epoch * 1000);

    printf("----------------------------------------------------------\n");


    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        start_time = omp_get_wtime();
        pai.decrypt(cz, ex);
        end_time = omp_get_wtime();
        average_time_dec += (end_time - start_time);
        if (mpz_cmp(cz, x) != 0) {
            gmp_printf("x = %Zd\n", y);
            gmp_printf("dec(x) = %Zd\n", cz);
            printf("Dec is error!\n");
            exit(-1);
        }
        mpz_clear(cz);
    }


    printf("Dec (%d average) time is  ------  %f ms\n", epoch, average_time_dec / epoch * 1000);

    printf("----------------------------------------------------------\n");
    // 创建临时副本来测PDec
    PaillierThdPrivateKey share_part1(sk1, pai.prikey.n, pai.prikey.nsquare);
    PaillierThdPrivateKey share_part2(sk2, pai.prikey.n, pai.prikey.nsquare);
    PaillierThdDec cp(share_part1, pai.pubkey);
    PaillierThdDec csp(share_part2, pai.pubkey);
    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        start_time = omp_get_wtime();
        cp.pdec(c1, ey);
        csp.pdec(c2, ey);
        cp.fdec(cz, c1, c2);
        end_time = omp_get_wtime();
        average_time_PDec_TDec += (end_time - start_time);
        if (mpz_cmp(cz, y) != 0) {
            gmp_printf("y = %Zd\n", y);
            gmp_printf("dec(y) = %Zd\n", cz);
            printf("PDec + TDec is error!\n");
            exit(-1);
        }
        mpz_clear(cz);

    }
    printf("PDec + TDec (%d average) time is  ------  %f ms\n", epoch, average_time_PDec_TDec / epoch * 1000);

    printf("----------------------------------------------------------\n");

    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        // x+y
        start_time = omp_get_wtime();
        mpz_mul(cz, ex, ey);
        mpz_mod(cz, cz, pai.pubkey.nsquare);
        end_time = omp_get_wtime();
        average_time_HE_Addition += (end_time - start_time);
        mpz_clear(cz);
    }

    printf("HE_Addition (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Addition / epoch * 1000);

    printf("----------------------------------------------------------\n");

    mpz_t ten;
    mpz_init(ten);
    mpz_set_si(ten, 10);
    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        // 10x
        start_time = omp_get_wtime();
        mpz_powm(cz, ex, ten, pai.pubkey.nsquare);
        end_time = omp_get_wtime();
        average_time_HE_ScalarMul += (end_time - start_time);
        mpz_clear(cz);
    }

    printf("HE_ScalarMul (%d average) time is  ------  %f ms\n", epoch, average_time_HE_ScalarMul / epoch * 1000);


    printf("----------------------------------------------------------\n");

    mpz_t neg_one;
    mpz_init(neg_one);
    mpz_set_si(neg_one, -1);
    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        // x-y
        start_time = omp_get_wtime();
        mpz_powm(cz, ey, neg_one, pai.pubkey.nsquare);
        mpz_mul(cz, ex, ey);
        mpz_mod(cz, cz, pai.pubkey.nsquare);
        end_time = omp_get_wtime();
        average_time_HE_Subtraction += (end_time - start_time);
        mpz_clear(cz);
    }

    printf("HE_Subtraction (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Subtraction / epoch * 1000);
    printf("----------------------------------------------------------\n");

    // 初始化CP
    int sock_to_cp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8081);

    // 发给CP <sk2, pai.prikey.n, ex, ey>
    gmp_printf("Send (PK,SK) to CP, sk1 = %Zd\n", sk1);
    send_mpz(sock_to_cp, sk1);
    send_mpz(sock_to_cp, pai.prikey.n);
    send_mpz(sock_to_cp, ex);
    send_mpz(sock_to_cp, ey);
    send_mpz(sock_to_cp, x); //用于验证
    send_mpz(sock_to_cp, y); //用于验证
    send_mpz(sock_to_cp, pai.prikey.lambda); //用于验证

    // 初始化CSP
    int sock_to_csp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8082);

    // 发给CSP <sk2, pai.prikey.n>
    gmp_printf("Send (PK,SK) to CSP, sk2 = %Zd\n", sk2);
    send_mpz(sock_to_csp, sk2);
    send_mpz(sock_to_csp, pai.prikey.n);

    gmp_printf("------ n = %Zd\n", pai.prikey.n);
    gmp_printf("------ nsquare = %Zd\n", pai.prikey.nsquare);
    gmp_printf("------ x = %Zd\n", x);
    gmp_printf("------ y = %Zd\n", y);
    gmp_printf("------ ex = %Zd\n", ex);
    gmp_printf("------ ey = %Zd\n", ey);

    // 收到结束标志
    close(sock_to_cp); // printf("关闭CP的连接\n");
    close(sock_to_csp); // printf("关闭CSP的连接\n");
    return 0;
}