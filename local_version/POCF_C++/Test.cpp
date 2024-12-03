/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-01
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include "include/pcpd.h"
#include "include/pocf.h"

using namespace std;
using namespace PHESPACE;

static int key_len = 512; // 256, 512, 768, 1024

int main()
{

    double start_time, end_time;
    double average_time_keygen = 0;
    double average_time_enc = 0;
    double average_time_dec = 0;
    double average_time_HE_Addition = 0;
    double average_time_HE_ScalarMul = 0;
    double average_time_HE_Subtraction = 0;
    double average_time_PDec_TDec = 0;
    double average_time_SMUL = 0;
    double average_time_SCMP = 0;
    double average_time_SSBA = 0;
    double average_time_SDIV = 0;
    double average_time_SXOR = 0;
    double average_time_SEQ = 0;

    Paillier pai;
    pai.KeyGen(key_len, sigma, 0); // k, sigma, yita

    start_time = omp_get_wtime();
    Paillier temp_pai;

    int epoch = 30;
    for (int i = 0; i < 5; i++)
    {
        temp_pai.KeyGen(key_len, sigma, 0); // k, sigma, yita
    }
    end_time = omp_get_wtime();
    average_time_keygen += (end_time - start_time);

    Paillier_Third_Party cp(pai.public_key, pai.private_key_1);
    Paillier_Third_Party csp(pai.public_key, pai.private_key_2);
    for (int i = 0; i < epoch; i++)
    {
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

        mpz_t x, y, ex, ey, eres, res;
        mpz_inits(x, y, ex, ey, eres, res, NULL);
        mpz_rrandomb(x, randstate, 6);
        mpz_rrandomb(y, randstate, 8);
//        gmp_printf("x = %Zd\n", x);
//        gmp_printf("y = %Zd\n", y);
        Enc(ex, pai.public_key, x);
        Enc(ey, pai.public_key, y);

        start_time = omp_get_wtime();
        SMUL(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        check_smul(i, eres, x, y, pai);
        average_time_SMUL += (end_time - start_time);
//        printf("%d [SMUL] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);


        start_time = omp_get_wtime();
        SCMP(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        check_scmp(i, eres, x, y, pai);
        average_time_SCMP += (end_time - start_time);
//        printf("%d [SCMP] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);


        mpz_t s_x, u_x, s, u, temp_s, temp_u;
        mpz_inits(s_x, u_x, s, u, temp_s, temp_u, NULL);
        start_time = omp_get_wtime();
        SSBA(s_x, u_x, ex, cp, csp);
        end_time = omp_get_wtime();
        check_ssba(i, s_x, u_x, x, pai);
        average_time_SSBA += (end_time - start_time);
//        printf("%d [SSBA] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t eq, ee, q, r, right_q, right_r;
        mpz_inits(eq, ee, q, r, right_q, right_r, NULL);
        start_time = omp_get_wtime();
        SDIV(eq, ee, ex, ey, 10, cp, csp, pai);
        end_time = omp_get_wtime();
        check_sdiv(i, eq, ee, x, y, pai);
        average_time_SDIV += (end_time - start_time);
//        printf("%d [SDIV] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);
//        printf("----------------------------------------------------------\n");


        start_time = omp_get_wtime();
        SXOR(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        average_time_SXOR += (end_time - start_time);

        start_time = omp_get_wtime();
        SEQ(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        check_seq(i, eres, x, y, pai);
        average_time_SEQ += (end_time - start_time);

        mpz_clears(a, ea, a1, a2, a_cp, a_csp, x, y, ex, ey, eres, res, NULL);
    }


    printf("\n\n");
    printf("KeyGen (%d average) time is  ------  %f ms\n", 5, average_time_keygen / 5 * 1000);
    printf("Enc (%d average) time is  ------  %f ms\n", epoch, average_time_enc / epoch * 1000);
    printf("Dec (%d average) time is  ------  %f ms\n", epoch, average_time_dec / epoch * 1000);
    printf("PDec + TDec (%d average) time is  ------  %f ms\n", epoch, average_time_PDec_TDec / epoch * 1000);
    printf("HE_Addition (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Addition / epoch * 1000);
    printf("HE_ScalarMul (%d average) time is  ------  %f ms\n", epoch, average_time_HE_ScalarMul / epoch * 1000);
    printf("HE_Subtraction (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Subtraction / epoch * 1000);
    printf("SMUL (%d average) time is  ------  %f ms\n", epoch, average_time_SMUL / epoch * 1000);
    printf("SCMP (%d average) time is  ------  %f ms\n", epoch, average_time_SCMP / epoch * 1000);
    printf("SSBA (%d average) time is  ------  %f ms\n", epoch, average_time_SSBA / epoch * 1000);
    printf("SDIV (10-bits) (%d average) time is  ------  %f ms\n", epoch, average_time_SDIV / epoch * 1000);
    printf("SXOR (%d average) time is  ------  %f ms\n", epoch, average_time_SXOR / epoch * 1000);
    printf("SEQ (%d average) time is  ------  %f ms\n", epoch, average_time_SEQ / epoch * 1000);

    return 0;
}
