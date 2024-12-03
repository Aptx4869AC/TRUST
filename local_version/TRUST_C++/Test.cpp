/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-01
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include "include/protocol.h"
#include "include/fastPai.h"
#include <gmp.h>
#include <omp.h>

using namespace PROTOCOLSPACE;
using namespace PHESPACE;

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
    double average_time_TEE_SMUL = 0;
    double average_time_TEE_SCMP = 0;
    double average_time_TEE_SABS = 0;
    double average_time_TEE_SDIV = 0;
    double average_time_Trust_FMUL = 0;
    double average_time_Trust_FCMP = 0;
    double average_time_Trust_FEQL = 0;
    double average_time_Trust_FABS = 0;
    double average_time_Trust_FTRN_v1 = 0;
    double average_time_Trust_FTRN_v2 = 0;

    double average_time_8 = 0;
    double average_time_9 = 0;
    double average_time_10 = 0;
    double average_time_11 = 0;
    double average_time_12 = 0;

    double average_bytes_PK = 0;
    double average_bytes_SK = 0;

    k = 112; // 控制安全强度
    sigma = 128; // 控制比特串大小，如随机数
    setrandom(&randstate);
    KeyGen keyGen(k, sigma);
    Paillier pai(keyGen.pk, keyGen.sk, sigma);

    int epoch = 10;
    int TEE_epoch = 10;
    int Trust_epoch = 10;
    start_time = omp_get_wtime();
    for (int i = 0; i < 10; i++)
    {
        KeyGen temp_keyGen(k, sigma);
    }
    end_time = omp_get_wtime();
    average_time_keygen += (end_time - start_time);

    PaillierThirdParty *psk = keyGen.pai_third_parties;
    PaillierThdDec cp = PaillierThdDec(psk[0].N, psk[0].partial_key, keyGen.pk, sigma);
    PaillierThdDec csp = PaillierThdDec(psk[1].N, psk[1].partial_key, keyGen.pk, sigma);

    /**
     * 开销测试：：FastPai下基础加密原语、WQ的SOCI+离线在线结合、Zhao的SOCI协议在线版本
     */
    for (int i = 0; i < epoch; i++)
    {
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
        mpz_set_si(ten, 1000);
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

        protocol sc;
        mpz_t x, y, ex, ey, eres, res;
        mpz_inits(x, y, ex, ey, eres, res, NULL);
        mpz_rrandomb(x, randstate, 10);
        mpz_rrandomb(y, randstate, 10);
//        gmp_printf("x = %Zd\n", x);
//        gmp_printf("y = %Zd\n", y);
        pai.encrypt(ex, x);
        pai.encrypt(ey, y);

        start_time = omp_get_wtime();
        sc.smul(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        check_fmul(i, eres, x, y, pai);
        average_time_SMUL += (end_time - start_time);
//        printf("%d [SMUL] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        /*
        start_time = omp_get_wtime();
        sc.online_smul(eres, ex, ey, cp, csp); // SOCI版本，无刷新
        end_time = omp_get_wtime();
        printf("time is  ------  %f ms\n",  i,  (end_time - start_time) * 1000);*/


        start_time = omp_get_wtime();
        sc.scmp(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        check_fcmp(i, eres, x, y, pai);
        average_time_SCMP += (end_time - start_time);
//        printf("%d [SCMP] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        /*
        start_time = omp_get_wtime();
        sc.online_scmp(eres, ex, ey, cp, csp); // SOCI版本，无刷新
        pai.decrypt(res, eres);
        end_time = omp_get_wtime();
        printf("time is  ------  %f ms\n",  i,  (end_time - start_time) * 1000);*/

        mpz_t s_x, u_x, s, u, temp_s, temp_u;
        mpz_inits(s_x, u_x, s, u, temp_s, temp_u, NULL);
        start_time = omp_get_wtime();
        sc.ssba(s_x, u_x, ex, cp, csp);
        end_time = omp_get_wtime();
        check_ssba(i, s_x, u_x, x, pai);
        average_time_SSBA += (end_time - start_time);
//        printf("%d [SSBA] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);
        /*
        start_time = omp_get_wtime();
        sc.online_ssba(s_x, u_x, ex, cp, csp); // SOCI版本，无刷新
        end_time = omp_get_wtime();
        printf("time is  ------  %f ms\n",  i,  (end_time - start_time) * 1000);*/

        mpz_t eq, ee, q, r, right_q, right_r;
        mpz_inits(eq, ee, q, r, right_q, right_r, NULL);

        start_time = omp_get_wtime();
        sc.sdiv(eq, ee, ex, ey, 10, cp, csp);
        end_time = omp_get_wtime();
        check_sdiv(i, eq, ee, x, y, pai);
        average_time_SDIV += (end_time - start_time);
//        printf("%d [SDIV] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);
        /*
        start_time = omp_get_wtime();
        sc.online_sdiv(eq, ee, ex, ey, 10, cp, csp); // SOCI版本，无刷新
        end_time = omp_get_wtime();
        printf("time is  ------  %f ms\n",  i,  (end_time - start_time) * 1000);*/

//        printf("----------------------------------------------------------\n");

        mpz_clears(a, ea, a1, a2, a_cp, a_csp, x, y, ex, ey, eres, res, NULL);
    }


    /**
     * 开销测试：：SOCI-TEE协议
     */
    for (int i = 0; i < TEE_epoch; i++)
    {
        protocol sc;
        mpz_t x, y, ex, ey, eres, res;
        mpz_inits(x, y, ex, ey, eres, res, NULL);
        mpz_rrandomb(x, randstate, 8);
        mpz_rrandomb(y, randstate, 8);
        pai.encrypt(ex, x);
        pai.encrypt(ey, y);

        start_time = omp_get_wtime();
        sc.tee_smul(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        check_fmul(i, eres, x, y, pai);
        average_time_TEE_SMUL += (end_time - start_time);

        start_time = omp_get_wtime();
        sc.tee_scmp(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        check_fcmp(i, eres, x, y, pai);
        average_time_TEE_SCMP += (end_time - start_time);

        start_time = omp_get_wtime();
        sc.tee_sabs(eres, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        check_sabs(i, eres, x, y, pai);
        average_time_TEE_SABS += (end_time - start_time);

        mpz_t eq, ee, q, r, right_q, right_r;
        mpz_inits(eq, ee, q, r, right_q, right_r, NULL);

        start_time = omp_get_wtime();
        sc.tee_sdiv(eq, ee, ex, ey, 10, cp, csp);
        end_time = omp_get_wtime();
        check_sdiv(i, eq, ee, x, y, pai);
        average_time_TEE_SDIV += (end_time - start_time);
    }

    /**
     * 开销测试：：Trust协议
     */
    for (int i = 0; i < Trust_epoch; i++)
    {
        protocol sc;
        mpz_t m_1, m_2, m_3, em_1, em_2, em_3, eres, res;
        mpz_inits(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
        mpz_rrandomb(m_1, randstate, 8);
        mpz_rrandomb(m_2, randstate, 8);
        do
        {
            mpz_rrandomb(m_3, randstate, 8);
        } while (mpz_cmp(m_2, m_3) == 0);
        pai.encrypt(em_1, m_1);
        pai.encrypt(em_2, m_2);
        pai.encrypt(em_3, m_3);
//        gmp_printf("m_1 = %Zd\n", m_1);
//        gmp_printf("m_2 = %Zd\n", m_2);
//        gmp_printf("m_3 = %Zd\n", m_3);

        // MUL
        start_time = omp_get_wtime();
        sc.trust_fmul(eres, em_1, em_2, cp, csp);
        end_time = omp_get_wtime();
        check_fmul(i, eres, m_1, m_2, pai);
        average_time_Trust_FMUL += (end_time - start_time);

        // CMP
        start_time = omp_get_wtime();
        sc.trust_fcmp(eres, em_1, em_2, cp, csp);
        end_time = omp_get_wtime();
        check_fcmp(i, eres, m_1, m_2, pai);
        average_time_Trust_FCMP += (end_time - start_time);

        // EQL
        start_time = omp_get_wtime();
        sc.trust_feql(eres, em_1, em_2, cp, csp);
        end_time = omp_get_wtime();
        check_feql(i, eres, m_1, m_2, pai);
        average_time_Trust_FEQL += (end_time - start_time);

        // ABS
        start_time = omp_get_wtime();
        sc.trust_fabs(eres, em_1, cp, csp);
        end_time = omp_get_wtime();
        check_fabs(i, eres, m_1, pai);
        average_time_Trust_FABS += (end_time - start_time);

//        mpz_set_si(m_1, 1);
//        pai.encrypt(em_1, m_1);

        // TRN，此版本是论文里面的
        start_time = omp_get_wtime();
        sc.trust_ftrn_version1(eres, em_1, em_2, em_3, cp, csp);
        end_time = omp_get_wtime();
        check_ftrn(i, eres, m_1, m_2, m_3, pai);
        average_time_Trust_FTRN_v1 += (end_time - start_time);

        int pi = (int) random() % 2;
        mpz_set_si(m_1, pi);
        pai.encrypt(em_1, m_1);

        // TRN，此版本是未写在论文里面的
        start_time = omp_get_wtime();
        sc.trust_ftrn_version2(eres, em_1, em_2, em_3, cp, csp);
        end_time = omp_get_wtime();
        check_ftrn(i, eres, m_1, m_2, m_3, pai);
        average_time_Trust_FTRN_v2 += (end_time - start_time);


//        printf("------------------\n");
        mpz_clears(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
    }

    /**
     * 以下所注释代码是本人参考 POCF 中的新协议进行的尝试，供开发使用。逆运算部分需要使用扩展欧几里得算法，目前尚未完成。SGCD部分已修成了安全版本。
     */

    /*
    int Nepoch = 1;
    for (int i = 0; i < Nepoch; i++) {

        protocol sc;
        mpz_t x, y, m, ex, ey, em, eres, res;
        mpz_inits(x, y, m, ex, ey, em, eres, res, NULL);
        mpz_rrandomb(x, randstate, 3);
        mpz_rrandomb(y, randstate, 10);
        mpz_rrandomb(m, randstate, 16);
        // gmp_printf("x = %Zd\n", x);
        // mpz_neg(x, x);
        // gmp_printf("x = %Zd\n", x);
        mpz_set_si(x, 5);
        mpz_set_si(m, 117);
        // mpz_set_si(y, 11);
        pai.encrypt(ex, x);
        pai.encrypt(ey, y);
        pai.encrypt(em, m);


        mpz_t xinv, exinv, right_xinv;
        mpz_inits(xinv, exinv, right_xinv, NULL);
        start_time = omp_get_wtime();
        sc.online_sinv(exinv, ex, em, 10, cp, csp, pai);
        end_time = omp_get_wtime();
        average_time_8 += (end_time - start_time);
        printf("[SINV] time is  ------  %f ms\n", (end_time - start_time) * 1000);
        pai.decrypt(xinv, exinv);
        if (mpz_sgn(xinv) == -1) {
            mpz_mod(xinv, xinv, m);
        }
        mpz_invert(right_xinv, x, m);
        if (mpz_cmp(xinv, right_xinv) != 0) {
            printf("error!\n");
            gmp_printf("N = %Zd\n", csp.pai.pk.N);

            gmp_printf("[SINV] The modular inverse of %Zd mod N is: %Zd\n", x, xinv);
            gmp_printf("[right] The modular inverse of %Zd mod N is: %Zd\n", x, right_xinv);
            exit(-1);
        } else {
            printf("[SINV] is correct\n");
        }

        mpz_t ex_y, x_y, right_x_y;
        mpz_inits(ex_y, x_y, right_x_y, NULL);
        start_time = omp_get_wtime();
        sc.online_sexp(ex_y, x, ey, cp, csp);
        end_time = omp_get_wtime();
        average_time_9 += (end_time - start_time);
        printf("[SEXP] time is  ------  %f ms\n", (end_time - start_time) * 1000);
        pai.decrypt(x_y, ex_y);

        if (mpz_sgn(x_y) == -1) {
            // printf("sign: %d\n", mpz_sgn(x_y));
            mpz_add(x_y, x_y, cp.pai.pk.N);
        }

        mpz_powm(right_x_y, x, y, cp.pai.pk.N);


        if (mpz_cmp(x_y, right_x_y) != 0) {
            printf("error!\n");
            gmp_printf("x^y = %Zd\n", x_y);
            gmp_printf("x^y = %Zd\n", right_x_y);

            exit(-1);
        } else {
            printf("[SEXP] is correct\n");
        }

        // mpz_set_si(x, 64);
        // mpz_set_si(y, 64);
        // pai.encrypt(ex, x);
        // pai.encrypt(ey, y);
        mpz_t f, ef, right_f;
        mpz_inits(f, ef, right_f, NULL);
        start_time = omp_get_wtime();
        sc.online_seqc(ef, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        average_time_10 += (end_time - start_time);
        printf("[SEQC] time is  ------  %f ms\n", (end_time - start_time) * 1000);
        pai.decrypt(f, ef);
        if (mpz_cmp(x, y) == 0 && mpz_cmp_si(f, 0) != 0) {
            printf("error!\n");
            gmp_printf("x = %Zd\n", x);
            gmp_printf("y = %Zd\n", y);
            gmp_printf("x==y? %Zd\n", f);

            exit(-1);
        } else {
            printf("[SEQC] is correct\n");
        }

        mpz_t A, I, eA, eI;
        mpz_inits(A, I, NULL);
        start_time = omp_get_wtime();
        sc.online_smms(eA, eI, ex, ey, cp, csp);
        end_time = omp_get_wtime();
        average_time_11 += (end_time - start_time);
        printf("[SMMS] time is  ------  %f ms\n", (end_time - start_time) * 1000);
        pai.decrypt(A, eA);
        pai.decrypt(I, eI);

        if (mpz_sgn(A) == -1) {
            mpz_add(A, A, cp.pai.pk.N);
        }
        if (mpz_sgn(I) == -1) {
            mpz_add(I, I, cp.pai.pk.N);
        }

        if (mpz_cmp(x, y) == 0 && mpz_cmp(A, x) == 0 && mpz_cmp(I, y) == 0) {
            printf("[SMMS] is correct\n");

        } else if (mpz_cmp(x, y) < 0 && mpz_cmp(A, y) == 0 && mpz_cmp(I, x) == 0) {
            printf("[SMMS] is correct\n");
        } else if (mpz_cmp(x, y) > 0 && mpz_cmp(A, x) == 0 && mpz_cmp(I, y) == 0) {
            printf("[SMMS] is correct\n");
        } else {
            printf("error!\n");
            gmp_printf("A = %Zd\n", A);
            gmp_printf("I = %Zd\n", I);

            exit(-1);
        }


        mpz_rrandomb(x, randstate, 10);
        mpz_rrandomb(y, randstate, 6);
        pai.encrypt(ex, x);
        pai.encrypt(ey, y);
        mpz_t c, ec, right_c;
        mpz_inits(c, ec, right_c, NULL);
        start_time = omp_get_wtime();
        sc.online_sgcd(ec, ex, ey, mpz_sizeinbase(x, 2), cp, csp);
        end_time = omp_get_wtime();
        average_time_12 += (end_time - start_time);
        printf("[SGCD] time is  ------  %f ms\n", (end_time - start_time) * 1000);
        pai.decrypt(c, ec);
        mpz_gcd(right_c, x, y);
        if (mpz_cmp(c, right_c) != 0) {
            printf("error!\n");
            gmp_printf("gcd(x,y) = %Zd\n", c);
            gmp_printf("gcd(x,y) = %Zd\n", right_c);

            exit(-1);
        } else {
            printf("[SGCD] is correct\n");
        }

        printf("----------------------------------------------------------\n");

    }
    */


    printf("\n\n");
    printf("KeyGen (%d average) time is  ------  %f ms\n", epoch, average_time_keygen / epoch * 1000);
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

    printf("[version1] TEE_SMUL (%d average) time is  ------  %f ms\n", TEE_epoch, average_time_TEE_SMUL / TEE_epoch * 1000);
    printf("[version1] TEE_SCMP (%d average) time is  ------  %f ms\n", TEE_epoch, average_time_TEE_SCMP / TEE_epoch * 1000);
    printf("[version1] TEE_SABS (%d average) time is  ------  %f ms\n", TEE_epoch, average_time_TEE_SABS / TEE_epoch * 1000);
    printf("[version1] TEE_SDIV (%d average) time is  ------  %f ms\n", TEE_epoch, average_time_TEE_SDIV / TEE_epoch * 1000);

    printf("[version2] Trust_FMUL (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FMUL / Trust_epoch * 1000);
    printf("[version2] Trust_FCMP (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FCMP / Trust_epoch * 1000);
    printf("[version2] Trust_FEQL (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FEQL / Trust_epoch * 1000);
    printf("[version2] Trust_FABS (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FABS / Trust_epoch * 1000);
    printf("[version2] Trust_FTRN_v1 (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FTRN_v1 / Trust_epoch * 1000);
    printf("[version2] Trust_FTRN_v2 (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FTRN_v2 / Trust_epoch * 1000);

    /*
    printf("SINV (%d average) time is  ------  %f ms\n", Nepoch, average_time_8 / Nepoch * 1000);
    printf("SEXP (%d average) time is  ------  %f ms\n", Nepoch, average_time_9 / Nepoch * 1000);
    printf("SEQC (%d average) time is  ------  %f ms\n", Nepoch, average_time_10 / Nepoch * 1000);
    printf("SMMS (%d average) time is  ------  %f ms\n", Nepoch, average_time_11 / Nepoch * 1000);
    printf("SGCD (10-bits) (%d average) time is  ------  %f ms\n", Nepoch, average_time_12 / Nepoch * 1000);*/


    return 0;
}
