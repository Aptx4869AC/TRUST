/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-31
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

    double average_time_Trust_FMUL = 0;
    double average_time_Trust_FCMP = 0;
    double average_time_Trust_FEQL = 0;
    double average_time_Trust_FABS = 0;
    double average_time_Trust_FTRN_v1 = 0;
    double average_time_Trust_FTRN_v2 = 0;

    k = 112; // 控制安全强度
    sigma = 128; // 控制比特串大小，如随机数
    setrandom(&randstate);
    KeyGen keyGen(k, sigma);
    Paillier pai(keyGen.pk, keyGen.sk, sigma);

    int epoch = 10;
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

        mpz_clears(a, ea, a1, a2, a_cp, a_csp, NULL);
    }



    int Trust_epoch = 10;
    /** 开销测试：Trust协议 **/
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

        /** 密文乘法 **/
        start_time = omp_get_wtime();
        sc.trust_fmul(eres, em_1, em_2, cp, csp);
        end_time = omp_get_wtime();
        check_fmul(i, eres, m_1, m_2, pai);
        average_time_Trust_FMUL += (end_time - start_time);

        /** 密文比较 **/
        start_time = omp_get_wtime();
        sc.trust_fcmp(eres, em_1, em_2, cp, csp);
        end_time = omp_get_wtime();
        check_fcmp(i, eres, m_1, m_2, pai);
        average_time_Trust_FCMP += (end_time - start_time);

        /** 密文判断相等 **/
        start_time = omp_get_wtime();
        sc.trust_feql(eres, em_1, em_2, cp, csp);
        end_time = omp_get_wtime();
        check_feql(i, eres, m_1, m_2, pai);
        average_time_Trust_FEQL += (end_time - start_time);

        /** 密文绝对值 **/
        start_time = omp_get_wtime();
        sc.trust_fabs(eres, em_1, cp, csp);
        end_time = omp_get_wtime();
        check_fabs(i, eres, m_1, pai);
        average_time_Trust_FABS += (end_time - start_time);

//        mpz_set_si(m_1, 1);
//        pai.encrypt(em_1, m_1);

        /** 三目运算，此版本是论文里面的 **/
        start_time = omp_get_wtime();
        sc.trust_ftrn_version1(eres, em_1, em_2, em_3, cp, csp);
        end_time = omp_get_wtime();
        check_ftrn(i, eres, m_1, m_2, m_3, pai);
        average_time_Trust_FTRN_v1 += (end_time - start_time);

        int pi = (int) random() % 2;
        mpz_set_si(m_1, pi);
        pai.encrypt(em_1, m_1);

        /** 三目运算，此版本是未写在论文里面的 **/
        start_time = omp_get_wtime();
        sc.trust_ftrn_version2(eres, em_1, em_2, em_3, cp, csp);
        end_time = omp_get_wtime();
        check_ftrn(i, eres, m_1, m_2, m_3, pai);
        average_time_Trust_FTRN_v2 += (end_time - start_time);


//        printf("------------------\n");
        mpz_clears(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
    }

    /** 特殊测试 **/
    mpz_t m, zero, one, neg_one, res, em, ezero, eone, eneg_one, eres;
    mpz_inits(m, zero, one, neg_one, res, em, ezero, eone, eneg_one, eres, NULL);
    mpz_rrandomb(m, randstate, 8);
    mpz_set_si(zero, 0);
    mpz_set_si(one, 1);
    mpz_set_si(neg_one, -1);
    pai.encrypt(em, m);
    pai.encrypt(ezero, zero);
    pai.encrypt(eone, one);
    pai.encrypt(eneg_one, neg_one);
    gmp_printf("m = %Zd\n", m);
    gmp_printf("em = %Zd\n", em);

    mpz_mul(eres, em, ezero);
    mpz_mod(eres, eres, pai.pk.N_Square);
    gmp_printf("[m] \cdot [0] --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    mpz_mul(eres, em, eone);
    mpz_mod(eres, eres, pai.pk.N_Square);
    gmp_printf("[m] \cdot [1] --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    mpz_mul(eres, em, eneg_one);
    mpz_mod(eres, eres, pai.pk.N_Square);
    gmp_printf("[m] \cdot [-1] --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);


    mpz_powm(eres, em, zero, pai.pk.N_Square);
    gmp_printf("[m]^0 --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);
    mpz_powm(eres, em, one, pai.pk.N_Square);
    gmp_printf("[m]^1 --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);
    mpz_powm(eres, em, neg_one, pai.pk.N_Square);
    gmp_printf("[m]^{-1} --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    printf("\n\n");
    printf("KeyGen (%d average) time is  ------  %f ms\n", 10, average_time_keygen / 10 * 1000);
    printf("Enc (%d average) time is  ------  %f ms\n", epoch, average_time_enc / epoch * 1000);
    printf("Dec (%d average) time is  ------  %f ms\n", epoch, average_time_dec / epoch * 1000);
    printf("PDec + TDec (%d average) time is  ------  %f ms\n", epoch, average_time_PDec_TDec / epoch * 1000);
    printf("HE_Addition (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Addition / epoch * 1000);
    printf("HE_ScalarMul (%d average) time is  ------  %f ms\n", epoch, average_time_HE_ScalarMul / epoch * 1000);
    printf("HE_Subtraction (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Subtraction / epoch * 1000);

    printf("Trust_FMUL (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FMUL / Trust_epoch * 1000);
    printf("Trust_FCMP (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FCMP / Trust_epoch * 1000);
    printf("Trust_FEQL (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FEQL / Trust_epoch * 1000);
    printf("Trust_FABS (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FABS / Trust_epoch * 1000);
    printf("Trust_FTRN_v1 (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FTRN_v1 / Trust_epoch * 1000);
    printf("Trust_FTRN_v2 (%d average) time is  ------  %f ms\n", Trust_epoch, average_time_Trust_FTRN_v2 / Trust_epoch * 1000);



    return 0;
}
