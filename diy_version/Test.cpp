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

    k = 80; // 控制安全强度
    sigma = 128; // 控制比特串大小，如随机数
    setrandom(&randstate);
    KeyGen keyGen(k, sigma);
    Paillier pai(keyGen.pk, keyGen.sk, sigma);
    PaillierThirdParty *psk = keyGen.pai_third_parties;
    PaillierThdDec cp = PaillierThdDec(psk[0].N, psk[0].partial_key, keyGen.pk, sigma);
    PaillierThdDec csp = PaillierThdDec(psk[1].N, psk[1].partial_key, keyGen.pk, sigma);

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

    printf("\n----------------------- 案例1 【加法同态性】 -----------------------\n");
    mpz_mul(eres, em, ezero);
    mpz_mod(eres, eres, pai.pk.N_Square);
    gmp_printf("[m] \cdot [0] --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    printf("\n----------------------- 案例2 【加法同态性】 -----------------------\n");
    mpz_mul(eres, em, eone);
    mpz_mod(eres, eres, pai.pk.N_Square);
    gmp_printf("[m] \cdot [1] --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    printf("\n----------------------- 案例3 【加法同态性】 -----------------------\n");
    mpz_mul(eres, em, eneg_one);
    mpz_mod(eres, eres, pai.pk.N_Square);
    gmp_printf("[m] \cdot [-1] --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);


    printf("\n----------------------- 案例4 【标量乘法同态性】 -----------------------\n");
    mpz_powm(eres, em, zero, pai.pk.N_Square);
    gmp_printf("[m]^0 --- eres = %Zd\n", eres); // 由此可见密文是1，于是需要补充 \cdot Enc([0]), 刷新密文
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);
    mpz_mul(eres, eres, ezero);
    mpz_mod(eres, eres, pai.pk.N_Square);
    gmp_printf("[res] \cdot [0] --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    printf("\n----------------------- 案例5 【标量乘法同态性】 -----------------------\n");
    mpz_powm(eres, em, one, pai.pk.N_Square);
    gmp_printf("[m]^1 --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    printf("\n----------------------- 案例6 【标量乘法同态性】 -----------------------\n");
    mpz_powm(eres, em, neg_one, pai.pk.N_Square);
    gmp_printf("[m]^{-1} --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    printf("\n----------------------- 案例7 【浮点数支持】 -----------------------\n");
    int eps = 10; // 精度等级

//   TODO
//    mpz_powm(eres, em, neg_one, pai.pk.N_Square);
//    gmp_printf("[m]^{-1} --- eres = %Zd\n", eres);
//    pai.decrypt(res, eres);
//    gmp_printf("res = %Zd\n", res);

    printf("\n\n");


    return 0;
}
