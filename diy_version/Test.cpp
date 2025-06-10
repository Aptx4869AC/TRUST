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

void division(double &res, mpz_t &value, long long multiple)
{
    if (multiple == 0)
    {
        cerr << "错误：除数不能为零！" << '\n';
        return;
    }
    mpz_t integer, remainder, fractional;
    mpz_inits(integer, remainder, fractional, NULL);

    if (mpz_cmp_si(value, 0) > 0)
    {
        mpz_div_ui(integer, value, multiple);

        mpz_t mod;
        mpz_init(mod);
        mpz_init_set_si(mod, multiple);
        mpz_mod(fractional, value, mod);

        // 转换为 int
        long long integer_val = mpz_get_si(integer);
        double fractional_val = mpz_get_si(fractional) / (multiple * 1.0);
        res = integer_val + fractional_val;
        mpz_clear(mod);
    } else
    {
        mpz_t neg_value;
        mpz_init(neg_value);
        mpz_neg(neg_value, value);

        mpz_div_ui(integer, neg_value, multiple);

        mpz_t mod;
        mpz_init(mod);
        mpz_init_set_si(mod, multiple);
        mpz_mod(fractional, neg_value, mod);

        // 转换为 int
        long long integer_val = mpz_get_si(integer);
        double fractional_val = mpz_get_si(fractional) / (multiple * 1.0);
        res = integer_val + fractional_val;
        res *= -1;
        mpz_clear(mod);
        mpz_clear(neg_value);
    }
}

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


    printf("\n----------------------- 案例7 【负数支持】 -----------------------\n");
    // 本质是模空间的代数性质（同余定理）
    mpz_t neg_m;
    mpz_init(neg_m);
    mpz_neg(neg_m, m);
    gmp_printf("-m = %Zd\n", neg_m);
    pai.encrypt(eres, neg_m); // 负数处理已经隐藏在 encrypt()
    gmp_printf("[-m] --- eres = %Zd\n", eres);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    printf("\n----------------------- 案例8 【浮点数支持】 -----------------------\n");
    // 本质是定点数编码
    int alpha = 4; // 精度等级
    mpz_t ten, eps, multiEps, inv_multiEps;
    mpz_inits(ten, eps, multiEps, inv_multiEps, NULL);
    mpz_set_si(ten, 10);
    mpz_set_si(eps, 2 * alpha);
    mpz_powm(multiEps, ten, eps, pai.pk.N_Square);
    mpz_invert(inv_multiEps, multiEps, pai.pk.N);
    printf("alpha = %d\n", alpha);
    gmp_printf("eps = %Zd\n", eps);
    gmp_printf("multiEps = %Zd\n", multiEps);
    gmp_printf("inv_multiEps = %Zd\n", inv_multiEps);

    mpz_t m1, m2;
    mpz_inits(m1, m2, NULL);
    double float_m1 = 4869.2345;
    mpz_set_si(m1, float_m1 * pow(10, 2 * alpha));
    gmp_printf("m1 = %Zd\n", m1);
    double float_m2 = 1412.6789;
    mpz_set_si(m2, float_m2 * pow(10, 2 * alpha));
    gmp_printf("m2 = %Zd\n", m2);

    mpz_t em1, em2;
    mpz_inits(em1, em2, NULL);
    pai.encrypt(em1, m1);
    pai.encrypt(em2, m2);
    gmp_printf("[m1] --- eres = %Zd\n", em1);
    gmp_printf("[m2] --- eres = %Zd\n", em2);
    mpz_mul(eres, em1, em2);
    mpz_mod(eres, eres, pai.pk.N_Square);
    pai.decrypt(res, eres);
    gmp_printf("res = %Zd\n", res);

    double value;
    division(value, res, pow(10, 2 * alpha));
    printf("m1 + m2 = %.4f\n", value);

    printf("\n\n");


    return 0;
}
