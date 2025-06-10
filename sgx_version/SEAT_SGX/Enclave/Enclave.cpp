/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-5-27
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include "Enclave.h"
#include "Enclave_t.h" /* print_string */
#include "sgx_trts.h"
#include <stdarg.h>
#include <stdio.h> /* vsnprintf */
#include <string.h>
#include "TrustedLibrary/sgx_tgmp.h"


gmp_randstate_t gmp_rand;
mpz_t ezero, eone;
mpz_t partial_key, N, N_Square, half_N, h_N;
int L_k;
int cnt;
mpz_t **table;
int rows;
int columns;

int printf(const char *fmt, ...)
{
    char buf[BUFSIZ] = {'\0'};
    va_list ap;
    va_start(ap, fmt);
    vsnprintf(buf, BUFSIZ, fmt, ap);
    va_end(ap);
    ocall_print_string(buf);
    return (int) strnlen(buf, BUFSIZ - 1) + 1;
}

void serialize_mpz(mpz_t num, unsigned char **buffer, size_t *size)
{
    *size = (mpz_sizeinbase(num, 2) + 7) / 8; // 计算所需字节数
    //*buffer = new unsigned char[*size];
    mpz_export(*buffer, size, 1, 1, 0, 0, num); // 导出为字节数组
}

void deserialize_mpz(mpz_t res, unsigned char *buffer, size_t size)
{
    mpz_import(res, size, 1, 1, 0, 0, buffer); // 从字节数组导入
}

void add(mpz_t dst, mpz_t ex, mpz_t ey)
{
    mpz_mul(dst, ex, ey);
    mpz_mod(dst, dst, N_Square);
}

void scl_mul(mpz_t dst, mpz_t c, mpz_t cons)
{
    mpz_powm(dst, c, cons, N_Square);
}

void sub(mpz_t dst, mpz_t c1, mpz_t c2)
{
    mpz_t neg_one;
    mpz_init(neg_one);
    mpz_set_si(neg_one, -1);
    scl_mul(dst, c2, neg_one);
    add(dst, dst, c1);
    mpz_clears(neg_one, NULL);
}

void L_function(mpz_t dst, mpz_t x, mpz_t N)
{
    mpz_sub_ui(dst, x, 1);
    mpz_div(dst, dst, N);
}

void PDec(mpz_t dst, mpz_t c)
{
    mpz_powm(dst, c, partial_key, N_Square);
}

void TDec(mpz_t dst, mpz_t c1, mpz_t c2)
{
    mpz_mul(dst, c1, c2);
    mpz_mod(dst, dst, N_Square);
    L_function(dst, dst, N);
}

void fastPower(mpz_t result, const mpz_t base, const mpz_t exponent, const mpz_t modulus)
{
    mpz_set_si(result, 1);  // 初始化结果为1
    mpz_t base_copy;
    mpz_init_set(base_copy, base);  // 复制底数用于计算

    mpz_t exp;
    mpz_init_set(exp, exponent);  // 用于迭代计算指数

    while (mpz_cmp_ui(exp, 0) > 0)
    {
        if (mpz_odd_p(exp))
        {
            mpz_mul(result, result, base_copy);  // 如果指数为奇数，累乘当前的底数到结果
            mpz_mod(result, result, modulus);  // 对结果取模
        }
        mpz_mul(base_copy, base_copy, base_copy);  // 底数平方
        mpz_mod(base_copy, base_copy, modulus);  // 对底数取模
        mpz_fdiv_q_2exp(exp, exp, 1);  // 指数除以2
    }

    mpz_clear(base_copy);
    mpz_clear(exp);
}


int blockSize = 5;

mpz_t **constructTable(mpz_t g, mpz_t nSquare, int l)
{
    mpz_t b, base;
    mpz_inits(b, base, NULL);
    mpz_set_si(b, blockSize);
    int lDivB = (l + blockSize - 1) / blockSize;
    mpz_set(base, g);
    auto **table = new mpz_t *[lDivB];
    int size = 1 << blockSize;
    for (int i = 0; i < lDivB; i++)
    {
        table[i] = new mpz_t[size];
        mpz_t tmp;
        mpz_init(tmp);
        for (int j = 0; j < size; j++)
        {
            mpz_init(table[i][j]);
            mpz_mul_ui(tmp, b, i);
            mpz_mod(tmp, tmp, nSquare);
            mpz_t two;
            mpz_init(two);
            mpz_set_si(two, 2);
            mpz_powm(tmp, two, tmp, nSquare);
            mpz_powm(tmp, base, tmp, nSquare);
            mpz_powm_ui(tmp, tmp, j, nSquare);

            mpz_set(table[i][j], tmp);
        }
        mpz_clear(tmp);
    }
    mpz_clears(b, base, NULL);
    rows = lDivB;
    columns = size;
    return table;
}

int *convertIntoBlock(mpz_t num)
{
    int b = blockSize;
    int l = mpz_sizeinbase(num, 2);
    cnt = (l + b - 1) / b;
    char *num_in_bit = mpz_get_str(nullptr, 2, num);
    int filllen = 0;
    if (l % b != 0) filllen = b - l % b;
    char *fill0 = new char[filllen + 1];
    for (int i = 0; i < filllen; i++) fill0[i] = '0';
    fill0[filllen] = '\0';
    char *tmp = new char[l + filllen + 1];

    // 使用循环逐字符拷贝 fill0 到 tmp
    int index = 0;
    for (int i = 0; fill0[i] != '\0'; ++i)
    {
        tmp[index++] = fill0[i];
    }

    // 使用循环逐字符拷贝 num_in_bit 到 tmp
    for (int i = 0; num_in_bit[i] != '\0'; ++i)
    {
        tmp[index++] = num_in_bit[i];
    }
    tmp[index] = '\0'; // 在末尾添加字符串结束符
    num_in_bit = tmp;

    int *blockList = new int[cnt];
    int left = (int) strlen(num_in_bit) - b;
    int right = (int) strlen(num_in_bit);
    for (int i = 0; i < cnt; i++)
    {
        char *blockStr = new char[b + 1];
        for (int j = 0; j < b; j++) blockStr[j] = num_in_bit[left + j];
        blockStr[b] = '\0';
        blockList[i] = (int) strtol(blockStr, nullptr, 2);
        left -= b;
        right -= b;
    }
    return blockList;
}

void compute(mpz_t dst, mpz_t x, mpz_t **table, mpz_t nSquare)
{
    int *blocks = convertIntoBlock(x);
    mpz_t tmp;
    mpz_inits(tmp, dst, NULL);
    mpz_set_si(tmp, 1);
    for (int i = 0; i < cnt; i++)
    {
        mpz_mul(tmp, tmp, table[i][blocks[i]]);
        mpz_mod(tmp, tmp, nSquare);
    }
    mpz_set(dst, tmp);
    mpz_clear(tmp);
}


void Encrypt(mpz_t dst, mpz_t m)
{
    if (mpz_cmp_ui(m, 0) < 0)
    {
        mpz_add(m, m, N);
    }//如果m是负数
    mpz_t r;
    mpz_init(r);
    mpz_urandomb(r, gmp_rand, L_k);
    if (table != nullptr)
    {
        // printf("[TEE] 加密过程-预计算分支!\n");
        mpz_t mask;
        mpz_init(mask);
        compute(mask, r, table, N_Square);//(h^N mod N^2)^r
        // fastPower(mask, pk.h_N, r, N_Square);
        // mpz_powm(mask, pk.h_N, r, N_Square);//正常操作
        mpz_mod(mask, mask, N_Square);//(h^ mod N^2)^r mod N^2

        mpz_mul(dst, m, N);//m*N
        mpz_add_ui(dst, dst, 1);//1+m*N
        mpz_mod(dst, dst, N_Square);//1+m*N mod N^2
        mpz_mul(dst, dst, mask);
        mpz_mod(dst, dst, N_Square);

        mpz_clears(mask, NULL);
    } else
    {
        // printf("[TEE] 加密过程-正常计算分支!\n");
        mpz_t tmp;
        mpz_init(tmp);
        mpz_powm(tmp, h_N, r, N_Square);//正常操作，mpz库自带快速幂，跟下一行的操作是一样的
        // fastPower(tmp, h_N, r, N_Square);

        mpz_mul(dst, m, N);
        mpz_add_ui(dst, dst, 1);
        mpz_mod(dst, dst, N_Square);
        mpz_mul(dst, dst, tmp);
        mpz_mod(dst, dst, N_Square);
        mpz_clears(tmp, NULL);
    }
    mpz_clears(r, NULL);
}

void Enclave_Init(unsigned char **seal_partial_key, size_t *partial_key_size,
                  unsigned char **seal_n, size_t *n_size,
                  unsigned char **seal_h_N, size_t *h_N_size, int temp_L_k)
{

    mpz_inits(partial_key, N, N_Square, half_N, h_N, NULL);
    deserialize_mpz(partial_key, *seal_partial_key, *partial_key_size);
    deserialize_mpz(N, *seal_n, *n_size);
    deserialize_mpz(h_N, *seal_h_N, *h_N_size); //真正花时间的地方是计算h_N

    uint32_t random_value = 0;
    sgx_read_rand((uint8_t * ) & random_value, sizeof(random_value));
    gmp_randinit_default(gmp_rand);
    gmp_randseed_ui(gmp_rand, (unsigned long int) random_value);

    mpz_mul(N_Square, N, N);
    mpz_set(half_N, N);
    mpz_div_ui(half_N, half_N, 2);// half_N = n / 2

    L_k = temp_L_k;
    table = constructTable(h_N, N_Square, L_k);

    mpz_t zero, one;
    mpz_inits(ezero, eone, zero, one, NULL);
    mpz_set_si(zero, 0);
    mpz_set_si(one, 1);
    Encrypt(ezero, zero);
    Encrypt(eone, one);
    mpz_clears(zero, one, NULL);
    printf("Enclave 初始化成功\n");
}

/**
 * 部分解密
 * @param seal_pc
 * @param pc_size
 * @param seal_c
 * @param c_size
 */
void Enclave_PDec(unsigned char **seal_pc, size_t *pc_size,
                  unsigned char **seal_c, size_t *c_size)
{
    mpz_t pc, c;
    mpz_inits(pc, c, NULL);
    deserialize_mpz(c, *seal_c, *c_size);

    PDec(pc, c);
    serialize_mpz(pc, seal_pc, pc_size);
    mpz_clears(pc, c, NULL);
}

void Enclave_Trust_Mask(unsigned char **seal_res, size_t *res_size,
                        unsigned char **seal_Mask_Pdec_1, size_t *Mask_Pdec_1_size,
                        unsigned char **seal_Mask, size_t *Mask_size,
                        size_t *currentBatchSize, size_t *multiExponent)
{

    mpz_t res, Mask, Mask_Pdec_1;
    mpz_inits(res, Mask, Mask_Pdec_1, NULL);
    deserialize_mpz(Mask, *seal_Mask, *Mask_size);
    deserialize_mpz(Mask_Pdec_1, *seal_Mask_Pdec_1, *Mask_Pdec_1_size);

    // Step 2, CSP computes
    mpz_t Mask_Pdec_2, Mask_add_r;
    mpz_inits(Mask_Pdec_2, Mask_add_r, NULL);
    // PDec + TDec
    PDec(Mask_Pdec_2, Mask);
    TDec(Mask_add_r, Mask_Pdec_1, Mask_Pdec_2);
    if (mpz_cmp(Mask_add_r, half_N) > 0)
    {
        mpz_sub(Mask_add_r, Mask_add_r, N);
    }
    mpz_div_ui(res, Mask_add_r, *currentBatchSize);
    mpz_div_ui(res, res, *multiExponent);
    Encrypt(res, res);
    serialize_mpz(res, seal_res, res_size);

    /* send [R] to cp */
    mpz_clears(res, Mask, Mask_Pdec_1, NULL);
    mpz_clears(Mask_Pdec_2, Mask_add_r, NULL);
}


/**
 * 密文乘
 * @param seal_res
 * @param res_size
 * @param seal_X
 * @param X_size
 * @param seal_X_PDec_1
 * @param X_PDec_1_size
 * @param seal_Y
 * @param Y_size
 * @param seal_em_2
 * @param em_2_size
 */
void Enclave_Trust_FMUL(unsigned char **seal_res, size_t *res_size,
                        unsigned char **seal_X, size_t *X_size,
                        unsigned char **seal_X_PDec_1, size_t *X_PDec_1_size,
                        unsigned char **seal_Y, size_t *Y_size,
                        unsigned char **seal_em_2, size_t *em_2_size)
{

    mpz_t res, X, X_PDec_1, Y, em_2;
    mpz_inits(res, X, X_PDec_1, Y, em_2, NULL);
    deserialize_mpz(X, *seal_X, *X_size);
    deserialize_mpz(X_PDec_1, *seal_X_PDec_1, *X_PDec_1_size);
    deserialize_mpz(Y, *seal_Y, *Y_size);
    deserialize_mpz(em_2, *seal_em_2, *em_2_size);

    // Step 2, CSP computes
    mpz_t X_PDec_2, m_1_add_r;
    mpz_inits(X_PDec_2, m_1_add_r, NULL);
    // PDec + TDec
    PDec(X_PDec_2, X);
    TDec(m_1_add_r, X_PDec_1, X_PDec_2);

    scl_mul(res, em_2, m_1_add_r);
    add(res, res, Y);
    add(res, res, ezero); //refresh
    serialize_mpz(res, seal_res, res_size);
    /* send [R] to cp */
    mpz_clears(res, X, X_PDec_1, Y, em_2, NULL);
    mpz_clears(X_PDec_2, m_1_add_r, NULL);

}

/**
 * 密文比较
 * @param seal_res
 * @param res_size
 * @param seal_ed
 * @param ed_size
 * @param seal_ed_1
 * @param ed_1_size
 * @param seal_epi
 * @param epi_size
 */
void Enclave_Trust_FCMP(unsigned char **seal_res, size_t *res_size,
                        unsigned char **seal_ed, size_t *ed_size,
                        unsigned char **seal_ed_1, size_t *ed_1_size,
                        unsigned char **seal_epi, size_t *epi_size)
{

    mpz_t res, ed, ed_1, epi;
    mpz_inits(res, ed, ed_1, epi, NULL);
    deserialize_mpz(ed, *seal_ed, *ed_size);
    deserialize_mpz(ed_1, *seal_ed_1, *ed_1_size);
    deserialize_mpz(epi, *seal_epi, *epi_size);


    // Step 2
    mpz_t ed_2, d, miu0, emiu0, one_sub_2miu0, half_N;
    mpz_inits(ed_2, d, miu0, emiu0, one_sub_2miu0, half_N, NULL);
    mpz_div_ui(half_N, N, 2);
    PDec(ed_2, ed);
    TDec(d, ed_1, ed_2);
    mpz_mod(d, d, N);
    if (mpz_cmp(d, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
        add(emiu0, ezero, ezero);   //refresh
    } else
    {
        mpz_set_si(miu0, 1);
        add(emiu0, eone, ezero);   //refresh
    }
    mpz_mul_si(one_sub_2miu0, miu0, 2);
    mpz_neg(one_sub_2miu0, one_sub_2miu0);
    mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);
    scl_mul(res, epi, one_sub_2miu0);
    add(res, emiu0, res);

    serialize_mpz(res, seal_res, res_size);
    /* send [R] to cp */
    mpz_clears(res, ed, ed_1, epi, NULL);
    mpz_clears(ed_2, d, miu0, emiu0, one_sub_2miu0, half_N, NULL);
}

/**
 * 密文相对比较
 * @param seal_res
 * @param res_size
 * @param seal_ed_1
 * @param ed_1_size
 * @param seal_ed_1_P1
 * @param ed_1_P1_size
 * @param seal_epi_1
 * @param epi_1_size
 * @param seal_ed_2
 * @param ed_2_size
 * @param seal_ed_2_P1
 * @param ed_2_P1_size
 * @param seal_epi_2
 * @param epi_2_size
 */
void Enclave_Trust_FEQL(unsigned char **seal_res, size_t *res_size,
                        unsigned char **seal_ed_1, size_t *ed_1_size,
                        unsigned char **seal_ed_1_P1, size_t *ed_1_P1_size,
                        unsigned char **seal_epi_1, size_t *epi_1_size,
                        unsigned char **seal_ed_2, size_t *ed_2_size,
                        unsigned char **seal_ed_2_P1, size_t *ed_2_P1_size,
                        unsigned char **seal_epi_2, size_t *epi_2_size)
{
    mpz_t res, ed_1, ed_1_P1, epi_1, ed_2, ed_2_P1, epi_2;
    mpz_inits(res, ed_1, ed_1_P1, epi_1, ed_2, ed_2_P1, epi_2, NULL);
    deserialize_mpz(ed_1, *seal_ed_1, *ed_1_size);
    deserialize_mpz(ed_1_P1, *seal_ed_1_P1, *ed_1_P1_size);
    deserialize_mpz(epi_1, *seal_epi_1, *epi_1_size);
    deserialize_mpz(ed_2, *seal_ed_2, *ed_2_size);
    deserialize_mpz(ed_2_P1, *seal_ed_2_P1, *ed_2_P1_size);
    deserialize_mpz(epi_2, *seal_epi_2, *epi_2_size);



    // Step 2
    mpz_t ed_1_P2, ed_2_P2, d1, d2, miu0, emiu0, miu0_prime, emiu0_prime, one_sub_2miu0, one_sub_2miu0_prime, half_N;
    mpz_inits(ed_1_P2, ed_2_P2, d1, d2, miu0, emiu0, miu0_prime, emiu0_prime, one_sub_2miu0, one_sub_2miu0_prime, half_N, NULL);
    mpz_div_ui(half_N, N, 2);

    PDec(ed_1_P2, ed_1);
    TDec(d1, ed_1_P1, ed_1_P2);
    PDec(ed_2_P2, ed_2);
    TDec(d2, ed_2_P1, ed_2_P2);
    mpz_mod(d1, d1, N);
    if (mpz_cmp(d1, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
        add(emiu0, ezero, ezero);   //refresh
    } else
    {
        mpz_set_si(miu0, 1);
        add(emiu0, eone, ezero);   //refresh
    }
    mpz_mul_si(one_sub_2miu0, miu0, 2);
    mpz_neg(one_sub_2miu0, one_sub_2miu0);
    mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);

    mpz_mod(d2, d2, N);
    if (mpz_cmp(d2, half_N) > 0)
    {
        mpz_set_si(miu0_prime, 0);
        add(emiu0_prime, ezero, ezero);   //refresh
    } else
    {
        mpz_set_si(miu0_prime, 1);
        add(emiu0_prime, eone, ezero);   //refresh
    }
    mpz_mul_si(one_sub_2miu0_prime, miu0_prime, 2);
    mpz_neg(one_sub_2miu0_prime, one_sub_2miu0_prime);
    mpz_add_ui(one_sub_2miu0_prime, one_sub_2miu0_prime, 1);

    mpz_t left, right;
    mpz_inits(left, right, NULL);
    scl_mul(left, epi_1, one_sub_2miu0);
    add(left, emiu0, left);
    scl_mul(right, epi_2, one_sub_2miu0_prime);
    add(right, emiu0_prime, right);
    add(res, left, right);
    serialize_mpz(res, seal_res, res_size);
    /* send [R] to cp */
    mpz_clears(res, ed_1, ed_1_P1, epi_1, ed_2, ed_2_P1, epi_2, NULL);
    mpz_clears(ed_1_P2, ed_2_P2, d1, d2, miu0, emiu0, miu0_prime, emiu0_prime, one_sub_2miu0, one_sub_2miu0_prime, half_N, NULL);
    mpz_clears(left, right, NULL);
}

/**
 * 密文绝对值
 * @param seal_res
 * @param res_size
 * @param seal_ed
 * @param ed_size
 * @param seal_ed_1
 * @param ed_1_size
 * @param seal_em_1pi
 * @param em_1pi_size
 * @param seal_em_1
 * @param em_1_size
 */
void Enclave_Trust_FABS(unsigned char **seal_res, size_t *res_size,
                        unsigned char **seal_ed, size_t *ed_size,
                        unsigned char **seal_ed_1, size_t *ed_1_size,
                        unsigned char **seal_em_1pi, size_t *em_1pi_size,
                        unsigned char **seal_em_1, size_t *em_1_size)
{

    mpz_t res, ed, ed_1, em_1pi, em_1;
    mpz_inits(res, ed, ed_1, em_1pi, em_1, NULL);
    deserialize_mpz(ed, *seal_ed, *ed_size);
    deserialize_mpz(ed_1, *seal_ed_1, *ed_1_size);
    deserialize_mpz(em_1pi, *seal_em_1pi, *em_1pi_size);
    deserialize_mpz(em_1, *seal_em_1, *em_1_size);

    // Step 2
    mpz_t ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N;
    mpz_inits(ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N, NULL);
    mpz_div_ui(half_N, N, 2);
    PDec(ed_2, ed);
    TDec(d, ed_1, ed_2);
    mpz_mod(d, d, N);
    if (mpz_cmp(d, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
    } else
    {
        mpz_set_si(miu0, 1);
    }
    mpz_mul_si(one_sub_2miu0, miu0, 2);
    mpz_neg(one_sub_2miu0, one_sub_2miu0);
    mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);
    mpz_mul_ui(neg_two_sub_4miu0, one_sub_2miu0, 2);
    mpz_neg(neg_two_sub_4miu0, neg_two_sub_4miu0);


    mpz_t left, right;
    mpz_inits(left, right, NULL);
    scl_mul(left, em_1, one_sub_2miu0);
    scl_mul(right, em_1pi, neg_two_sub_4miu0);
    add(res, left, right);
    add(res, res, ezero);   //refresh

    serialize_mpz(res, seal_res, res_size);
    /* send [R] to cp */
    mpz_clears(res, ed, ed_1, em_1pi, em_1, NULL);
    mpz_clears(ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N, NULL);
    mpz_clears(left, right, NULL);
}

/**
 * 密文运算符version1
 * @param seal_res
 * @param res_size
 * @param seal_em_3_sub_m_2
 * @param em_3_sub_m_2_size
 * @param seal_em_2
 * @param em_2_size
 * @param seal_ed_1
 * @param ed_1_size
 * @param seal_ed_1_P1
 * @param ed_1_P1_size
 * @param seal_epi_1_m_3_sub_m_2
 * @param epi_1_m_3_sub_m_2_size
 * @param seal_ed_2
 * @param ed_2_size
 * @param seal_ed_2_P1
 * @param ed_2_P1_size
 * @param seal_epi_2_m_3_sub_m_2
 * @param epi_2_m_3_sub_m_2_size
 */
void Enclave_Trust_FTRN_v1(unsigned char **seal_res, size_t *res_size,
                           unsigned char **seal_em_3_sub_m_2, size_t *em_3_sub_m_2_size,
                           unsigned char **seal_em_2, size_t *em_2_size,
                           unsigned char **seal_ed_1, size_t *ed_1_size,
                           unsigned char **seal_ed_1_P1, size_t *ed_1_P1_size,
                           unsigned char **seal_epi_1_m_3_sub_m_2, size_t *epi_1_m_3_sub_m_2_size,
                           unsigned char **seal_ed_2, size_t *ed_2_size,
                           unsigned char **seal_ed_2_P1, size_t *ed_2_P1_size,
                           unsigned char **seal_epi_2_m_3_sub_m_2, size_t *epi_2_m_3_sub_m_2_size)
{
    mpz_t res, em_3_sub_m_2, em_2, ed_1, ed_1_P1, epi_1_m_3_sub_m_2, ed_2, ed_2_P1, epi_2_m_3_sub_m_2;
    mpz_inits(res, em_3_sub_m_2, em_2, ed_1, ed_1_P1, epi_1_m_3_sub_m_2, ed_2, ed_2_P1, epi_2_m_3_sub_m_2, NULL);
    deserialize_mpz(em_3_sub_m_2, *seal_em_3_sub_m_2, *em_3_sub_m_2_size);
    deserialize_mpz(em_2, *seal_em_2, *em_2_size);
    deserialize_mpz(ed_1, *seal_ed_1, *ed_1_size);
    deserialize_mpz(ed_1_P1, *seal_ed_1_P1, *ed_1_P1_size);
    deserialize_mpz(epi_1_m_3_sub_m_2, *seal_epi_1_m_3_sub_m_2, *epi_1_m_3_sub_m_2_size);
    deserialize_mpz(ed_2, *seal_ed_2, *ed_2_size);
    deserialize_mpz(ed_2_P1, *seal_ed_2_P1, *ed_2_P1_size);
    deserialize_mpz(epi_2_m_3_sub_m_2, *seal_epi_2_m_3_sub_m_2, *epi_2_m_3_sub_m_2_size);


    // Step 2
    mpz_t ed_1_P2, ed_2_P2, d1, d2, miu0, miu0_prime, one_sub_2miu0, one_sub_2miu0_prime, miu0_add_miu0_prime, half_N;
    mpz_inits(ed_1_P2, ed_2_P2, d1, d2, miu0, miu0_prime, one_sub_2miu0, one_sub_2miu0_prime, miu0_add_miu0_prime, half_N, NULL);
    mpz_div_ui(half_N, N, 2);

    PDec(ed_1_P2, ed_1);
    TDec(d1, ed_1_P1, ed_1_P2);
    PDec(ed_2_P2, ed_2);
    TDec(d2, ed_2_P1, ed_2_P2);
    if (mpz_cmp(d1, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
    } else if (mpz_cmp(d1, half_N) < 0)
    {
        mpz_set_si(miu0, 1);
    } else
    {
        printf(" ERROR !\n");
        throw;
    }
    mpz_mul_si(one_sub_2miu0, miu0, 2);
    mpz_neg(one_sub_2miu0, one_sub_2miu0);
    mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);
    if (mpz_cmp(d2, half_N) > 0)
    {
        mpz_set_si(miu0_prime, 0);
    } else if (mpz_cmp(d2, half_N) < 0)
    {
        mpz_set_si(miu0_prime, 1);
    } else
    {
        printf(" ERROR !\n");
        throw;
    }
    mpz_mul_si(one_sub_2miu0_prime, miu0_prime, 2);
    mpz_neg(one_sub_2miu0_prime, one_sub_2miu0_prime);
    mpz_add_ui(one_sub_2miu0_prime, one_sub_2miu0_prime, 1);
    mpz_add(miu0_add_miu0_prime, miu0, miu0_prime);

    // [R] compute
    mpz_t left, mid, right;
    mpz_inits(left, mid, right, NULL);
    scl_mul(left, em_3_sub_m_2, miu0_add_miu0_prime); // left
    add(left, em_2, left); // left
    scl_mul(mid, epi_1_m_3_sub_m_2, one_sub_2miu0); // mid
    scl_mul(right, epi_2_m_3_sub_m_2, one_sub_2miu0_prime); // right

    add(res, left, mid);
    add(res, res, right);
    add(res, res, ezero);   //refresh

    serialize_mpz(res, seal_res, res_size);
    /* send [R] to cp */
    mpz_clears(res, em_3_sub_m_2, em_2, ed_1, ed_1_P1, epi_1_m_3_sub_m_2, ed_2, ed_2_P1, epi_2_m_3_sub_m_2, NULL);
    mpz_clears(ed_1_P2, ed_2_P2, d1, d2, miu0, miu0_prime, one_sub_2miu0, one_sub_2miu0_prime, miu0_add_miu0_prime, half_N, NULL);
    mpz_clears(left, mid, right, NULL);
}

/**
 * 密文运算符version2
 * @param seal_res
 * @param res_size
 * @param seal_em_3
 * @param em_3_size
 * @param seal_em_2
 * @param em_2_size
 * @param seal_ed
 * @param ed_size
 * @param seal_ed_1
 * @param ed_1_size
 * @param seal_epi_m_3_sub_m_2
 * @param epi_m_3_sub_m_2_size
 */
void Enclave_Trust_FTRN_v2(unsigned char **seal_res, size_t *res_size,
                           unsigned char **seal_em_3, size_t *em_3_size,
                           unsigned char **seal_em_2, size_t *em_2_size,
                           unsigned char **seal_ed, size_t *ed_size,
                           unsigned char **seal_ed_1, size_t *ed_1_size,
                           unsigned char **seal_epi_m_3_sub_m_2, size_t *epi_m_3_sub_m_2_size)
{
    mpz_t res, em_3, em_2, ed, ed_1, epi_m_3_sub_m_2;
    mpz_inits(res, em_3, em_2, ed, ed_1, epi_m_3_sub_m_2, NULL);
    deserialize_mpz(em_3, *seal_em_3, *em_3_size);
    deserialize_mpz(em_2, *seal_em_2, *em_2_size);
    deserialize_mpz(ed, *seal_ed, *ed_size);
    deserialize_mpz(ed_1, *seal_ed_1, *ed_1_size);
    deserialize_mpz(epi_m_3_sub_m_2, *seal_epi_m_3_sub_m_2, *epi_m_3_sub_m_2_size);


    // Step 2
    mpz_t ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N;
    mpz_inits(ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N, NULL);
    mpz_div_ui(half_N, N, 2);

    PDec(ed_2, ed);
    TDec(d, ed_1, ed_2);
    if (mpz_cmp(d, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
    } else
    {
        mpz_set_si(miu0, 1);
    }

    mpz_neg(one_sub_miu0, miu0);
    mpz_add_ui(one_sub_miu0, one_sub_miu0, 1);

    mpz_mul_si(one_sub_2miu0, miu0, 2);
    mpz_neg(one_sub_2miu0, one_sub_2miu0);
    mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);

    // [R] compute
    mpz_t left, mid, right;
    mpz_inits(left, mid, right, NULL);
    scl_mul(left, em_2, one_sub_miu0); // left
    scl_mul(mid, em_3, miu0); // mid
    scl_mul(right, epi_m_3_sub_m_2, one_sub_2miu0); // right

    add(res, left, mid);
    add(res, res, right);
    add(res, res, ezero);   //refresh

    serialize_mpz(res, seal_res, res_size);
    /* send [R] to cp */
    mpz_clears(res, em_3, em_2, ed, ed_1, epi_m_3_sub_m_2, NULL);
    mpz_clears(ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N, NULL);
    mpz_clears(left, mid, right, NULL);
}

void Enclave_Release()
{

    mpz_clears(ezero, eone, NULL);
    mpz_clears(partial_key, N, N_Square, half_N, h_N, NULL);
    for (int i = 0; i < rows; ++i)
    {
        for (int j = 0; j < columns; ++j)
        {
            mpz_clear(table[i][j]);
        }
        free(table[i]); // 释放第 i 行的内存
    }
    free(table); // 最后释放指向表格的指针
    printf("Enclave Release is ok!\n");

}



