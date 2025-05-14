#ifndef OFFLINE_H
#define OFFLINE_H

#pragma once

namespace OFFLINESPACE {
    class Offline_val {
    public:
        Offline_val();

        ~Offline_val();

        mpz_t r1, r2, er1, er2, neg_er1r2, r3, r4, er3r4, er4, e0, e1;

        mpz_t trust_r, trust_r1, trust_r2, trust_r1_prime, trust_r2_prime;
        mpz_t e_trust_r, e_trust_r1, e_trust_r2, e_trust_r1_prime, e_trust_r2_prime;
        mpz_t e_trust_r1_add_r2, e_trust_r1_prime_add_r2_prime;
        mpz_t e_trust_2r1_prime_add_r2_prime, e_trust_r2_prime_sub_r1_prime;

    };

    Offline_val::Offline_val() {}

    Offline_val::~Offline_val() {}

    //tupleS0
    //e0, e1
    class Pre_Compute {
    private:
        static int blockSize;
    public:
        Pre_Compute();

        ~Pre_Compute();

        static mpz_t **constructTable(mpz_t g, mpz_t nSquare, int l);

        static int *convertIntoBlock(mpz_t num, int &cnt);

        static void compute(mpz_t dst, mpz_t x, mpz_t **table, mpz_t nSquare);

        static void fastPower(mpz_t result, const mpz_t base, const mpz_t exponent, const mpz_t modulus);
    };

    Pre_Compute::Pre_Compute() {}

    Pre_Compute::~Pre_Compute() {}

    int Pre_Compute::blockSize = 5;

    mpz_t **Pre_Compute::constructTable(mpz_t g, mpz_t nSquare, int l) {
        mpz_t b, base;
        mpz_inits(b, base, NULL);
        mpz_set_si(b, blockSize);
        int lDivB = (l + blockSize - 1) / blockSize;
        mpz_set(base, g);
        auto **table = new mpz_t *[lDivB];
        for (int i = 0; i < lDivB; i++) {
            int size = 1 << blockSize;
            table[i] = new mpz_t[size];
            mpz_t tmp;
            mpz_init(tmp);
            for (int j = 0; j < size; j++) {
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
        return table;
    }

    int *Pre_Compute::convertIntoBlock(mpz_t num, int &cnt) {
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
        strcpy(tmp, fill0);
        strcat(tmp, num_in_bit);
        delete[] fill0;
        delete[] num_in_bit;
        num_in_bit = tmp;

        int *blockList = new int[cnt];
        int left = (int) strlen(num_in_bit) - b;
        int right = (int) strlen(num_in_bit);
        for (int i = 0; i < cnt; i++) {
            char *blockStr = new char[b + 1];
            for (int j = 0; j < b; j++) blockStr[j] = num_in_bit[left + j];
            blockStr[b] = '\0';
            blockList[i] = (int) strtol(blockStr, nullptr, 2);
            delete[] blockStr;
            left -= b;
            right -= b;
        }
        delete[] num_in_bit;
        return blockList;
    }

    void Pre_Compute::compute(mpz_t dst, mpz_t x, mpz_t **table, mpz_t nSquare) {
        int cnt;
        int *blocks = convertIntoBlock(x, cnt);
        mpz_t tmp;
        mpz_inits(tmp, dst, NULL);
        mpz_set_si(tmp, 1);
        for (int i = 0; i < cnt; i++) {
            mpz_mul(tmp, tmp, table[i][blocks[i]]);
            mpz_mod(tmp, tmp, nSquare);
        }
        mpz_set(dst, tmp);
        mpz_clear(tmp);
    }

    void Pre_Compute::fastPower(mpz_t result, const mpz_t base, const mpz_t exponent, const mpz_t modulus) {
        mpz_set_ui(result, 1);  // 初始化结果为1
        mpz_t base_copy;
        mpz_init_set(base_copy, base);  // 复制底数用于计算

        mpz_t exp;
        mpz_init_set(exp, exponent);  // 用于迭代计算指数

        while (mpz_cmp_ui(exp, 0) > 0) {
            if (mpz_odd_p(exp)) {
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

}

#endif //OFFLINE_H
