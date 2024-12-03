/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-01
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#ifndef POCF_H
#define POCF_H
#pragma once

#include "pcpd.h"
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>

using namespace std;
using namespace PHESPACE;

/**
 * POCF里面的乘法协议
 * @param result
 * @param e_x
 * @param e_y
 * @param cp
 * @param csp
 */
void SMUL(mpz_t result, const mpz_t &e_x, const mpz_t &e_y, const Paillier_Third_Party &cp,
          const Paillier_Third_Party &csp) {

    mpz_t N, N_square, N_half, flag, temp, e_r_x, e_r_y, X, Y, X1, Y1, X2, Y2, h, H, neg_one, S_1, S_2, S_3, e_x_y, r_x, r_y;
    mpz_inits(N, N_square, flag, N_half,
              temp, e_r_x, e_r_y, X, Y, X1, Y1, X2, Y2, h, H, neg_one, S_1, S_2, S_3, e_x_y, r_x, r_y,
              NULL);

    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);
    mpz_set(N_half, cp.public_key.N_half);

    // step 1 (CP)
    // 选取两个随机数 r_x 和 r_y
    while (true) {
        mpz_urandomm(r_x, randstate, N);
        mpz_gcd(flag, r_x, N);
        if (mpz_cmp_ui(flag, 1) == 0) {
            break;
        }
    }
    while (true) {
        mpz_urandomm(r_y, randstate, N);
        mpz_gcd(flag, r_y, N);
        if (mpz_cmp_ui(flag, 1) == 0) {
            break;
        }
    }

    // 加密随机数 r_x 和 r_y
    Enc(e_r_x, cp.public_key, r_x);
    Enc(e_r_y, cp.public_key, r_y);

    // 计算 X = e_x * e_r_x mod N_square 和 Y = e_y * e_r_y mod N_square
    mpz_mul(X, e_x, e_r_x);
    mpz_mod(X, X, N_square);
    mpz_mul(Y, e_y, e_r_y);
    mpz_mod(Y, Y, N_square);
    // 部分解密 X 和 Y
    PDec(X1, cp.private_key, X);
    PDec(Y1, cp.private_key, Y);
    // 发X,Y,X1,Y1给CSP

    // step 2 (CSP)
    // 计算 h = TDec(X1, PDec(csp. private_key, X), N) * TDec(Y1, PDec(csp. private_key, Y), N)
    PDec(X2, csp.private_key, X);
    PDec(Y2, csp.private_key, Y);
    mpz_t tdec_temp_1, tdec_temp_2;
    mpz_inits(tdec_temp_1, tdec_temp_2, NULL);
    TDec(tdec_temp_1, X1, X2, N_square, N, N_half);
    TDec(tdec_temp_2, Y1, Y2, N_square, N, N_half);
    mpz_mul(h, tdec_temp_1, tdec_temp_2);
    // 加密 h
    Enc(H, csp.public_key, h);
    mpz_clears(tdec_temp_1, tdec_temp_2, NULL);
    // 发H给CP

    // step 3 (CP)
    // 计算 S_1, S_2, S_3
    mpz_set_si(neg_one, -1);
    mpz_mul(temp, r_x, r_y);
    mpz_mod(temp, temp, N);
    Enc(S_1, cp.public_key, temp);
    mpz_powm(S_1, S_1, neg_one, N_square);

    mpz_neg(temp, r_y);
    mpz_powm(S_2, e_x, temp, N_square);

    mpz_neg(temp, r_x);
    mpz_powm(S_3, e_y, temp, N_square);

    // 计算 e_x_y
    mpz_mul(temp, H, S_1);
    mpz_mod(temp, temp, N_square);
    mpz_mul(result, temp, S_2);
    mpz_mod(result, result, N_square);
    mpz_mul(result, result, S_3);
    mpz_mod(result, result, N_square);

    // 清理内存
    mpz_clears(N, N_square, N_half,
               flag, temp, e_r_x, e_r_y, X, Y, X1, Y1, X2, Y2, h, H, neg_one, S_1, S_2, S_3, r_x, r_y,
               NULL);
}


/**
 * POCF的比较，与SOCI的输出有区别，在这里是返回一个加密的u，如果x>=y,则u等于0；否则为u=1
 * @param e_u
 * @param e_x
 * @param e_y
 * @param cp
 * @param csp
 */
void SCMP(mpz_t e_u, const mpz_t &e_x, const mpz_t &e_y,
          const Paillier_Third_Party &cp, const Paillier_Third_Party &csp) {
    mpz_t one, eone, zero, neg_one, two, temp, N, N_square, N_half, e_x_1, e_y_1, new_modulus, e_l, K, l, u_another;
    mpz_inits(one, eone, zero, temp, neg_one, two, N, N_square, N_half, e_x_1, e_y_1, new_modulus, e_l, K, l, u_another,
              NULL);

    mpz_set_si(zero, 0);
    mpz_set_si(one, 1);
    mpz_set_si(two, 2);
    mpz_set_si(neg_one, -1);
    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);
    mpz_set(N_half, cp.public_key.N_half);


    // step1 cp
    // e_x_1 = 2x+1, e_y_1=2y
    Enc(eone, cp.public_key, one);
    mpz_powm(e_x_1, e_x, two, N_square);
    mpz_mul(e_x_1, e_x_1, eone);
    mpz_mod(e_x_1, e_x_1, N_square);
    mpz_powm(e_y_1, e_y, two, N_square);

    srand(time(NULL));
    int s = rand() % 2; // 生成0或1的随机数
    int N_bits = mpz_sizeinbase(N, 2);
    mpz_powm_ui(new_modulus, two, (N_bits / 4 - 1), N_square);
    mpz_t r;
    mpz_init(r);
    mpz_urandomm(r, randstate, new_modulus);
    if (mpz_cmp_si(r, 0) == 0) {
        mpz_add(r, r, one);
    }

    if (s == 1) {
        mpz_powm(temp, e_y_1, neg_one, N_square); // -y1
        mpz_mul(temp, e_x_1, temp); // x1-y1
        mpz_mod(temp, temp, N_square);
        mpz_powm(e_l, temp, r, N_square); // r(x1-y1)
    } else {
        mpz_powm(temp, e_x_1, neg_one, N_square); // -x1
        mpz_mul(temp, e_y_1, temp); // y1-x1
        mpz_mod(temp, temp, N_square);
        mpz_powm(e_l, temp, r, N_square); // r(y1-x1)
    }
    PDec(K, cp.private_key, e_l);
    // send K and e_l to CSP

    // step2 CSP
    PDec(temp, csp.private_key, e_l);
    TDec(l, K, temp, N_square, N, N_half);
    if (mpz_cmp_si(l, 0) < 0) {
        mpz_add(l, l, cp.public_key.N);
    }
    if (mpz_sizeinbase(l, 2) > (mpz_sizeinbase(N, 2) / 2)) {
        Enc(u_another, csp.public_key, one);
    } else {
        Enc(u_another, csp.public_key, zero);
    }
    // send u_another to CP

    // step3 CP
    if (s == 1) {
        CR(u_another, N, N_square);
        mpz_set(e_u, u_another);
    } else if (s == 0) {
        mpz_set(temp, eone);
        mpz_powm(u_another, u_another, neg_one, N_square);
        mpz_mul(e_u, temp, u_another);
        mpz_mod(e_u, e_u, N_square);
    }
    mpz_clears(one, zero, neg_one, two, temp, N, N_square, N_half, e_x_1, e_y_1, new_modulus, e_l, K, l, u_another,
               NULL);
}

/**
 * POCF里面的异或运算，限定范围 x, y ∈{0, 1}
 * @param result
 * @param e_x
 * @param e_y
 * @param cp
 * @param csp
 */
void SXOR(mpz_t result, const mpz_t &e_x, const mpz_t &e_y,
          const Paillier_Third_Party &cp, const Paillier_Third_Party &csp) {
    mpz_t one, neg_one, N, N_square, temp, e_1, f_1, f_2;
    mpz_inits(one, neg_one, N, N_square, temp, e_1, f_1, f_2, NULL);

    mpz_set_si(one, 1);
    mpz_set_si(neg_one, -1);
    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);

    // 计算 e_1
    Enc(e_1, cp.public_key, one);

    // 计算 f_1
    mpz_t inv_e_x;
    mpz_init(inv_e_x);
    mpz_powm(inv_e_x, e_x, neg_one, N_square);
    mpz_mul(inv_e_x, inv_e_x, e_1);
    mpz_mod(inv_e_x, inv_e_x, N_square);
    SMUL(f_1, inv_e_x, e_y, cp, csp);
    mpz_clear(inv_e_x);


    // 计算 f_2
    mpz_t inv_e_y;
    mpz_init(inv_e_y);
    mpz_powm(inv_e_y, e_y, neg_one, N_square);
    mpz_mul(inv_e_y, inv_e_y, e_1);
    mpz_mod(inv_e_y, inv_e_y, N_square);
    SMUL(f_2, e_x, inv_e_y, cp, csp);
    mpz_clear(inv_e_y);


    // 计算 f
    mpz_mul(result, f_1, f_2);
    mpz_mod(result, result, N_square);

    mpz_clears(one, neg_one, N, N_square, temp, e_1, f_1, f_2, NULL);
}

/**
 * POCF里面的相等判断
 * @param result
 * @param e_x
 * @param e_y
 * @param cp
 * @param csp
 */
void SEQ(mpz_t result, const mpz_t &e_x, const mpz_t &e_y,
         const Paillier_Third_Party &cp, const Paillier_Third_Party &csp) {
    mpz_t u_1, u_2;
    mpz_inits(u_1, u_2, NULL);

    // 计算 u_1
    SCMP(u_1, e_x, e_y, cp, csp);

    // 计算 u_2
    SCMP(u_2, e_y, e_x, cp, csp);

    // 计算 f
    SXOR(result, u_1, u_2, cp, csp);

    mpz_clears(u_1, u_2, NULL);
}



/**
 * 来源论文：An Efficient and Probabilistic Secure Bit-Decomposition
 * @param e_x_i
 * @param T
 * @param cp
 * @param csp
 */
void Encrypted_LSB(mpz_t e_x_i, const mpz_t &T, const Paillier_Third_Party &cp, const Paillier_Third_Party &csp) {
    mpz_t temp, zero, ezero, one, eone, neg_one, N, N_square, N_half, r, er, Y, Y1, Y2, y, alpha;
    mpz_inits(temp, zero, ezero, one, eone, neg_one, N, N_square, N_half, r, er, Y, Y1, Y2, y, alpha, NULL);


    mpz_set_si(zero, 0);
    mpz_set_si(one, 1);
    mpz_set_si(neg_one, -1);
    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);
    mpz_set(N_half, cp.public_key.N_half);
    Enc(eone, cp.public_key, one);

    // Bob:
    mpz_urandomm(r, randstate, N);
    Enc(er, cp.public_key, r);
    mpz_mul(Y, T, er);
    mpz_mod(Y, Y, N_square);
    PDec(Y1, cp.private_key, Y);

    // Alice:
    PDec(Y2, csp.private_key, Y);
    TDec(y, Y1, Y2, N_square, N, N_half);
    if (mpz_cmp_si(y, 0) < 0) {
        mpz_add(y, y, cp.public_key.N);
    }
    if (mpz_even_p(y)) {
        // 如果 y 是偶数
        Enc(alpha, csp.public_key, zero);
    } else {
        Enc(alpha, csp.public_key, one);
    }

    // Bob:
    if (mpz_even_p(r)) {
        mpz_set(e_x_i, alpha);
    } else {
        mpz_powm(Y1, alpha, neg_one, N_square);
        mpz_mul(e_x_i, eone, Y1);
        mpz_mod(e_x_i, e_x_i, N_square);
    }

    mpz_clears(temp, zero, one, eone, neg_one, N, N_square, N_half, r, er, Y, Y1, Y2, y, alpha, NULL);
}

int SVR(const mpz_t &e_x, const vector<mpz_t *> &enc_bits,
        const Paillier_Third_Party &cp, const Paillier_Third_Party &csp) {
    mpz_t N, N_square, N_half, one, two, temp, neg_one, U, V, r, W, W1, W2, D_w;
    mpz_inits(N, N_square, N_half, one, two, temp, neg_one, U, V, r, W, W1, W2, D_w, NULL);


    mpz_set_si(two, 2);
    mpz_set_si(neg_one, -1);
    mpz_set_si(one, 1);
    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);
    mpz_set(N_half, cp.public_key.N_half);

    // Bob
    mpz_set(U, *enc_bits[0]);
    for (int i = 1; i < enc_bits.size(); ++i) {
        mpz_set_si(temp, i);
        mpz_powm(temp, two, temp, N_square);
        mpz_powm(temp, *enc_bits[i], temp, N_square);
        mpz_mul(U, U, temp);
        mpz_mod(U, U, N_square);
    }

    mpz_powm(V, e_x, neg_one, N_square);
    mpz_mul(V, V, U);
    mpz_mod(V, V, N_square);

    mpz_sub(temp, N, one);
    mpz_urandomm(r, randstate, temp);
    mpz_powm(W, V, r, N_square);
    PDec(W1, cp.private_key, W);

    // Alice
    mpz_t part_result1, part_result2, result;
    mpz_inits(part_result1, part_result2, result, NULL);
    PDec(W2, csp.private_key, W);
    TDec(D_w, W1, W2, N_square, N, N_half);
    if (mpz_cmp_si(D_w, 0) == 0) {
        mpz_set_si(temp, 1);
    } else {
        mpz_set_si(temp, 0);
    }
    Enc(result, csp.public_key, temp);
    PDec(part_result2, csp.private_key, result);

    // Bob
    PDec(part_result1, cp.private_key, result);
    TDec(D_w, part_result1, part_result2, N_square, N, N_half);
    if (mpz_cmp_si(D_w, 0) == 0) {
        return 0;
    } else {
        return 1;
    }
    mpz_clears(N, N_square, N_half, one, two, temp, neg_one, U, V, r, W, W1, W2, D_w, part_result1, part_result2,
               result, NULL);
}


/**
 * secure bit-decomposition(SBD) 来源论文 An Efficient and Probabilistic Secure Bit-Decomposition 将解密更换了部分解密
 * @param enc_bits
 * @param e_x
 * @param m
 * @param cp
 * @param csp
 */
void SBD(vector<mpz_t *> &enc_bits, const mpz_t &e_x, const int m,
         const Paillier_Third_Party &cp, const Paillier_Third_Party &csp) {
    mpz_t two, neg_one, N, N_square, temp, l, T, Z;
    mpz_inits(two, neg_one, N, N_square, temp, l, T, Z, NULL);

    mpz_set_si(two, 2);
    mpz_set_si(neg_one, -1);
    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);

    mpz_invert(l, two, N);

    while (true) {
        mpz_set(T, e_x);
        enc_bits.clear();
        for (int i = 0; i < m; ++i) {
            mpz_t e_x_i;
            mpz_init(e_x_i);
            Encrypted_LSB(e_x_i, T, cp, csp);
            mpz_t *ptr = (mpz_t *) malloc(sizeof(mpz_t));
            mpz_init_set(*ptr, e_x_i);
            enc_bits.push_back(ptr);
            mpz_powm(Z, e_x_i, neg_one, N_square);
            mpz_mul(Z, Z, T);
            mpz_mod(Z, Z, N_square);
            mpz_powm(T, Z, l, N_square);
        }

        int gama = SVR(e_x, enc_bits, cp, csp);
        if (gama == 1) {
            break;
        }
    }
    reverse(enc_bits.begin(), enc_bits.end());  // 高位在前面
}


/**
 * X. Liu, R. H. Deng, K.-K. R. Choo, and J. Weng,
 * “An efficient privacypreserving outsourced calculation toolkit with multiple keys”
 * IEEE Transactions on Information Forensics and Security, vol. 11, no. 11, pp. 2401–2414, 2016
 * @param result_e_s
 * @param result_x_another
 * @param e_x
 * @param cp
 * @param csp
 */
void SSBA(mpz_t result_e_s, mpz_t result_x_another, const mpz_t &e_x,
          const Paillier_Third_Party &cp, const Paillier_Third_Party &csp) {
    mpz_t one, two, neg_one, eone, N, N_square, N_half, new_modulus,
            r, temp, e_l, L, e_u, e_s, e_y, x_another, l;
    mpz_inits(one, two, neg_one, eone, N, N_square, N_half, new_modulus,
              r, temp, e_l, L, e_u, e_s, e_y, x_another, l, NULL);

    mpz_set_si(one, 1);
    mpz_set_si(neg_one, -1);
    mpz_set_si(two, 2);
    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);
    mpz_set(N_half, cp.public_key.N_half);
    Enc(eone, cp.public_key, one);

    // step 1 CP
    int s = rand() % 2; // 生成0或1的随机数
    int N_bits = mpz_sizeinbase(N, 2);
    mpz_powm_ui(new_modulus, two, N_bits / 4, N_square);
    mpz_urandomm(r, randstate, new_modulus);
    if (mpz_cmp_si(r, 0) == 0) {
        mpz_add(r, r, one);
    }
    if (s == 1) {
        mpz_powm(e_l, e_x, two, N_square);
        mpz_mul(e_l, e_l, eone);
        mpz_mod(e_l, e_l, N_square);
        mpz_powm(e_l, e_l, r, N_square);
    } else {
        mpz_powm(e_l, e_x, two, N_square);
        mpz_mul(e_l, e_l, eone);
        mpz_mod(e_l, e_l, N_square);
        mpz_t neg_r;
        mpz_init(neg_r);
        mpz_neg(neg_r, r);
        mpz_powm(e_l, e_l, neg_r, N_square);
    }

    PDec(L, cp.private_key, e_l);

    // 将这个 e_l 和 L 发给 csp

    // step 2 csp
    PDec(temp, csp.private_key, e_l);
    TDec(l, L, temp, N_square, N, N_half);
    if (mpz_cmp_si(l, 0) < 0) {
        mpz_add(l, l, cp.public_key.N);
    }
    int bit_length_l = mpz_sizeinbase(l, 2);
    int threshold = 3 * mpz_sizeinbase(N, 2) / 8;
    mpz_t u;
    mpz_init(u);
    if (bit_length_l < threshold) {
        mpz_set_si(u, 1);
    } else {
        mpz_set_si(u, 0);
    }
    Enc(e_u, csp.public_key, u);
    // 将这个 e_u 发给 cp

    // step 3 cp
    if (s == 1) {
        CR(e_u, N, N_square);
        mpz_set(e_s, e_u);
    } else {
        mpz_powm(e_s, e_u, neg_one, N_square);
        mpz_mul(e_s, eone, e_s);
        mpz_mod(e_s, e_s, N_square);
        CR(e_s, N, N_square);
    }
    mpz_powm(e_y, e_s, two, N_square);
    mpz_powm(temp, eone, neg_one, N_square);
    mpz_mul(e_y, e_y, temp);
    mpz_mod(e_y, e_y, N_square);
    SMUL(x_another, e_x, e_y, cp, csp);

    mpz_set(result_e_s, e_s);
    mpz_set(result_x_another, x_another);

    mpz_clears(one, two, neg_one, eone, N, N_square, N_half, new_modulus,
               r, temp, e_l, L, e_u, e_s, e_y,
               x_another, l, NULL);
}

/**
 * X. Liu, R. H. Deng, K.-K. R. Choo, and J. Weng,
 * “An efficient privacypreserving outsourced calculation toolkit with multiple keys”
 * IEEE Transactions on Information Forensics and Security, vol. 11, no. 11, pp. 2401–2414, 2016
 * @param result_q_xing
 * @param result_r_xing
 * @param e_x
 * @param e_y
 * @param l
 * @param cp
 * @param csp
 * @param pai
 */
void SDIV(mpz_t result_q_xing, mpz_t result_r_xing, const mpz_t &e_x, const mpz_t &e_y, const int l,
          const Paillier_Third_Party &cp, const Paillier_Third_Party &csp, Paillier pai) {
    mpz_t zero, one, neg_one, two, N, N_square, N_half, s_x, x_xing, s_y, y_xing, e_0, e_1, e_r, e_q, e_f, e_1_f, e_fx, e_K_1, e_y_another, e_x_another;
    // 初始化mpz_t变量
    mpz_inits(zero, one, neg_one, two, N, N_square, N_half, s_x, x_xing, s_y, y_xing, e_0, e_1, e_r, e_q, e_f, e_1_f,
              e_fx, e_K_1, e_y_another, e_x_another, NULL);

    mpz_set_si(zero, 0);
    mpz_set_si(one, 1);
    mpz_set_si(neg_one, -1);
    mpz_set_si(two, 2);
    mpz_set(N, cp.public_key.N);
    mpz_set(N_square, cp.public_key.N_square);
    mpz_set(N_half, cp.public_key.N_half);

    // 计算e_0和e_1
    Enc(e_0, cp.public_key, zero);
    Enc(e_1, cp.public_key, one);

    // cp and csp
    SEQ(e_f, e_x, e_0, cp, csp);

    // cp
    mpz_powm(e_1_f, e_f, neg_one, N_square);
    mpz_mul(e_1_f, e_1, e_1_f);
    mpz_mod(e_1_f, e_1_f, N_square);

    // cp and csp
    SMUL(e_fx, e_f, e_x, cp, csp);
    SMUL(e_y_another, e_f, e_y, cp, csp);

    // cp
    mpz_mul(e_x_another, e_fx, e_1_f);
    mpz_mod(e_x_another, e_x_another, N_square);

    // cp and csp
    SSBA(s_x, x_xing, e_x_another, cp, csp);
    SSBA(s_y, y_xing, e_y_another, cp, csp);

    // 计算enc_q_bits
    vector < mpz_t * > enc_q_bits;
    SBD(enc_q_bits, y_xing, l, cp, csp);

    // 初始化enc_a_bits
    vector < mpz_t * > enc_a_bits;
    for (int i = 0; i < l; ++i) {
        mpz_t *ptr = (mpz_t *) malloc(sizeof(mpz_t));
        mpz_init_set(*ptr, e_0);
        enc_a_bits.push_back(ptr);
    }

    int exponent = l - 1;
    // 循环迭代
    for (int _ = 0; _ < l; ++_) {
        for (int i = 0; i < l - 1; i++) {
            mpz_set(*enc_a_bits[i], *enc_a_bits[i + 1]);
        }
        mpz_set(*enc_a_bits[l - 1], *enc_q_bits[0]);
        for (int i = 0; i < l - 1; i++) {
            mpz_set(*enc_q_bits[i], *enc_q_bits[i + 1]);
        }

        mpz_t e_A, e_Q, e_B;
        mpz_inits(e_A, e_Q, e_B, NULL);
        mpz_set(e_A, *enc_a_bits[l - 1]); // 最后一个e_a，也就是最低位的比特
        for (int i = 0; i < l - 1; ++i) {
            mpz_t temp_e_a;
            mpz_init(temp_e_a);
            mpz_powm_ui(temp_e_a, two, exponent - i, N_square);
            mpz_powm(temp_e_a, *enc_a_bits[i], temp_e_a, N_square);
            mpz_mul(e_A, temp_e_a, e_A);
            mpz_mod(e_A, e_A, N_square);
            mpz_clear(temp_e_a);
        }

        SCMP(e_Q, e_A, x_xing, cp, csp);
        mpz_powm(*enc_q_bits[l - 1], e_Q, neg_one, N_square);
        mpz_mul(*enc_q_bits[l - 1], e_1, *enc_q_bits[l - 1]);
        mpz_mod(*enc_q_bits[l - 1], *enc_q_bits[l - 1], N_square);

        mpz_powm(e_B, x_xing, neg_one, N_square);
        SMUL(e_B, e_B, *enc_q_bits[l - 1], cp, csp);

        mpz_mul(e_A, e_A, e_B);
        mpz_mod(e_A, e_A, N_square);
        SBD(enc_a_bits, e_A, l, cp, csp);
        mpz_clears(e_A, e_Q, e_B, NULL);
    }

    // 计算e_r和e_q
    mpz_set(e_r, *enc_a_bits[l - 1]);
    for (int i = 0; i < l - 1; ++i) {
        mpz_t temp_e_a;
        mpz_init(temp_e_a);
        mpz_powm_ui(temp_e_a, two, exponent - i, N_square);
        mpz_powm(temp_e_a, *enc_a_bits[i], temp_e_a, N_square);
        mpz_mul(e_r, temp_e_a, e_r);
        mpz_mod(e_r, e_r, N_square);
        mpz_clear(temp_e_a);
    }
    mpz_set(e_q, *enc_q_bits[l - 1]);
    for (int i = 0; i < l - 1; ++i) {
        mpz_t temp_e_q;
        mpz_init(temp_e_q);
        mpz_powm_ui(temp_e_q, two, exponent - i, N_square);
        mpz_powm(temp_e_q, *enc_q_bits[i], temp_e_q, N_square);
        mpz_mul(e_q, temp_e_q, e_q);
        mpz_mod(e_q, e_q, N_square);
        mpz_clear(temp_e_q);
    }

    // 计算e_K_1、r_xing和q_xing
    SMUL(e_K_1, s_x, s_y, cp, csp);
    SMUL(result_r_xing, e_r, s_y, cp, csp);
    SMUL(result_q_xing, e_q, e_K_1, cp, csp);

    enc_a_bits.clear();
    enc_q_bits.clear();

    // 清理内存
    mpz_clears(zero, one, neg_one, two, N, N_square, N_half, s_x, x_xing, s_y, y_xing, e_0, e_1, e_r, e_q, e_f, e_1_f,
               e_fx, e_K_1, e_y_another, e_x_another, NULL);
}

/**
 * 乘法正确性检验
 * @param i
 * @param result_cxy
 * @param x
 * @param y
 * @param pai
 */
void check_smul(int i, mpz_t result_cxy, mpz_t x, mpz_t y, Paillier pai) {
    mpz_t z, xy;
    mpz_inits(z, xy, NULL);
    Dec(z, pai.private_key, result_cxy);
    mpz_mul(xy, x, y);
    if (mpz_cmp(z, xy) != 0) {
        gmp_printf("Z = %Zd\n", z);
        gmp_printf("%Zd * %Zd = %Zd\n", x, y, xy);
        printf("SMUL error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(z, xy, NULL);
}

/**
 * 比较正确性检验
 * @param i
 * @param result_cmp_x_y
 * @param x
 * @param y
 * @param pai
 */
void check_scmp(int i, mpz_t result_cmp_x_y, mpz_t x, mpz_t y, Paillier pai) {
    mpz_t z, answer;
    mpz_inits(z, answer, NULL);
    if (mpz_cmp(x, y) >= 0) {
        mpz_set_si(answer, 0);
    } else mpz_set_si(answer, 1);

    Dec(z, pai.private_key, result_cmp_x_y);
    if (mpz_cmp(z, answer) != 0) {
        gmp_printf("dec answer x>=y? = %Zd\n", z);
        gmp_printf("right answer %Zd\n", answer);
        printf("SCMP error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(z, answer, NULL);
}

/**
 * 相等判断正确性检验
 * @param i
 * @param result_eq_x_y
 * @param x
 * @param y
 * @param pai
 */
void check_seq(int i, mpz_t result_eq_x_y, mpz_t x, mpz_t y, Paillier pai) {
    mpz_t z, answer;
    mpz_inits(z, answer, NULL);
    if (mpz_cmp(x, y) == 0) {
        mpz_set_si(answer, 0);
    } else mpz_set_si(answer, 1);

    Dec(z, pai.private_key, result_eq_x_y);
    if (mpz_cmp(z, answer) != 0) {
        gmp_printf("dec answer x==y? = %Zd\n", z);
        gmp_printf("right answer %Zd\n", answer);
        printf("SEQ error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(z, answer, NULL);
}


/**
 * 符号位获取正确性检验，结果与SOCI中SSBA的表示相反
 * @param i
 * @param result_cs
 * @param result_cu
 * @param x
 * @param pai
 */
void check_ssba(int i, mpz_t result_cs, mpz_t result_cu, mpz_t x, Paillier pai) {
    mpz_t s, u, answer_s, answer_u;
    mpz_inits(s, u, answer_s, answer_u, NULL);

    Dec(s, pai.private_key, result_cs);
    Dec(u, pai.private_key, result_cu);

    if (mpz_cmp_si(x, 0) >= 0) {
        mpz_set_si(answer_s, 1);
        mpz_set(answer_u, x);
    } else {
        mpz_set_si(answer_s, 0);
        mpz_set(answer_u, x);
        mpz_neg(answer_u, answer_u); // answer_s = -answer_s
    }

    if (mpz_cmp(answer_u, u) != 0 || mpz_cmp(answer_s, s) != 0) {
        gmp_printf("dec answer s = %Zd u = %Zd\n", s, u);
        gmp_printf("right answer s = %Zd u = %Zd\n", answer_s, answer_u);
        printf("SSBA error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(s, u, answer_s, answer_u, NULL);
}

/**
 * 除法正确性检验
 * @param i
 * @param result_cq
 * @param result_ce
 * @param x
 * @param y
 * @param pai
 */
void check_sdiv(int i, mpz_t result_cq, mpz_t result_ce, mpz_t x, mpz_t y, Paillier pai) {
    mpz_t q, e, right_q, right_e;
    mpz_inits(q, e, right_q, right_e, NULL);
    Dec(q, pai.private_key, result_cq);
    Dec(e, pai.private_key, result_ce);

    mpz_div(right_q, y, x);
    mpz_mod(right_e, y, x);
    if (mpz_cmp(q, right_q) != 0 && mpz_cmp(e, right_e) != 0) {
        gmp_printf("q = %Zd e = %Zd\n", q, e);
        gmp_printf("right_q = %Zd right_e = %Zd\n", right_q, right_e);
        printf("SDIV error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(q, e, right_q, right_e, NULL);
}




void fundamental(const Paillier &pai, const Paillier_Third_Party &cp, const Paillier_Third_Party &csp) {
    mpz_t m, result, ciphertext;
    mpz_inits(m, result, ciphertext, NULL);
    mpz_set_si(m, 4869);
    gmp_printf("m = %Zd\n", m);
    Enc(result, pai.public_key, m);
    mpz_set(ciphertext, result);
    gmp_printf("Enc(m) = %Zd\n", result);
    Dec(result, pai.private_key, result);
    gmp_printf("Dec(m) = %Zd\n", result);

    CR(ciphertext, cp.public_key.N, cp.public_key.N_square);
    gmp_printf("CR(m) = %Zd\n", ciphertext);

    mpz_t result_1, result_2;
    mpz_inits(result_1, result_2, NULL);
    PDec(result_1, cp.private_key, ciphertext);
    PDec(result_2, csp.private_key, ciphertext);
    TDec(result, result_1, result_2, cp.public_key.N_square, cp.public_key.N, cp.public_key.N_half);
    gmp_printf("TDec(m) = %Zd\n", result);
}

#endif

