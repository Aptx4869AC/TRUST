/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#ifndef CSP_PROTOCOL_H
#define CSP_PROTOCOL_H
#pragma once

#include "pcpd.h"
#include "network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;
static mpz_t zero, one, two, neg_one, N, N_square, N_half;

/**
 * 部分解密
 * @param csp
 * @param cp_sock
 */
void CSP_PDec(const Paillier_Third_Party &csp, const int &cp_sock) {
    mpz_t Y, Y2;
    mpz_inits(Y, Y2, NULL);
    recv_mpz(cp_sock, Y);
    PDec(Y2, csp.private_key, Y);
    send_mpz(cp_sock, Y2);
    mpz_clears(Y, Y2, NULL);
}

/**
 * POCF乘法
 * @param csp
 * @param cp_sock
 */
void CSP_SMUL(const Paillier_Third_Party &csp, const int &cp_sock) {
    mpz_t flag, temp, e_r_x, e_r_y, X, Y, X1, Y1, X2, Y2, h, H, S_1, S_2, S_3, e_x_y, r_x, r_y;
    mpz_inits(flag, temp, e_r_x, e_r_y, X, Y, X1, Y1, X2, Y2, h, H, S_1, S_2, S_3, e_x_y, r_x, r_y, NULL);

    recv_mpz(cp_sock, X);
    recv_mpz(cp_sock, Y);
    recv_mpz(cp_sock, X1);
    recv_mpz(cp_sock, Y1);

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

    // 发H给CP
    send_mpz(cp_sock, H);

    // 清理内存
    mpz_clears(tdec_temp_1, tdec_temp_2, NULL);
    mpz_clears(flag, temp, e_r_x, e_r_y, X, Y, X1, Y1, X2, Y2, h, H, S_1, S_2, S_3, r_x, r_y, NULL);
}

/**
 * POCF比较
 * @param csp
 * @param cp_sock
 */
void CSP_SCMP(const Paillier_Third_Party &csp, const int &cp_sock) {
    mpz_t temp, e_l, K, l, u_another;
    mpz_inits(temp, e_l, K, l, u_another, NULL);

    // step2 CSP
    recv_mpz(cp_sock, K);
    recv_mpz(cp_sock, e_l);
    PDec(temp, csp.private_key, e_l);
    TDec(l, K, temp, N_square, N, N_half);
    if (mpz_cmp_si(l, 0) < 0) {
        mpz_add(l, l, N);
    }
    if (mpz_sizeinbase(l, 2) > (mpz_sizeinbase(N, 2) / 2)) {
        Enc(u_another, csp.public_key, one);
    } else {
        Enc(u_another, csp.public_key, zero);
    }
    // send u_another to CP
    send_mpz(cp_sock, u_another);

    mpz_clears(temp, e_l, K, l, u_another, NULL);
}

/**
 * 符号位获取
 * @param csp
 * @param cp_sock
 */
void CSP_SSBA(const Paillier_Third_Party &csp, const int &cp_sock) {
    mpz_t l, temp, e_l, L, e_u;
    mpz_inits(l, temp, e_l, L, e_u, NULL);

    // step 2 csp
    recv_mpz(cp_sock, e_l);
    recv_mpz(cp_sock, L);
    PDec(temp, csp.private_key, e_l);
    TDec(l, L, temp, N_square, N, N_half);
    if (mpz_cmp_si(l, 0) < 0) {
        mpz_add(l, l, N);
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
    send_mpz(cp_sock, e_u);

    mpz_clears(l, temp, e_l, L, e_u, NULL);
}


void CSP_Encrypted_LSB(const Paillier_Third_Party &csp, const int &cp_sock) {

    mpz_t temp, ezero, eone, Y, Y1, Y2, y, alpha;
    mpz_inits(temp, ezero, eone, Y, Y1, Y2, y, alpha, NULL);

    recv_mpz(cp_sock, Y1);
    recv_mpz(cp_sock, Y);
    // Alice:
    PDec(Y2, csp.private_key, Y);
    TDec(y, Y1, Y2, N_square, N, N_half);
    if (mpz_cmp_si(y, 0) < 0) {
        mpz_add(y, y, N);
    }
    if (mpz_even_p(y)) {
        // 如果 y 是偶数
        Enc(alpha, csp.public_key, zero);
    } else {
        Enc(alpha, csp.public_key, one);
    }
    send_mpz(cp_sock, alpha);
    mpz_clears(temp, ezero, eone, Y, Y1, Y2, y, alpha, NULL);
}


void CSP_SVR(const Paillier_Third_Party &csp, const int &cp_sock) {
    mpz_t temp, W, W1, W2, D_w, result, part_result;
    mpz_inits(temp, W, W1, W2, D_w, result, part_result, NULL);

    recv_mpz(cp_sock, W1);
    recv_mpz(cp_sock, W);

    // Alice
    PDec(W2, csp.private_key, W);
    TDec(D_w, W1, W2, N_square, N, N_half);
    if (mpz_cmp_si(D_w, 0) == 0) {
        mpz_set_si(result, 1);
    } else {
        mpz_set_si(result, 0);
    }
    Enc(result, csp.public_key, result);
    PDec(part_result, csp.private_key, result);
    send_mpz(cp_sock, result);
    send_mpz(cp_sock, part_result);
    mpz_clears(temp, W, W1, W2, D_w, result, part_result, NULL);
}


#endif