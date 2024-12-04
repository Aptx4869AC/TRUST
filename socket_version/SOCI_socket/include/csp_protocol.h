/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#ifndef CSP_PROTOCOL_H
#define CSP_PROTOCOL_H
#pragma once

#include "paillier.h"
#include "network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

mpz_t zero, one, neg_one;

/**
 * SOCI乘法
 * @param csp
 * @param cp_sock
 */
void CSP_SMUL(PaillierThdDec &csp, int &cp_sock)
{

    // Step-2 CSP
    mpz_t X, Y, X1, Y1, X2, Y2, maskedX, maskedY, maskedXY, exy;
    mpz_inits(X2, Y2, maskedX, maskedY, maskedXY, exy, NULL);

    // 从CP接收 <X,Y,X1,Y1>
    recv_mpz(cp_sock, X);
    recv_mpz(cp_sock, Y);
    recv_mpz(cp_sock, X1);
    recv_mpz(cp_sock, Y1);

    csp.pdec(X2, X);
    csp.pdec(Y2, Y);
    csp.fdec(maskedX, X1, X2);
    csp.fdec(maskedY, Y1, Y2);
    mpz_mul(maskedXY, maskedX, maskedY);
    mpz_mod(maskedXY, maskedXY, csp.pai.pubkey.n);
    csp.pai.encrypt(exy, maskedXY);

    // CSP发送 <exy> 给CP
    send_mpz(cp_sock, exy);
    mpz_clears(X, Y, X1, Y1, X2, Y2, maskedX, maskedY, maskedXY, exy, NULL);
}

/**
 * SOCI比较
 * @param csp
 * @param cp_sock
 */
void CSP_SCMP(PaillierThdDec &csp, int &cp_sock)
{

    //Step-2
    mpz_t d, D2, D1, D, res;
    mpz_inits(d, D2, D1, D, res, NULL);
    recv_mpz(cp_sock, D1);
    recv_mpz(cp_sock, D);

    csp.pdec(D2, D);
    csp.fdec(d, D1, D2);

    if (mpz_cmp(d, csp.pai.pubkey.half_n) > 0)
    {
        csp.pai.encrypt(res, zero);
    } else
    {
        csp.pai.encrypt(res, one);
    }
    send_mpz(cp_sock, res);

    mpz_clears(d, D2, D1, D, res, NULL);
}

/**
 * TRUST相等判断
 * @param csp
 * @param cp_sock
 */
void CSP_Trust_FEQL(PaillierThdDec &csp, int &cp_sock)
{
    mpz_t ed_1, ed_1_P1, epi_1, ed_2, ed_2_P1, epi_2;
    mpz_inits(ed_1, ed_1_P1, epi_1, ed_2, ed_2_P1, epi_2, NULL);
    recv_mpz(cp_sock, ed_1);
    recv_mpz(cp_sock, ed_1_P1);
    recv_mpz(cp_sock, epi_1);
    recv_mpz(cp_sock, ed_2);
    recv_mpz(cp_sock, ed_2_P1);
    recv_mpz(cp_sock, epi_2);


    // Step 2
    mpz_t ed_1_P2, ed_2_P2, d1, d2, miu0, emiu0, miu0_prime, emiu0_prime, one_sub_2miu0, one_sub_2miu0_prime, half_N;
    mpz_inits(ed_1_P2, ed_2_P2, d1, d2, miu0, emiu0, miu0_prime, emiu0_prime, one_sub_2miu0, one_sub_2miu0_prime, half_N, NULL);
    mpz_div_ui(half_N, csp.pai.pubkey.n, 2);

    csp.pdec(ed_1_P2, ed_1);
    csp.fdec(d1, ed_1_P1, ed_1_P2);
    csp.pdec(ed_2_P2, ed_2);
    csp.fdec(d2, ed_2_P1, ed_2_P2);
    mpz_mod(d1, d1, csp.pai.pubkey.n);
    if (mpz_cmp(d1, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
        csp.pai.encrypt(emiu0, miu0);
    } else if (mpz_cmp(d1, half_N) < 0)
    {
        mpz_set_si(miu0, 1);
        csp.pai.encrypt(emiu0, miu0);
    } else
    {
        cout << " ERROR !\n";
        throw;
    }

    mpz_mul_si(one_sub_2miu0, miu0, 2);
    mpz_neg(one_sub_2miu0, one_sub_2miu0);
    mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);
    if (mpz_cmp(d2, half_N) > 0)
    {
        mpz_set_si(miu0_prime, 0);
        csp.pai.encrypt(emiu0_prime, miu0_prime);
    } else if (mpz_cmp(d2, half_N) < 0)
    {
        mpz_set_si(miu0_prime, 1);
        csp.pai.encrypt(emiu0_prime, miu0_prime);
    } else
    {
        cout << " ERROR !\n";
        throw;
    }
    mpz_mul_si(one_sub_2miu0_prime, miu0_prime, 2);
    mpz_neg(one_sub_2miu0_prime, one_sub_2miu0_prime);
    mpz_add_ui(one_sub_2miu0_prime, one_sub_2miu0_prime, 1);

    mpz_t left, right, res;
    mpz_inits(left, right, res, NULL);
    csp.pai.scl_mul(left, epi_1, one_sub_2miu0);
    csp.pai.add(left, emiu0, left);
    csp.pai.scl_mul(right, epi_2, one_sub_2miu0_prime);
    csp.pai.add(right, emiu0_prime, right);
    csp.pai.add(res, left, right);

    /* send [R] to cp */
    send_mpz(cp_sock, res);

    mpz_clears(ed_1, ed_1_P1, epi_1, ed_2, ed_2_P1, epi_2, NULL);
    mpz_clears(ed_1_P2, ed_2_P2, d1, d2, miu0, emiu0, miu0_prime, emiu0_prime, one_sub_2miu0, one_sub_2miu0_prime, half_N, NULL);
    mpz_clears(left, right, res, NULL);
}

/**
 * TRUST绝对值
 * @param csp
 * @param cp_sock
 */
void CSP_Trust_FABS(PaillierThdDec &csp, int &cp_sock)
{

    mpz_t ed, ed_1, em_1, em_1pi;
    mpz_inits(ed, ed_1, em_1, em_1pi, NULL);
    recv_mpz(cp_sock, ed);
    recv_mpz(cp_sock, ed_1);
    recv_mpz(cp_sock, em_1);
    recv_mpz(cp_sock, em_1pi);

    // Step 2
    mpz_t ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N;
    mpz_inits(ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N, NULL);
    mpz_div_ui(half_N, csp.pai.pubkey.n, 2);
    csp.pdec(ed_2, ed);
    csp.fdec(d, ed_1, ed_2);
    if (mpz_cmp(d, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
    } else if (mpz_cmp(d, half_N) < 0)
    {
        mpz_set_si(miu0, 1);
    } else
    {
        cout << " ERROR !\n";
        throw;
    }
    mpz_mul_si(one_sub_2miu0, miu0, 2);
    mpz_neg(one_sub_2miu0, one_sub_2miu0);
    mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);
    mpz_mul_ui(neg_two_sub_4miu0, one_sub_2miu0, 2);
    mpz_neg(neg_two_sub_4miu0, neg_two_sub_4miu0);


    // ezero
    mpz_t ezero;
    mpz_init(ezero);
    csp.pai.encrypt(ezero, zero);

    mpz_t left, right, res;
    mpz_inits(left, right, res, NULL);
    csp.pai.scl_mul(left, em_1, one_sub_2miu0);
    csp.pai.scl_mul(right, em_1pi, neg_two_sub_4miu0);
    csp.pai.add(res, left, right);
    csp.pai.add(res, res, ezero);   //refresh

    /* send [R] to cp */
    send_mpz(cp_sock, res);

    mpz_clears(ed, ed_1, em_1, em_1pi, NULL);
    mpz_clears(ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N, NULL);
    mpz_clears(left, right, res, NULL);
}

/**
 * TRUST三目运算
 * @param csp
 * @param cp_sock
 */
void CSP_Trust_FTRN_version1(PaillierThdDec &csp, int &cp_sock)
{

    mpz_t em_3_sub_m_2, em_2, ed_1, ed_1_P1, epi_1_m_3_sub_m_2, ed_2, ed_2_P1, epi_2_m_3_sub_m_2;
    mpz_inits(em_3_sub_m_2, em_2, ed_1, ed_1_P1, epi_1_m_3_sub_m_2, ed_2, ed_2_P1, epi_2_m_3_sub_m_2, NULL);
    recv_mpz(cp_sock, em_3_sub_m_2);
    recv_mpz(cp_sock, em_2);
    recv_mpz(cp_sock, ed_1);
    recv_mpz(cp_sock, ed_1_P1);
    recv_mpz(cp_sock, epi_1_m_3_sub_m_2);
    recv_mpz(cp_sock, ed_2);
    recv_mpz(cp_sock, ed_2_P1);
    recv_mpz(cp_sock, epi_2_m_3_sub_m_2);

    // Step 2
    mpz_t ed_1_P2, ed_2_P2, d1, d2, miu0, miu0_prime, one_sub_2miu0, one_sub_2miu0_prime, miu0_add_miu0_prime, half_N;
    mpz_inits(ed_1_P2, ed_2_P2, d1, d2, miu0, miu0_prime, one_sub_2miu0, one_sub_2miu0_prime, miu0_add_miu0_prime, half_N, NULL);
    mpz_div_ui(half_N, csp.pai.pubkey.n, 2);

    csp.pdec(ed_1_P2, ed_1);
    csp.fdec(d1, ed_1_P1, ed_1_P2);
    csp.pdec(ed_2_P2, ed_2);
    csp.fdec(d2, ed_2_P1, ed_2_P2);

    if (mpz_cmp(d1, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
    } else if (mpz_cmp(d1, half_N) < 0)
    {
        mpz_set_si(miu0, 1);
    } else
    {
        cout << " ERROR !\n";
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
        cout << " ERROR !\n";
        throw;
    }
    mpz_mul_si(one_sub_2miu0_prime, miu0_prime, 2);
    mpz_neg(one_sub_2miu0_prime, one_sub_2miu0_prime);
    mpz_add_ui(one_sub_2miu0_prime, one_sub_2miu0_prime, 1);
    mpz_add(miu0_add_miu0_prime, miu0, miu0_prime);

    mpz_t ezero;
    mpz_init(ezero);
    csp.pai.encrypt(ezero, zero);

    // [R] compute
    mpz_t left, mid, right, res;
    mpz_inits(left, mid, right, res, NULL);
    csp.pai.scl_mul(left, em_3_sub_m_2, miu0_add_miu0_prime); // left
    csp.pai.add(left, em_2, left); // left
    csp.pai.scl_mul(mid, epi_1_m_3_sub_m_2, one_sub_2miu0); // mid
    csp.pai.scl_mul(right, epi_2_m_3_sub_m_2, one_sub_2miu0_prime); // right

    csp.pai.add(res, left, mid);
    csp.pai.add(res, res, right);
    csp.pai.add(res, res, ezero);   //refresh

    /* send [R] to cp */
    send_mpz(cp_sock, res);

    mpz_clears(em_3_sub_m_2, em_2, ed_1, ed_1_P1, epi_1_m_3_sub_m_2, ed_2, ed_2_P1, epi_2_m_3_sub_m_2, NULL);
    mpz_clears(ed_1_P2, ed_2_P2, d1, d2, miu0, miu0_prime, one_sub_2miu0, one_sub_2miu0_prime, miu0_add_miu0_prime, half_N, NULL);
    mpz_clears(left, mid, right, res, NULL);
}

/**
 * TRUST三目运算
 * @param csp
 * @param cp_sock
 */
void CSP_Trust_FTRN_version2(PaillierThdDec &csp, int &cp_sock)
{
    mpz_t ed, ed_1, em_2, em_3, epi_m_3_sub_m_2;
    mpz_inits(ed, ed_1, em_2, em_3, epi_m_3_sub_m_2, NULL);
    recv_mpz(cp_sock, ed);
    recv_mpz(cp_sock, ed_1);
    recv_mpz(cp_sock, em_2);
    recv_mpz(cp_sock, em_3);
    recv_mpz(cp_sock, epi_m_3_sub_m_2);

    // Step 2
    mpz_t ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N;
    mpz_inits(ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N, NULL);
    mpz_div_ui(half_N, csp.pai.pubkey.n, 2);

    csp.pdec(ed_2, ed);
    csp.fdec(d, ed_1, ed_2);
    if (mpz_cmp(d, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
    } else if (mpz_cmp(d, half_N) < 0)
    {
        mpz_set_si(miu0, 1);
    } else
    {
        cout << " ERROR !\n";
        throw;
    }

    mpz_neg(one_sub_miu0, miu0);
    mpz_add_ui(one_sub_miu0, one_sub_miu0, 1);

    mpz_mul_si(one_sub_2miu0, miu0, 2);
    mpz_neg(one_sub_2miu0, one_sub_2miu0);
    mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);

    mpz_t ezero;
    mpz_init(ezero);
    csp.pai.encrypt(ezero, zero);

    // [R] compute
    mpz_t left, mid, right, res;
    mpz_inits(left, mid, right, res, NULL);
    csp.pai.scl_mul(left, em_2, one_sub_miu0); // left
    csp.pai.scl_mul(mid, em_3, miu0); // mid
    csp.pai.scl_mul(right, epi_m_3_sub_m_2, one_sub_2miu0); // right

    csp.pai.add(res, left, mid);
    csp.pai.add(res, res, right);
    csp.pai.add(res, res, ezero);   //refresh

    /* send [R] to cp */
    send_mpz(cp_sock, res);

    mpz_clears(ed, ed_1, em_2, em_3, epi_m_3_sub_m_2, NULL);
    mpz_clears(ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N, NULL);
    mpz_clears(left, mid, right, res, NULL);
}

#endif