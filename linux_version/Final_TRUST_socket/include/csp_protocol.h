/******************************************************************************
 * Author: Aptx4869AC 
 * Created: 2024-12-01 
 * GitHub: https://github.com/Aptx4869AC  
 *****************************************************************************/

#ifndef CSP_PROTOCOL_H
#define CSP_PROTOCOL_H
#pragma once

#include "fastPai.h"
#include "network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

/**
 * 部分解密
 * @param csp
 * @param cp_sock
 */
void CSP_PDec(PaillierThdDec &csp, int &cp_sock)
{
    mpz_t Y, Y2;
    mpz_inits(Y, Y2, NULL);
    recv_mpz(cp_sock, Y);
    csp.PDec(Y2, Y);
    send_mpz(cp_sock, Y2);
    mpz_clears(Y, Y2, NULL);
}

/**
 * SOCI-TEE乘法
 * @param csp
 * @param cp_sock
 */
void CSP_TEE_SMUL(PaillierThdDec &csp, int &cp_sock)
{
    mpz_t C, C1, C2, x_add_r1, ey, Y, exy;
    mpz_inits(C, C1, C2, x_add_r1, ey, Y, exy, NULL);
    recv_mpz(cp_sock, C1);
    recv_mpz(cp_sock, C);
    recv_mpz(cp_sock, ey);
    recv_mpz(cp_sock, Y);

    // Step 2, CSP computes
    csp.PDec(C2, C);
    csp.TDec(x_add_r1, C1, C2);
    csp.pai.scl_mul(exy, ey, x_add_r1);
    csp.pai.add(exy, exy, Y);
    /* send <e_xr1_mul_yr2> to cp */
    send_mpz(cp_sock, exy);

    mpz_clears(C, C1, C2, x_add_r1, ey, Y, exy, NULL);
}


/**
 * SOCI-TEE比较
 * @param csp
 * @param cp_sock
 */
void CSP_TEE_SCMP_version2(PaillierThdDec &csp, int &cp_sock)
{
    Offline_val tuple_S1 = csp.pai.offlineVal;
    mpz_t d, D, D1, epi, D2, u0, one_sub_2u0;
    mpz_inits(d, D, D1, epi, D2, u0, one_sub_2u0, NULL);
    recv_mpz(cp_sock, D1);
    recv_mpz(cp_sock, D);
    recv_mpz(cp_sock, epi);

    // Step 2
    csp.PDec(D2, D);
    csp.TDec(d, D1, D2);
    if (mpz_cmp(d, csp.pk.mid) > 0)
    {
        mpz_init_set(u0, tuple_S1.e0);
        mpz_init_set_si(one_sub_2u0, 1);
    } else
    {
        mpz_init_set(u0, tuple_S1.e1);
        mpz_init_set_si(one_sub_2u0, -1);
    }

    mpz_t res;
    mpz_init(res);
    csp.pai.scl_mul(res, epi, one_sub_2u0);
    csp.pai.add(res, res, u0);

    send_mpz(cp_sock, res);

    mpz_clears(d, D, D1, epi, D2, u0, one_sub_2u0, res, NULL);

}

/**
 * SOCI-TEE绝对值
 * @param csp
 * @param cp_sock
 */
void CSP_TEE_SABS(PaillierThdDec &csp, int &cp_sock)
{
    mpz_t d, D, D1, epi, e_x_sub_y, D2, u0, one_sub_twou0, res;
    mpz_inits(d, D, D1, epi, e_x_sub_y, D2, u0, one_sub_twou0, res, NULL);
    recv_mpz(cp_sock, D1);
    recv_mpz(cp_sock, D);
    recv_mpz(cp_sock, e_x_sub_y);

    // Step 2
    csp.PDec(D2, D);
    csp.TDec(d, D1, D2);
    if (mpz_cmp(d, csp.pk.mid) > 0)
    {
        mpz_set_si(u0, 0);
        mpz_set_si(one_sub_twou0, 1);
    } else
    {
        mpz_set_si(u0, 1);
        mpz_set_si(one_sub_twou0, -1);
    }
    csp.pai.scl_mul(res, e_x_sub_y, one_sub_twou0);

    send_mpz(cp_sock, res);

    mpz_clears(d, D, D1, epi, e_x_sub_y, D2, u0, one_sub_twou0, res, NULL);

}

/**
 * SOCI+乘法
 * @param csp
 * @param cp_sock
 */
void CSP_SMUL(PaillierThdDec &csp, int &cp_sock)
{
    mpz_t C, C1, C2, L_mul_xr1_yr2, xr1, yr2, e_xr1_mul_yr2;
    mpz_inits(C, C1, C2, L_mul_xr1_yr2, xr1, yr2, e_xr1_mul_yr2, NULL);
    recv_mpz(cp_sock, C1);
    recv_mpz(cp_sock, C);

    // Step 2, CSP computes
    csp.PDec(C2, C);
    csp.TDec(L_mul_xr1_yr2, C1, C2);
    mpz_divmod(xr1, yr2, L_mul_xr1_yr2, csp.pk.L);
    mpz_mul(e_xr1_mul_yr2, xr1, yr2);
    mpz_mod(e_xr1_mul_yr2, e_xr1_mul_yr2, csp.pk.N);
    csp.pai.encrypt(e_xr1_mul_yr2, e_xr1_mul_yr2);
    /* send <e_xr1_mul_yr2> to cp */
    send_mpz(cp_sock, e_xr1_mul_yr2);

    mpz_clears(C, C1, C2, L_mul_xr1_yr2, xr1, yr2, e_xr1_mul_yr2, NULL);
}

/**
 * SOCI+比较
 * @param csp
 * @param cp_sock
 */
void CSP_SCMP(PaillierThdDec &csp, int &cp_sock)
{

    Offline_val tuple_S1 = csp.pai.offlineVal;

    mpz_t D, D1, D2, d, eu0, e0_csp, e1_csp;
    mpz_inits(D, D1, D2, d, eu0, e0_csp, e1_csp, NULL);
    recv_mpz(cp_sock, D1);
    recv_mpz(cp_sock, D);

    // Step 2
    csp.PDec(D2, D);
    csp.TDec(d, D1, D2);
    if (mpz_cmp(d, csp.pk.mid) > 0)
    {
        csp.pai.add(e0_csp, tuple_S1.e0, tuple_S1.e0);
        mpz_set(eu0, e0_csp);
    } else
    {
        csp.pai.add(e1_csp, tuple_S1.e1, tuple_S1.e0);
        mpz_set(eu0, e1_csp);
    }
    /* send <eu0> to cp */
    send_mpz(cp_sock, eu0);
    mpz_clears(D, D1, D2, d, eu0, e0_csp, e1_csp, NULL);
}

/**
 * TRUST相等判断
 * @param csp
 * @param cp_sock
 */
void CSP_Trust_FEQL(PaillierThdDec &csp, int &cp_sock)
{
    Offline_val tuple_S1 = csp.pai.offlineVal;

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
    mpz_div_ui(half_N, csp.pk.N, 2);

    csp.PDec(ed_1_P2, ed_1);
    csp.TDec(d1, ed_1_P1, ed_1_P2);
    csp.PDec(ed_2_P2, ed_2);
    csp.TDec(d2, ed_2_P1, ed_2_P2);
    mpz_mod(d1, d1, csp.pk.N);
    if (mpz_cmp(d1, half_N) > 0)
    {
        mpz_set_si(miu0, 0);
        csp.pai.add(emiu0, tuple_S1.e0, tuple_S1.e0);   //refresh
    } else if (mpz_cmp(d1, half_N) < 0)
    {
        mpz_set_si(miu0, 1);
        csp.pai.add(emiu0, tuple_S1.e1, tuple_S1.e0);   //refresh
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
        csp.pai.add(emiu0_prime, tuple_S1.e0, tuple_S1.e0);   //refresh
    } else if (mpz_cmp(d2, half_N) < 0)
    {
        mpz_set_si(miu0_prime, 1);
        csp.pai.add(emiu0_prime, tuple_S1.e1, tuple_S1.e0);   //refresh
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
    Offline_val tuple_S1 = csp.pai.offlineVal;


    mpz_t ed, ed_1, em_1, em_1pi;
    mpz_inits(ed, ed_1, em_1, em_1pi, NULL);
    recv_mpz(cp_sock, ed);
    recv_mpz(cp_sock, ed_1);
    recv_mpz(cp_sock, em_1);
    recv_mpz(cp_sock, em_1pi);

    // Step 2
    mpz_t ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N;
    mpz_inits(ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N, NULL);
    mpz_div_ui(half_N, csp.pk.N, 2);
    csp.PDec(ed_2, ed);
    csp.TDec(d, ed_1, ed_2);
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


    mpz_t left, right, res;
    mpz_inits(left, right, res, NULL);
    csp.pai.scl_mul(left, em_1, one_sub_2miu0);
    csp.pai.scl_mul(right, em_1pi, neg_two_sub_4miu0);
    csp.pai.add(res, left, right);
    csp.pai.add(res, res, tuple_S1.e0);   //refresh

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
    Offline_val tuple_S1 = csp.pai.offlineVal;


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
    mpz_div_ui(half_N, csp.pk.N, 2);

    csp.PDec(ed_1_P2, ed_1);
    csp.TDec(d1, ed_1_P1, ed_1_P2);
    csp.PDec(ed_2_P2, ed_2);
    csp.TDec(d2, ed_2_P1, ed_2_P2);

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

    // [R] compute
    mpz_t left, mid, right, res;
    mpz_inits(left, mid, right, res, NULL);
    csp.pai.scl_mul(left, em_3_sub_m_2, miu0_add_miu0_prime); // left
    csp.pai.add(left, em_2, left); // left
    csp.pai.scl_mul(mid, epi_1_m_3_sub_m_2, one_sub_2miu0); // mid
    csp.pai.scl_mul(right, epi_2_m_3_sub_m_2, one_sub_2miu0_prime); // right

    csp.pai.add(res, left, mid);
    csp.pai.add(res, res, right);
    csp.pai.add(res, res, tuple_S1.e0);   //refresh

    /* send [R] to cp */
    send_mpz(cp_sock, res);

    mpz_clears(em_3_sub_m_2, em_2, ed_1, ed_1_P1, epi_1_m_3_sub_m_2, ed_2, ed_2_P1, epi_2_m_3_sub_m_2, NULL);
    mpz_clears(ed_1_P2, ed_2_P2, d1, d2, miu0, miu0_prime, one_sub_2miu0, one_sub_2miu0_prime, miu0_add_miu0_prime, half_N, NULL);
    mpz_clears(left, mid, right, res, NULL);
}


#endif