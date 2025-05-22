/******************************************************************************
 * Author: Aptx4869AC 
 * Created: 2024-12-01 
 * GitHub: https://github.com/Aptx4869AC  
 *****************************************************************************/

#ifndef CP_PROTOCOL_H
#define CP_PROTOCOL_H
#pragma once

#include "fastPai.h"
#include "network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

ssize_t total_bytes = 0;

/**
 * 部分解密
 * @param result
 * @param Y
 */
void CSP_PDec(mpz_t result, mpz_t &Y)
{
    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);
    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 3);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    send_mpz(csp_sock, Y);
    recv_mpz(csp_sock, result);
    close(csp_sock);
}

/**
 * SOCI+乘法
 * @param res
 * @param ex
 * @param ey
 * @param cp
 */
void CP_SMUL(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp)
{

    Offline_val tuple_S0 = cp.pai.offlineVal;

    // Step 1, CP computes
    mpz_t er1, er2, neg_er1r2;
    mpz_inits(er1, er2, neg_er1r2, NULL);
    cp.pai.add(er1, tuple_S0.er1, tuple_S0.e0);
    cp.pai.add(er2, tuple_S0.er2, tuple_S0.e0);
    cp.pai.add(neg_er1r2, tuple_S0.neg_er1r2, tuple_S0.e0);
    mpz_t X, Y, C, C1;
    mpz_inits(X, Y, C, C1, NULL);
    cp.pai.add(X, ex, er1);
    cp.pai.add(Y, ey, er2);
    cp.pai.scl_mul(C, X, cp.pk.L);
    cp.pai.add(C, C, Y);
    cp.PDec(C1, C);
    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 1);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    /* send <C, C1> to csp */
    send_mpz(csp_sock, C1);
    send_mpz(csp_sock, C);


    // Step 3, CP computes
    mpz_t neg_er2x, neg_er1y, exy, neg_r1, neg_r2, e_xr1_mul_yr2;
    mpz_inits(neg_er2x, neg_er1y, exy, neg_r1, neg_r2, e_xr1_mul_yr2, NULL);
    mpz_neg(neg_r1, tuple_S0.r1);
    mpz_neg(neg_r2, tuple_S0.r2);
    cp.pai.scl_mul(neg_er2x, ex, neg_r2);
    cp.pai.scl_mul(neg_er1y, ey, neg_r1);
    recv_mpz(csp_sock, e_xr1_mul_yr2);

    cp.pai.add(exy, e_xr1_mul_yr2, neg_er2x);
    cp.pai.add(exy, exy, neg_er1y);
    cp.pai.add(exy, exy, neg_er1r2);
    mpz_set(res, exy);

    mpz_clears(er1, er2, neg_er1r2, X, Y, C, C1, neg_er2x, neg_er1y, exy, neg_r1, neg_r2, e_xr1_mul_yr2, NULL);
    close(csp_sock);
}

/**
 * SOCI+比较
 * @param res
 * @param ex
 * @param ey
 * @param cp
 */
void CP_SCMP(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

    mpz_t er4, er3r4, D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1;
    mpz_inits(er4, er3r4, D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1, NULL);
    cp.pai.add(er4, tuple_S0.er4, tuple_S0.e0); //refresh
    cp.pai.add(er3r4, tuple_S0.er3r4, tuple_S0.e0); //refresh

    int pi = (int) rand() % 2;
    mpz_neg(neg_r1, tuple_S0.r3);

    // Step 1
    if (pi == 0)
    {
        cp.pai.scl_mul(exr1, ex, tuple_S0.r3);
        cp.pai.scl_mul(neg_eyr1, ey, neg_r1);
        cp.pai.add(D, exr1, neg_eyr1);
        cp.pai.add(D, D, er3r4);
    } else        // pi == 1
    {
        cp.pai.scl_mul(eyr1, ey, tuple_S0.r3);
        cp.pai.scl_mul(neg_exr1, ex, neg_r1);
        cp.pai.add(D, eyr1, neg_exr1);
        cp.pai.add(D, D, er4);
    }
    cp.PDec(D1, D);

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 2);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    send_mpz(csp_sock, D1);
    send_mpz(csp_sock, D);

    // Step 3
    if (pi == 1)
    {
        mpz_t e1_cp;
        mpz_init(e1_cp);
        cp.pai.add(e1_cp, tuple_S0.e1, tuple_S0.e0);
        recv_mpz(csp_sock, res);
        close(csp_sock);
        cp.pai.sub(res, tuple_S0.e1, res);
        mpz_clear(e1_cp);
    } else
    {
        recv_mpz(csp_sock, res);
        close(csp_sock);
    }

    mpz_clears(er4, er3r4, D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1, NULL);

}

/**
 * SOCI+符号位获取
 * @param s_x
 * @param u_x
 * @param ex
 * @param cp
 */
void SSBA(mpz_t s_x, mpz_t u_x, mpz_t ex, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

    mpz_t e0, e1;
    mpz_inits(e0, e1, NULL);
    cp.pai.add(e0, tuple_S0.e0, tuple_S0.e0);//refresh
    cp.pai.add(e1, tuple_S0.e1, tuple_S0.e0);//refresh

    // Step-1
    CP_SCMP(s_x, ex, e0, cp);

    // Step-2
    mpz_t n_sub_2;
    mpz_init(n_sub_2);
    mpz_sub_ui(n_sub_2, cp.pai.pk.N, 2);

    mpz_t sign;
    mpz_init(sign);
    cp.pai.scl_mul(sign, s_x, n_sub_2);   // [s_x]^(N-2)
    cp.pai.add(sign, e1, sign);    // [1]*[s_x^(N-2)]

    // Step-3
    CP_SMUL(u_x, sign, ex, cp);

    mpz_clears(e0, e1, n_sub_2, sign, NULL);

}

/**
 * SOCI+除法
 * @param eq
 * @param ee
 * @param ex
 * @param ey
 * @param ell
 * @param cp
 */
void SDIV(mpz_t eq, mpz_t ee, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

    mpz_t e0, e1;
    mpz_inits(e0, e1, NULL);
    cp.pai.add(e0, tuple_S0.e0, tuple_S0.e0);//refresh
    cp.pai.add(e1, tuple_S0.e1, tuple_S0.e0);//refresh

    // Step-1
    mpz_set(eq, e0);

    mpz_t c, u, e, ue, m, two, neg_one;
    mpz_inits(c, u, e, ue, m, two, neg_one, NULL);
    mpz_set_si(two, 2);
    mpz_set(ee, ex);
    mpz_set_si(neg_one, -1);

    for (int i = ell; i >= 0; i--)
    {
        // Step-2
        mpz_pow_ui(e, two, i);      // e=2^i
        cp.pai.scl_mul(c, ey, e);   // [y]^{2^i}

        // Step-3
        CP_SCMP(u, ee, c, cp);

        // Step-4
        cp.pai.scl_mul(u, u, neg_one);
        cp.pai.add(u, e1, u);
        cp.pai.scl_mul(ue, u, e);
        cp.pai.add(eq, eq, ue);

        // Step-5
        CP_SMUL(m, u, c, cp);

        // Step-6
        cp.pai.scl_mul(m, m, neg_one);
        cp.pai.add(ee, ee, m);
    }
    mpz_clears(c, u, e, ue, m, two, neg_one, NULL);
}

/**
 * SOCI-TEE乘法
 * @param res
 * @param ex
 * @param ey
 * @param cp
 */
void CP_TEE_SMUL(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp)
{

    Offline_val tuple_S0 = cp.pai.offlineVal;

    // Step 1, CP computes
    mpz_t neg_r1, er1;
    mpz_inits(neg_r1, er1, NULL);
    cp.pai.add(er1, tuple_S0.er1, tuple_S0.e0);
    mpz_neg(neg_r1, tuple_S0.r1);
    mpz_t X, Y, C, C1;
    mpz_inits(X, Y, C, C1, NULL);
    cp.pai.add(X, ex, er1);
    cp.pai.scl_mul(Y, ey, neg_r1);
    mpz_set(C, X);
    cp.PDec(C1, C);

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 4);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    /* send <C, C1, [y], [-ry]> to csp */
    send_mpz(csp_sock, C1);
    send_mpz(csp_sock, C);
    send_mpz(csp_sock, ey);
    send_mpz(csp_sock, Y);


    // Step 3, CP computes
    recv_mpz(csp_sock, res);
    close(csp_sock);

    mpz_clears(neg_r1, X, Y, C, C1, er1, NULL);

}


/**
 * SOCI-TEE比较
 * @param res
 * @param ex
 * @param ey
 * @param cp
 */
void CP_TEE_SCMP_version2(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

    mpz_t er4, er3r4;
    mpz_inits(er4, er3r4, NULL);

    int pi = (int) rand() % 2;
    mpz_t D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1;
    mpz_inits(D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1, NULL);
    mpz_neg(neg_r1, tuple_S0.r3);

    // Step 1
    if (pi == 0)
    {
        cp.pai.scl_mul(exr1, ex, tuple_S0.r3);
        cp.pai.scl_mul(neg_eyr1, ey, neg_r1);
        cp.pai.add(D, exr1, neg_eyr1);
        cp.pai.add(er3r4, tuple_S0.er3r4, tuple_S0.e0);//refresh
        cp.pai.add(D, D, er3r4);
    } else        // pi == 1
    {
        cp.pai.scl_mul(eyr1, ey, tuple_S0.r3);
        cp.pai.scl_mul(neg_exr1, ex, neg_r1);
        cp.pai.add(D, eyr1, neg_exr1);
        cp.pai.add(er4, tuple_S0.er4, tuple_S0.e0);//refresh
        cp.pai.add(D, D, er4);
    }
    cp.PDec(D1, D);
    mpz_t epi;
    mpz_init(epi);
    if (pi == 0)
    {
        mpz_set(epi, tuple_S0.e0);
        cp.pai.add(epi, tuple_S0.e0, tuple_S0.e0);//refresh
    } else
    {
        mpz_set(epi, tuple_S0.e1);
        cp.pai.add(epi, tuple_S0.e1, tuple_S0.e0);//refresh
    }

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 6);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    send_mpz(csp_sock, D1);
    send_mpz(csp_sock, D);
    send_mpz(csp_sock, epi);

    // Step 3
    recv_mpz(csp_sock, res);
    close(csp_sock);

    mpz_clears(er4, er3r4, D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1, epi, NULL);

}



/**
 * SOCI-TEE除法
 * @param eq
 * @param ee
 * @param ex
 * @param ey
 * @param ell
 * @param cp
 */
void TEE_SDIV(mpz_t eq, mpz_t ee, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

    mpz_t e0, e1;
    mpz_inits(e0, e1, NULL);
    cp.pai.add(e0, tuple_S0.e0, tuple_S0.e0);//refresh
    cp.pai.add(e1, tuple_S0.e1, tuple_S0.e0);//refresh

    // Step-1
    mpz_set(eq, e0);

    mpz_t c, u, e, ue, m, two, neg_one;
    mpz_inits(c, u, e, ue, m, two, neg_one, NULL);
    mpz_set_si(two, 2);
    mpz_set(ee, ex);
    mpz_set_si(neg_one, -1);

    for (int i = ell; i >= 0; i--)
    {
        // Step-2
        mpz_pow_ui(e, two, i);      // e=2^i
        cp.pai.scl_mul(c, ey, e);   // [y]^{2^i}

        // Step-3
        CP_SCMP(u, ee, c, cp);

        // Step-4
        cp.pai.scl_mul(u, u, neg_one);
        cp.pai.add(u, e1, u);
        cp.pai.scl_mul(ue, u, e);
        cp.pai.add(eq, eq, ue);

        // Step-5
        CP_TEE_SMUL(m, u, c, cp);

        // Step-6
        cp.pai.scl_mul(m, m, neg_one);
        cp.pai.add(ee, ee, m);
    }
    mpz_clears(c, u, e, ue, m, two, neg_one, NULL);
}


/**
 * TRUST相等判断
 * @param res
 * @param em_1
 * @param em_2
 * @param cp
 */
void CP_Trust_FEQL(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp)
{

    Offline_val tuple_S0 = cp.pai.offlineVal;

    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_t ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1;
    mpz_inits(ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);
    mpz_set(r1, tuple_S0.trust_r1);
    mpz_neg(neg_r1, r1);

    // Step 1
    int pi_1 = (int) rand() % 2;
    if (pi_1 == 0)
    {
        // [r1(m1-m2) + r1+r2]
        cp.pai.scl_mul(em_1r1, em_1, r1);
        cp.pai.scl_mul(neg_em_2r1, em_2, neg_r1);
        cp.pai.add(ed_1, em_1r1, neg_em_2r1);
        cp.pai.add(e_trust_r1_add_r2, tuple_S0.e_trust_r1_add_r2, tuple_S0.e0);  //refresh
        cp.pai.add(ed_1, ed_1, e_trust_r1_add_r2);
    } else
    {
        // [r1(m2-m1) + r2]
        cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
        cp.pai.scl_mul(em_2r1, em_2, r1);
        cp.pai.add(ed_1, neg_em_1r1, em_2r1);
        cp.pai.add(e_trust_r2, tuple_S0.e_trust_r2, tuple_S0.e0);  //refresh
        cp.pai.add(ed_1, ed_1, e_trust_r2);
    }
    cp.PDec(ed_1_P1, ed_1);
    mpz_t epi_1;
    mpz_inits(epi_1, NULL);
    if (pi_1 == 0)
    {
        cp.pai.add(epi_1, tuple_S0.e0, tuple_S0.e0);   //refresh

    } else
    {
        cp.pai.add(epi_1, tuple_S0.e1, tuple_S0.e0);   //refresh
    }

    mpz_t e_trust_r2_prime, e_trust_r1_prime_add_r2_prime;
    mpz_inits(e_trust_r2_prime, e_trust_r1_prime_add_r2_prime, NULL);
    mpz_t ed_2, ed_2_P1, em_1r1_prime, em_2r1_prime, r1_prime, neg_r1_prime, neg_em_2r1_prime, neg_em_1r1_prime;
    mpz_inits(ed_2, ed_2_P1, em_1r1_prime, em_2r1_prime, r1_prime, neg_r1_prime, neg_em_2r1_prime, neg_em_1r1_prime, NULL);
    mpz_set(r1_prime, tuple_S0.trust_r1_prime);
    mpz_neg(neg_r1_prime, r1_prime);

    int pi_2 = (int) rand() % 2;
    // Step 1
    if (pi_2 == 0)
    {
        // [r1_prime(m2-m1) + r1_prime+r2_prime]
        cp.pai.scl_mul(em_2r1_prime, em_2, r1_prime);
        cp.pai.scl_mul(neg_em_1r1_prime, em_1, neg_r1_prime);
        cp.pai.add(ed_2, neg_em_1r1_prime, em_2r1_prime);
        cp.pai.add(e_trust_r1_prime_add_r2_prime, tuple_S0.e_trust_r1_prime_add_r2_prime, tuple_S0.e0);  //refresh
        cp.pai.add(ed_2, ed_2, e_trust_r1_prime_add_r2_prime);
    } else
    {
        // [r1_prime(m1-m2) + r2_prime]
        cp.pai.scl_mul(neg_em_2r1_prime, em_2, neg_r1_prime);
        cp.pai.scl_mul(em_1r1_prime, em_1, r1_prime);
        cp.pai.add(ed_2, em_1r1_prime, neg_em_2r1_prime);
        cp.pai.add(e_trust_r2_prime, tuple_S0.e_trust_r2_prime, tuple_S0.e0);  //refresh
        cp.pai.add(ed_2, ed_2, e_trust_r2_prime);
    }
    cp.PDec(ed_2_P1, ed_2);
    mpz_t epi_2;
    mpz_inits(epi_2, NULL);
    if (pi_2 == 0)
    {
        cp.pai.add(epi_2, tuple_S0.e0, tuple_S0.e0);   //refresh

    } else
    {
        cp.pai.add(epi_2, tuple_S0.e1, tuple_S0.e0);   //refresh
    }

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 8);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    /* send <[d_1], (PDec)[d_1], [pi_1], [d_2], (PDec)[d_2], [pi_2]> to csp */
    send_mpz(csp_sock, ed_1);
    send_mpz(csp_sock, ed_1_P1);
    send_mpz(csp_sock, epi_1);
    send_mpz(csp_sock, ed_2);
    send_mpz(csp_sock, ed_2_P1);
    send_mpz(csp_sock, epi_2);

    // Step 3, CP get
    recv_mpz(csp_sock, res);
    close(csp_sock);


    mpz_clears(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_clears(ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);
    mpz_clears(epi_1, NULL);
    mpz_clears(e_trust_r2_prime, e_trust_r1_prime_add_r2_prime, NULL);
    mpz_clears(ed_2, ed_2_P1, em_1r1_prime, em_2r1_prime, r1_prime, neg_r1_prime, neg_em_2r1_prime, neg_em_1r1_prime, NULL);
    mpz_clears(epi_2, NULL);
}

/**
 * TRUST绝对值
 * @param eres
 * @param em_1
 * @param cp
 */
void CP_Trust_FABS(mpz_t eres, mpz_t em_1, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;
    
    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);

    mpz_t ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1;
    mpz_inits(ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
    mpz_set(r1, tuple_S0.trust_r1);
    mpz_neg(neg_r1, r1);

    int pi = (int) rand() % 2;
    // Step 1
    if (pi == 0)
    {
        cp.pai.scl_mul(em_1r1, em_1, r1);
        cp.pai.add(e_trust_r1_add_r2, tuple_S0.e_trust_r1_add_r2, tuple_S0.e0);  //refresh
        cp.pai.add(ed, em_1r1, e_trust_r1_add_r2);
    } else
    {
        cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
        cp.pai.add(e_trust_r2, tuple_S0.e_trust_r2, tuple_S0.e0);  //refresh
        cp.pai.add(ed, neg_em_1r1, e_trust_r2);
    }
    cp.PDec(ed_1, ed);

    mpz_t em_1pi;
    mpz_inits(em_1pi, NULL);
    cp.pai.scl_mul(em_1pi, em_1, pi);
    cp.pai.add(em_1pi, em_1pi, tuple_S0.e0);   //refresh


    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 9);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    /* send <[d], (PDec)[d], [pai * m_1], [m_1]> to csp */
    send_mpz(csp_sock, ed);
    send_mpz(csp_sock, ed_1);
    send_mpz(csp_sock, em_1);
    send_mpz(csp_sock, em_1pi);

    // Step 3, CP get
    recv_mpz(csp_sock, eres);
    close(csp_sock);

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_clears(ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
    mpz_clears(em_1pi, NULL);
}

/**
 * TRUST三目运算
 * @param eres
 * @param em_1
 * @param em_2
 * @param em_3
 * @param cp
 */
void CP_Trust_FTRN_version1(mpz_t eres, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;


    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_t ed_1, ed_1_P1, em_1r1, r1, neg_r1, neg_em_1r1;
    mpz_inits(ed_1, ed_1_P1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
    mpz_set(r1, tuple_S0.trust_r1);
    mpz_neg(neg_r1, r1);
    // Step 1

    int pi_1 = rand() % 2;
    if (pi_1 == 0)
    {
        // [r1*m1 + r2]
        cp.pai.scl_mul(em_1r1, em_1, r1);
        cp.pai.add(e_trust_r2, tuple_S0.e_trust_r2, tuple_S0.e0);  //refresh
        cp.pai.add(ed_1, em_1r1, e_trust_r2);
    } else
    {
        // [-r1*m1 + r1+r2]
        cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
        cp.pai.add(e_trust_r1_add_r2, tuple_S0.e_trust_r1_add_r2, tuple_S0.e0);  //refresh
        cp.pai.add(ed_1, neg_em_1r1, e_trust_r1_add_r2);
    }
    cp.PDec(ed_1_P1, ed_1);
    mpz_t e_trust_2r1_prime_add_r2_prime, e_trust_r2_prime_sub_r1_prime;
    mpz_inits(e_trust_2r1_prime_add_r2_prime, e_trust_r2_prime_sub_r1_prime, NULL);

    mpz_t ed_2, ed_2_P1, em_1r1_prime, r1_prime, neg_r1_prime, neg_em_1r1_prime;
    mpz_inits(ed_2, ed_2_P1, em_1r1_prime, r1_prime, neg_r1_prime, neg_em_1r1_prime, NULL);
    mpz_set(r1_prime, tuple_S0.trust_r1_prime);
    mpz_neg(neg_r1_prime, r1_prime);

    int pi_2 = rand() % 2;
    if (pi_2 == 0)
    {
        cp.pai.scl_mul(neg_em_1r1_prime, em_1, neg_r1_prime);
        cp.pai.add(e_trust_2r1_prime_add_r2_prime, tuple_S0.e_trust_2r1_prime_add_r2_prime, tuple_S0.e0);  //refresh
        cp.pai.add(ed_2, neg_em_1r1_prime, e_trust_2r1_prime_add_r2_prime);
    } else
    {

        cp.pai.scl_mul(em_1r1_prime, em_1, r1_prime);
        cp.pai.add(e_trust_r2_prime_sub_r1_prime, tuple_S0.e_trust_r2_prime_sub_r1_prime, tuple_S0.e0);  //refresh
        cp.pai.add(ed_2, em_1r1_prime, e_trust_r2_prime_sub_r1_prime);
    }

    cp.PDec(ed_2_P1, ed_2);
    mpz_t em_3_sub_m_2, epi_1_m_3_sub_m_2, epi_2_m_3_sub_m_2;
    mpz_inits(em_3_sub_m_2, epi_1_m_3_sub_m_2, epi_2_m_3_sub_m_2, NULL);

    // [m3-m2]
    cp.pai.scl_mul(em_3_sub_m_2, em_2, -1);
    cp.pai.add(em_3_sub_m_2, em_3, em_3_sub_m_2);
    // [pi_1 * (m3-m2)]
    cp.pai.scl_mul(epi_1_m_3_sub_m_2, em_3_sub_m_2, pi_1);
    cp.pai.add(epi_1_m_3_sub_m_2, epi_1_m_3_sub_m_2, tuple_S0.e0);   //refresh

    // [pi_2 * (m3-m2)]
    cp.pai.scl_mul(epi_2_m_3_sub_m_2, em_3_sub_m_2, pi_2);
    cp.pai.add(epi_2_m_3_sub_m_2, epi_2_m_3_sub_m_2, tuple_S0.e0);   //refresh

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 10);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    /* send <[m3-m2], [m2], [d_1], (PDec)[d_1], [pi_1 * (m3-m2)], [d_2], (PDec)[d_2], [pi_2 * (m3-m2)]> to csp */
    send_mpz(csp_sock, em_3_sub_m_2);
    send_mpz(csp_sock, em_2);
    send_mpz(csp_sock, ed_1);
    send_mpz(csp_sock, ed_1_P1);
    send_mpz(csp_sock, epi_1_m_3_sub_m_2);
    send_mpz(csp_sock, ed_2);
    send_mpz(csp_sock, ed_2_P1);
    send_mpz(csp_sock, epi_2_m_3_sub_m_2);


    // Step 3, CP get
    recv_mpz(csp_sock, eres);
    close(csp_sock);

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_clears(ed_1, ed_1_P1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
    mpz_clears(e_trust_2r1_prime_add_r2_prime, e_trust_r2_prime_sub_r1_prime, NULL);
    mpz_clears(ed_2, ed_2_P1, em_1r1_prime, r1_prime, neg_r1_prime, neg_em_1r1_prime, NULL);
    mpz_clears(em_3_sub_m_2, epi_1_m_3_sub_m_2, epi_2_m_3_sub_m_2, NULL);
}

#endif