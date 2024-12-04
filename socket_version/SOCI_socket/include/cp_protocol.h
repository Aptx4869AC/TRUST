/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#ifndef CP_PROTOCOL_H
#define CP_PROTOCOL_H
#pragma once

#include "paillier.h"
#include "network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

mpz_t zero, one, neg_one;

ssize_t total_bytes = 0;

void getRandom(mpz_t r, int len, gmp_randstate_t state)
{
    do
    {
        mpz_urandomb(r, state, len);
    } while (mpz_cmp_ui(r, 0) == 0);
}

void getRandomwithUpper(mpz_t r, int len, mpz_t up, gmp_randstate_t state)
{
    do
    {
        getRandom(r, len, state);
    } while (mpz_cmp_ui(r, 0) == 0 || mpz_cmp(r, up) >= 0);
}

/**
 * SOCI乘法
 * @param res
 * @param ex
 * @param ey
 * @param cp
 * @param record_bytes_flag
 */
void CP_SMUL(mpz_t res, mpz_t &ex, mpz_t &ey, PaillierThdDec &cp, bool record_bytes_flag)
{
    // Step-1 CP
    mpz_t r1, r2, er1, er2, X, Y, X1, Y1;
    mpz_inits(r1, r2, er1, er2, X, Y, X1, Y1, NULL);
    mpz_t exr2, eyr1, r1r2, er1r2, exy;
    mpz_inits(exr2, eyr1, r1r2, er1r2, exy, NULL);

    mpz_urandomb(r1, randstate, sigma); // 生成 SIGMA 位数的随机数
    mpz_urandomb(r2, randstate, sigma); // 生成 SIGMA 位数的随机数
    cp.pai.encrypt(er1, r1);
    cp.pai.encrypt(er2, r2);
    cp.pai.add(X, ex, er1);
    cp.pai.add(Y, ey, er2);
    cp.pdec(X1, X);
    cp.pdec(Y1, Y);

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);
    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 1);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);
    // CP发送 <X,Y,X1,Y1> 给CSP

    send_mpz(csp_sock, X);
    send_mpz(csp_sock, Y);
    send_mpz(csp_sock, X1);
    send_mpz(csp_sock, Y1);

    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(X);
        src_vec.push_back(Y);
        src_vec.push_back(X1);
        src_vec.push_back(Y1);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }


    // Step-3 CP
    mpz_powm(exr2, ex, r2, cp.pai.pubkey.nsquare);
    mpz_powm(eyr1, ey, r1, cp.pai.pubkey.nsquare);
    mpz_mul(r1r2, r1, r2);
    cp.pai.encrypt(er1r2, r1r2);
    cp.pai.add(res, exr2, eyr1);
    cp.pai.add(res, res, er1r2);
    cp.pai.scl_mul(res, res, neg_one);
    recv_mpz(csp_sock, exy);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(exy);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }


    cp.pai.add(res, exy, res);

    mpz_clears(r1, r2, er1, er2, X, Y, X1, Y1, NULL);
    mpz_clears(exr2, eyr1, r1r2, er1r2, exy, NULL);
}

/**
 * SOCI比较
 * @param res
 * @param ex
 * @param ey
 * @param cp
 * @param record_bytes_flag
 */
void CP_SCMP(mpz_t res, mpz_t &ex, mpz_t &ey, PaillierThdDec &cp, bool record_bytes_flag)
{

    // Step-1
    mpz_t r1, r2, r0, r1r2, er2, D, D1, exr, eyr;
    mpz_inits(r1, r2, r0, r1r2, er2, D, D1, exr, eyr, NULL);

    mpz_t mid;
    mpz_init(mid);
    mpz_set(mid, cp.pai.pubkey.half_n); // mid = n
    getRandom(r1, sigma, randstate); // r1 = random(sigma)
    getRandomwithUpper(r2, sigma, r1, randstate); // r2 = random([0, r1-1])
    mpz_sub(r2, mid, r2); // r2 = mid - r2

    mpz_add(r1r2, r1, r2);        // r1r2 = r1 + r2
    if (mpz_cmp(r2, cp.pai.pubkey.half_n) > 0)
    {
        printf("r2 <= N/2\n");
        exit(-1);
    }
    if (mpz_cmp(r1r2, cp.pai.pubkey.half_n) <= 0)
    {
        printf("r1+r2 > N/2\n");
        exit(-1);
    }
    int pi = rand() % 2;
    if (pi == 0)         // D = [r_1*(x-y+1)+r2]
    {
        mpz_add(r2, r1, r2);        // r2 = r1 + r2
        cp.pai.encrypt(er2, r2);    // er2 = [r1+r2]
        cp.pai.scl_mul(exr, ex, r1); // ex = [r1 * x]
        //mpz_neg(r1, r1);
        cp.pai.scl_mul(eyr, ey, neg_one); // ey = [-r1 * y]
        cp.pai.scl_mul(eyr, eyr, r1); // ey = [-r1 * y]
        cp.pai.add(D, exr, eyr);
        cp.pai.add(D, D, er2);


    } else                              // D = [r_1*(y-x)+r2]
    {
        cp.pai.encrypt(er2, r2);
        cp.pai.scl_mul(eyr, ey, r1);
        cp.pai.scl_mul(exr, ex, neg_one);
        cp.pai.scl_mul(exr, exr, r1);

        cp.pai.add(D, eyr, exr);
        cp.pai.add(D, D, er2);
    }
    cp.pdec(D1, D);

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 2);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);
    send_mpz(csp_sock, D1);
    send_mpz(csp_sock, D);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(D1);
        src_vec.push_back(D);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }


    // Step-3
    recv_mpz(csp_sock, res);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(res);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    if (pi == 0)
    {
        mpz_set(res, res);
    } else
    {
        cp.pai.scl_mul(res, res, neg_one);
        mpz_t eone;
        mpz_init(eone);
        cp.pai.encrypt(eone, one);
        cp.pai.add(res, eone, res);
        mpz_clear(eone);
    }

    mpz_clears(r1, r2, r0, r1r2, er2, D, D1, exr, eyr, mid, NULL);
}


/**
 * SOCI符号位获取
 * @param s_x
 * @param u_x
 * @param ex
 * @param cp
 * @param record_bytes_flag
 */
void SSBA(mpz_t s_x, mpz_t u_x, mpz_t &ex, PaillierThdDec &cp, bool record_bytes_flag)
{
    // Step-1
    mpz_t ezero;
    mpz_init(ezero);
    cp.pai.encrypt(ezero, zero);
    CP_SCMP(s_x, ex, ezero, cp, record_bytes_flag);

    // Step-2
    mpz_t sign;
    mpz_init(sign);

    mpz_t n_sub_2;
    mpz_init(n_sub_2);
    mpz_sub_ui(n_sub_2, cp.pai.pubkey.n, 2);

    cp.pai.scl_mul(sign, s_x, n_sub_2);   // [s_x]^(N-2)
    mpz_t eone;
    mpz_init(eone);
    cp.pai.encrypt(eone, one);
    cp.pai.add(sign, eone, sign);    // [1]*[s_x^(N-2)]

    // Step-3
    CP_SMUL(u_x, sign, ex, cp, record_bytes_flag);
    mpz_clear(sign);
    mpz_clears(eone, ezero, NULL);
}

/**
 * SOCI除法
 * @param eq
 * @param er
 * @param ex
 * @param ey
 * @param ell
 * @param cp
 * @param record_bytes_flag
 */
void SDIV(mpz_t eq, mpz_t &er, mpz_t &ex, mpz_t &ey, int ell, PaillierThdDec &cp, bool record_bytes_flag)
{
    // Step-1
    mpz_t ezero;
    mpz_init(ezero);
    cp.pai.encrypt(ezero, zero);
    mpz_set(eq, ezero);

    mpz_t c, u, e, ue, m, two;
    mpz_inits(c, u, e, ue, m, two, NULL);
    mpz_set_si(two, 2);
    mpz_set(er, ex);

    mpz_t eone;
    mpz_init(eone);
    cp.pai.encrypt(eone, one);

    for (int i = ell; i >= 0; i--)
    {
        // Step-2
        mpz_pow_ui(e, two, i);      // e=2^i
        cp.pai.scl_mul(c, ey, e);   // [y]^{2^i}

        // Step-3
        CP_SCMP(u, er, c, cp, record_bytes_flag);

        // Step-4
        cp.pai.scl_mul(u, u, neg_one);
        cp.pai.add(u, eone, u);
        cp.pai.scl_mul(ue, u, e);
        cp.pai.add(eq, eq, ue);

        // Step-5
        CP_SMUL(m, u, c, cp, record_bytes_flag);

        // Step-6
        cp.pai.scl_mul(m, m, neg_one);
        cp.pai.add(er, er, m);
    }
    mpz_clears(c, u, e, ue, m, two, NULL);
    mpz_clears(eone, ezero, NULL);
}

/**
 * TRUST相等判断
 * @param res
 * @param em_1
 * @param em_2
 * @param cp
 * @param record_bytes_flag
 */
void CP_Trust_FEQL(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, bool record_bytes_flag)
{
    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_t ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1;
    mpz_inits(ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);

    mpz_t half_N;
    mpz_init(half_N);
    mpz_div_ui(half_N, cp.pai.pubkey.n, 2);

    mpz_t tmp;
    mpz_init(tmp);
    mpz_urandomb(tmp, randstate, sigma);

    // e_trust_r1_add_r2、e_trust_r2
    mpz_t trust_r1, trust_r2;
    mpz_inits(trust_r1, trust_r2, NULL);
    mpz_urandomb(trust_r1, randstate, sigma);
    mpz_set(r1, trust_r1);
    mpz_neg(neg_r1, r1);

    mpz_urandomb(tmp, randstate, sigma);
    while (true)
    {
        if (mpz_cmp(tmp, trust_r1) < 0)
        {
            break;
        }
        mpz_urandomb(tmp, randstate, sigma);
    }
    mpz_sub(trust_r2, half_N, tmp);
    mpz_add(e_trust_r1_add_r2, trust_r1, trust_r2);

    // Step 1
    int pi_1 = (int) random() % 2;
    if (pi_1 == 0)
    {
        // [r1(m1-m2) + r1+r2]
        cp.pai.scl_mul(em_1r1, em_1, r1);
        cp.pai.scl_mul(neg_em_2r1, em_2, neg_r1);
        cp.pai.add(ed_1, em_1r1, neg_em_2r1);
        cp.pai.encrypt(e_trust_r1_add_r2, e_trust_r1_add_r2);
        cp.pai.add(ed_1, ed_1, e_trust_r1_add_r2);
    } else
    {
        // [r1(m2-m1) + r2]
        cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
        cp.pai.scl_mul(em_2r1, em_2, r1);
        cp.pai.add(ed_1, neg_em_1r1, em_2r1);
        cp.pai.encrypt(e_trust_r2, trust_r2);
        cp.pai.add(ed_1, ed_1, e_trust_r2);
    }
    cp.pdec(ed_1_P1, ed_1);
    mpz_t epi_1;
    mpz_inits(epi_1, NULL);
    if (pi_1 == 0)
    {
        mpz_set_si(epi_1, 0);
        cp.pai.encrypt(epi_1, epi_1);

    } else
    {
        mpz_set_si(epi_1, 1);
        cp.pai.encrypt(epi_1, epi_1);
    }

    mpz_t e_trust_r2_prime, e_trust_r1_prime_add_r2_prime;
    mpz_inits(e_trust_r2_prime, e_trust_r1_prime_add_r2_prime, NULL);
    mpz_t ed_2, ed_2_P1, em_1r1_prime, em_2r1_prime, r1_prime, neg_r1_prime, neg_em_2r1_prime, neg_em_1r1_prime;
    mpz_inits(ed_2, ed_2_P1, em_1r1_prime, em_2r1_prime, r1_prime, neg_r1_prime, neg_em_2r1_prime, neg_em_1r1_prime, NULL);


    // e_trust_r1_prime_add_r2_prime、e_trust_r2_prime
    mpz_t trust_r1_prime, trust_r2_prime;
    mpz_inits(trust_r1_prime, trust_r2_prime, NULL);
    mpz_urandomb(trust_r1_prime, randstate, sigma);
    mpz_set(r1_prime, trust_r1_prime);
    mpz_neg(neg_r1_prime, r1_prime);

    mpz_urandomb(tmp, randstate, sigma);
    while (true)
    {
        if (mpz_cmp(tmp, trust_r1_prime) < 0)
        {
            break;
        }
        mpz_urandomb(tmp, randstate, sigma);
    }
    mpz_sub(trust_r2_prime, half_N, tmp);
    mpz_add(e_trust_r1_prime_add_r2_prime, trust_r1_prime, trust_r2_prime);

    int pi_2 = (int) random() % 2;
    // Step 1
    if (pi_2 == 0)
    {
        // [r1_prime(m2-m1) + r1_prime+r2_prime]
        cp.pai.scl_mul(em_2r1_prime, em_2, r1_prime);
        cp.pai.scl_mul(neg_em_1r1_prime, em_1, neg_r1_prime);
        cp.pai.add(ed_2, neg_em_1r1_prime, em_2r1_prime);
        cp.pai.encrypt(e_trust_r1_prime_add_r2_prime, e_trust_r1_prime_add_r2_prime);
        cp.pai.add(ed_2, ed_2, e_trust_r1_prime_add_r2_prime);
    } else
    {
        // [r1_prime(m1-m2) + r2_prime]
        cp.pai.scl_mul(neg_em_2r1_prime, em_2, neg_r1_prime);
        cp.pai.scl_mul(em_1r1_prime, em_1, r1_prime);
        cp.pai.add(ed_2, em_1r1_prime, neg_em_2r1_prime);
        cp.pai.encrypt(e_trust_r2_prime, trust_r2_prime);
        cp.pai.add(ed_2, ed_2, e_trust_r2_prime);
    }
    cp.pdec(ed_2_P1, ed_2);
    mpz_t epi_2;
    mpz_inits(epi_2, NULL);
    if (pi_2 == 0)
    {
        mpz_set_si(epi_2, 0);
        cp.pai.encrypt(epi_2, epi_2);

    } else
    {
        mpz_set_si(epi_2, 1);
        cp.pai.encrypt(epi_2, epi_2);
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
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(ed_1);
        src_vec.push_back(ed_1_P1);
        src_vec.push_back(epi_1);
        src_vec.push_back(ed_2);
        src_vec.push_back(ed_2_P1);
        src_vec.push_back(epi_2);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    // Step 3, CP get
    recv_mpz(csp_sock, res);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(res);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_clears(ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);
    mpz_clears(epi_1, NULL);
    mpz_clears(e_trust_r2_prime, e_trust_r1_prime_add_r2_prime, NULL);
    mpz_clears(ed_2, ed_2_P1, em_1r1_prime, em_2r1_prime, r1_prime, neg_r1_prime, neg_em_2r1_prime, neg_em_1r1_prime, NULL);
    mpz_clears(epi_2, NULL);
}


/**
 * TRUST绝对值判断
 * @param eres
 * @param em_1
 * @param cp
 * @param record_bytes_flag
 */
void CP_Trust_FABS(mpz_t eres, mpz_t em_1, PaillierThdDec cp, bool record_bytes_flag)
{

    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);

    mpz_t ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1;
    mpz_inits(ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);


    mpz_t half_N;
    mpz_init(half_N);
    mpz_div_ui(half_N, cp.pai.pubkey.n, 2);

    mpz_t tmp;
    mpz_init(tmp);
    mpz_urandomb(tmp, randstate, sigma);

    // e_trust_r1_add_r2、e_trust_r2
    mpz_t trust_r1, trust_r2;
    mpz_inits(trust_r1, trust_r2, NULL);
    mpz_urandomb(trust_r1, randstate, sigma);
    mpz_set(r1, trust_r1);
    mpz_neg(neg_r1, r1);

    mpz_urandomb(tmp, randstate, sigma);
    while (true)
    {
        if (mpz_cmp(tmp, trust_r1) < 0)
        {
            break;
        }
        mpz_urandomb(tmp, randstate, sigma);
    }
    mpz_sub(trust_r2, half_N, tmp);
    mpz_add(e_trust_r1_add_r2, trust_r1, trust_r2);

    int pi = (int) random() % 2;
    // Step 1
    if (pi == 0)
    {
        cp.pai.scl_mul(em_1r1, em_1, r1);
        cp.pai.encrypt(e_trust_r1_add_r2, e_trust_r1_add_r2);
        cp.pai.add(ed, em_1r1, e_trust_r1_add_r2);
    } else
    {
        cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
        cp.pai.encrypt(e_trust_r2, trust_r2);
        cp.pai.add(ed, neg_em_1r1, e_trust_r2);
    }
    cp.pdec(ed_1, ed);

    // ezero
    mpz_t ezero;
    mpz_init(ezero);
    cp.pai.encrypt(ezero, zero);

    mpz_t em_1pi;
    mpz_inits(em_1pi, NULL);
    cp.pai.scl_mul(em_1pi, em_1, pi);
    cp.pai.add(em_1pi, em_1pi, ezero);


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
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(ed);
        src_vec.push_back(ed_1);
        src_vec.push_back(em_1);
        src_vec.push_back(em_1pi);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    // Step 3, CP get
    recv_mpz(csp_sock, eres);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(eres);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

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
 * @param record_bytes_flag
 */
void CP_Trust_FTRN_version1(mpz_t eres, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, bool record_bytes_flag)
{

    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_t ed_1, ed_1_P1, em_1r1, r1, neg_r1, neg_em_1r1;
    mpz_inits(ed_1, ed_1_P1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
    mpz_urandomb(r1, randstate, sigma);
    mpz_neg(neg_r1, r1);


    mpz_t half_N;
    mpz_init(half_N);
    mpz_div_ui(half_N, cp.pai.pubkey.n, 2);

    mpz_t tmp;
    mpz_init(tmp);
    mpz_urandomb(tmp, randstate, sigma);

    // e_trust_r1_add_r2、e_trust_r2
    mpz_t trust_r1, trust_r2;
    mpz_inits(trust_r1, trust_r2, NULL);
    mpz_urandomb(trust_r1, randstate, sigma);
    mpz_set(r1, trust_r1);
    mpz_neg(neg_r1, r1);

    mpz_urandomb(tmp, randstate, sigma);
    while (true)
    {
        if (mpz_cmp(tmp, trust_r1) < 0)
        {
            break;
        }
        mpz_urandomb(tmp, randstate, sigma);
    }
    mpz_sub(trust_r2, half_N, tmp);
    mpz_add(e_trust_r1_add_r2, trust_r1, trust_r2);



    // Step 1
    int pi_1 = rand() % 2;
    if (pi_1 == 0)
    {
        // [r1*m1 + r2]
        cp.pai.scl_mul(em_1r1, em_1, r1);
        cp.pai.encrypt(e_trust_r2, trust_r2);
        cp.pai.add(ed_1, em_1r1, e_trust_r2);
    } else
    {
        // [-r1*m1 + r1+r2]
        cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
        cp.pai.encrypt(e_trust_r1_add_r2, e_trust_r1_add_r2);
        cp.pai.add(ed_1, neg_em_1r1, e_trust_r1_add_r2);
    }
    cp.pdec(ed_1_P1, ed_1);
    mpz_t e_trust_2r1_prime_add_r2_prime, e_trust_r2_prime_sub_r1_prime;
    mpz_inits(e_trust_2r1_prime_add_r2_prime, e_trust_r2_prime_sub_r1_prime, NULL);

    mpz_t ed_2, ed_2_P1, em_1r1_prime, r1_prime, neg_r1_prime, neg_em_1r1_prime;
    mpz_inits(ed_2, ed_2_P1, em_1r1_prime, r1_prime, neg_r1_prime, neg_em_1r1_prime, NULL);
    mpz_urandomb(r1_prime, randstate, sigma);
    mpz_neg(neg_r1_prime, r1_prime);


    // e_trust_2r1_prime_add_r2_prime、e_trust_r2_prime_sub_r1_prime
    mpz_t trust_r1_prime, trust_r2_prime;
    mpz_inits(trust_r1_prime, trust_r2_prime, NULL);
    mpz_urandomb(trust_r1_prime, randstate, sigma);
    mpz_set(r1_prime, trust_r1_prime);
    mpz_neg(neg_r1_prime, r1_prime);

    mpz_urandomb(tmp, randstate, sigma);
    while (true)
    {
        if (mpz_cmp(tmp, trust_r1_prime) < 0)
        {
            break;
        }
        mpz_urandomb(tmp, randstate, sigma);
    }
    mpz_sub(trust_r2_prime, half_N, tmp);

    mpz_add(e_trust_2r1_prime_add_r2_prime, trust_r1_prime, trust_r1_prime);
    mpz_add(e_trust_2r1_prime_add_r2_prime, e_trust_2r1_prime_add_r2_prime, trust_r2_prime);
    mpz_sub(e_trust_r2_prime_sub_r1_prime, trust_r2_prime, trust_r1_prime);


    int pi_2 = rand() % 2;
    if (pi_2 == 0)
    {
        cp.pai.scl_mul(neg_em_1r1_prime, em_1, neg_r1_prime);
        cp.pai.encrypt(e_trust_2r1_prime_add_r2_prime, e_trust_2r1_prime_add_r2_prime);
        cp.pai.add(ed_2, neg_em_1r1_prime, e_trust_2r1_prime_add_r2_prime);
    } else
    {

        cp.pai.scl_mul(em_1r1_prime, em_1, r1_prime);
        cp.pai.encrypt(e_trust_r2_prime_sub_r1_prime, e_trust_r2_prime_sub_r1_prime);
        cp.pai.add(ed_2, em_1r1_prime, e_trust_r2_prime_sub_r1_prime);
    }

    // ezero
    mpz_t ezero;
    mpz_init(ezero);
    cp.pai.encrypt(ezero, zero);

    cp.pdec(ed_2_P1, ed_2);
    mpz_t em_3_sub_m_2, epi_1_m_3_sub_m_2, epi_2_m_3_sub_m_2;
    mpz_inits(em_3_sub_m_2, epi_1_m_3_sub_m_2, epi_2_m_3_sub_m_2, NULL);

    // [m3-m2]
    cp.pai.scl_mul(em_3_sub_m_2, em_2, -1);
    cp.pai.add(em_3_sub_m_2, em_3, em_3_sub_m_2);
    // [pi_1 * (m3-m2)]
    cp.pai.scl_mul(epi_1_m_3_sub_m_2, em_3_sub_m_2, pi_1);
    cp.pai.add(epi_1_m_3_sub_m_2, epi_1_m_3_sub_m_2, ezero);

    // [pi_2 * (m3-m2)]
    cp.pai.scl_mul(epi_2_m_3_sub_m_2, em_3_sub_m_2, pi_2);
    cp.pai.add(epi_2_m_3_sub_m_2, epi_2_m_3_sub_m_2, ezero);

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

    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(em_3_sub_m_2);
        src_vec.push_back(em_2);
        src_vec.push_back(ed_1);
        src_vec.push_back(ed_1_P1);
        src_vec.push_back(epi_1_m_3_sub_m_2);
        src_vec.push_back(ed_2);
        src_vec.push_back(ed_2_P1);
        src_vec.push_back(epi_2_m_3_sub_m_2);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    // Step 3, CP get
    recv_mpz(csp_sock, eres);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(eres);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_clears(ed_1, ed_1_P1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
    mpz_clears(e_trust_2r1_prime_add_r2_prime, e_trust_r2_prime_sub_r1_prime, NULL);
    mpz_clears(ed_2, ed_2_P1, em_1r1_prime, r1_prime, neg_r1_prime, neg_em_1r1_prime, NULL);
    mpz_clears(em_3_sub_m_2, epi_1_m_3_sub_m_2, epi_2_m_3_sub_m_2, NULL);
}


/**
 * TRUST三目运算
 * @param eres
 * @param em_1
 * @param em_2
 * @param em_3
 * @param cp
 * @param record_bytes_flag
 */
void CP_Trust_FTRN_version2(mpz_t eres, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, bool record_bytes_flag)
{
    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);

    mpz_t ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1;
    mpz_inits(ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
    mpz_urandomb(r1, randstate, sigma);
    mpz_neg(neg_r1, r1);

    mpz_t half_N;
    mpz_init(half_N);
    mpz_div_ui(half_N, cp.pai.pubkey.n, 2);

    mpz_t tmp;
    mpz_init(tmp);
    mpz_urandomb(tmp, randstate, sigma);

    // e_trust_r1_add_r2、e_trust_r2
    mpz_t trust_r1, trust_r2;
    mpz_inits(trust_r1, trust_r2, NULL);
    mpz_urandomb(trust_r1, randstate, sigma);
    mpz_set(r1, trust_r1);
    mpz_neg(neg_r1, r1);

    mpz_urandomb(tmp, randstate, sigma);
    while (true)
    {
        if (mpz_cmp(tmp, trust_r1) < 0)
        {
            break;
        }
        mpz_urandomb(tmp, randstate, sigma);
    }
    mpz_sub(trust_r2, half_N, tmp);
    mpz_add(e_trust_r1_add_r2, trust_r1, trust_r2);

    int pi = rand() % 2;
    if (pi == 0)
    {
        // [r1*m1 + r2]
        cp.pai.scl_mul(em_1r1, em_1, r1);
        cp.pai.encrypt(e_trust_r2, trust_r2);
        cp.pai.add(ed, em_1r1, e_trust_r2);
    } else
    {
        // [-r1*m1 + r1+r2]
        cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
        cp.pai.encrypt(e_trust_r1_add_r2, e_trust_r1_add_r2);
        cp.pai.add(ed, neg_em_1r1, e_trust_r1_add_r2);
    }
    cp.pdec(ed_1, ed);

    mpz_t em_3_sub_m_2, epi_m_3_sub_m_2;
    mpz_inits(em_3_sub_m_2, epi_m_3_sub_m_2, NULL);

    // [m3-m2]
    cp.pai.scl_mul(em_3_sub_m_2, em_2, -1);
    cp.pai.add(em_3_sub_m_2, em_3, em_3_sub_m_2);
    // [pi * (m3-m2)]
    cp.pai.scl_mul(epi_m_3_sub_m_2, em_3_sub_m_2, pi);

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 11);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    /* send <[m3], [m2], [d_1], (PDec)[d_1], [pi * (m3-m2)] > to csp */
    send_mpz(csp_sock, ed);
    send_mpz(csp_sock, ed_1);
    send_mpz(csp_sock, em_2);
    send_mpz(csp_sock, em_3);
    send_mpz(csp_sock, epi_m_3_sub_m_2);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(ed);
        src_vec.push_back(ed_1);
        src_vec.push_back(em_2);
        src_vec.push_back(em_3);
        src_vec.push_back(epi_m_3_sub_m_2);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    // Step 3, CP get
    recv_mpz(csp_sock, eres);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(eres);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }


    mpz_clears(e_trust_r2, e_trust_r1_add_r2, NULL);
    mpz_clears(ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
    mpz_clears(em_3_sub_m_2, epi_m_3_sub_m_2, NULL);
}

/**
 * 乘法正确性检验
 * @param i
 * @param result_cxy
 * @param x
 * @param y
 * @param pai
 */
void check_smul(int i, mpz_t result_cxy, mpz_t x, mpz_t y, Paillier pai)
{
    mpz_t z, xy;
    mpz_inits(z, xy, NULL);
    pai.decrypt(z, result_cxy);
    mpz_mul(xy, x, y);
    if (mpz_cmp(z, xy) != 0)
    {
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
void check_scmp(int i, mpz_t result_cmp_x_y, mpz_t x, mpz_t y, Paillier pai)
{
    mpz_t z, answer;
    mpz_inits(z, answer, NULL);
    if (mpz_cmp(x, y) >= 0)
    {
        mpz_set_si(answer, 0);
    } else mpz_set_si(answer, 1);

    pai.decrypt(z, result_cmp_x_y);
    if (mpz_cmp(z, answer) != 0)
    {
        gmp_printf("dec answer x>=y? = %Zd\n", z);
        gmp_printf("right answer %Zd\n", answer);
        printf("SCMP error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(z, answer, NULL);
}

/**
 * 符号位获取正确性检验
 * @param i
 * @param result_cs
 * @param result_cu
 * @param x
 * @param pai
 */
void check_ssba(int i, mpz_t result_cs, mpz_t result_cu, mpz_t x, Paillier pai)
{
    mpz_t s, u, answer_s, answer_u;
    mpz_inits(s, u, answer_s, answer_u, NULL);

    pai.decrypt(s, result_cs);
    pai.decrypt(u, result_cu);

    if (mpz_cmp_si(x, 0) >= 0)
    {
        mpz_set_si(answer_s, 0);
        mpz_set(answer_u, x);
    } else
    {
        mpz_set_si(answer_s, 1);
        mpz_set(answer_u, x);
        mpz_neg(answer_u, answer_u); // answer_s = -answer_s
    }
    if (mpz_cmp(answer_u, u) != 0 || mpz_cmp(answer_s, s) != 0)
    {
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
void check_sdiv(int i, mpz_t result_cq, mpz_t result_ce, mpz_t x, mpz_t y, Paillier pai)
{
    mpz_t q, e, right_q, right_e;
    mpz_inits(q, e, right_q, right_e, NULL);
    pai.decrypt(q, result_cq);
    pai.decrypt(e, result_ce);

    mpz_div(right_q, x, y);
    mpz_mod(right_e, x, y);
    if (mpz_cmp(q, right_q) != 0 && mpz_cmp(e, right_e) != 0)
    {
        gmp_printf("q = %Zd e = %Zd\n", q, e);
        gmp_printf("right_q = %Zd right_e = %Zd\n", right_q, right_e);
        printf("SDIV error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(q, e, right_q, right_e, NULL);
}

/**
 * 检查Trust_FEQL正确性
 * @param i
 * @param eres
 * @param em_1
 * @param em_2
 * @param pai
 */
void check_feql(int i, mpz_t eres, mpz_t m_1, mpz_t m_2, Paillier pai)
{
    mpz_t half_N;
    mpz_init(half_N);
    mpz_div_ui(half_N, pai.pubkey.n, 2);
    if (mpz_cmp(m_1, half_N) > 0)
    {
        mpz_sub(m_1, m_1, pai.pubkey.n);
    }
    if (mpz_cmp(m_2, half_N) > 0)
    {
        mpz_sub(m_2, m_2, pai.pubkey.n);
    }


    mpz_t res, m_1_eql_m_2;
    mpz_inits(res, m_1_eql_m_2, NULL);
    if (mpz_cmp(m_1, m_2) == 0)
    {
        mpz_set_si(m_1_eql_m_2, 0);
    } else mpz_set_si(m_1_eql_m_2, 1);

    pai.decrypt(res, eres);
    if (mpz_cmp(res, m_1_eql_m_2) != 0)
    {
        gmp_printf("m_1 = %Zd\n", m_1);
        gmp_printf("m_2 = %Zd\n", m_2);
        gmp_printf("dec answer x==y? = %Zd\n", res);
        gmp_printf("right answer %Zd\n", m_1_eql_m_2);
        printf("FEQL error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(res, m_1_eql_m_2, NULL);
}

/**
 * 检查Trust_FABS正确性
 * @param i
 * @param eres
 * @param em_1
 * @param pai
 */
void check_fabs(int i, mpz_t eres, mpz_t m_1, Paillier pai)
{
    mpz_t half_N;
    mpz_init(half_N);
    mpz_div_ui(half_N, pai.pubkey.n, 2);
    if (mpz_cmp(m_1, half_N) > 0)
    {
        mpz_sub(m_1, m_1, pai.pubkey.n);
    }

    mpz_t res, m_1_abs;
    mpz_inits(res, m_1_abs, NULL);
    if (mpz_cmp_si(m_1, 0) > 0)
    {
        mpz_set(m_1_abs, m_1);
    } else
    {
        mpz_set(m_1_abs, m_1);
        mpz_neg(m_1_abs, m_1_abs);
    }

    pai.decrypt(res, eres);
    if (mpz_cmp(res, m_1_abs) != 0)
    {
        gmp_printf("m_1 = %Zd\n", m_1);
        gmp_printf("dec answer x = %Zd\n", res);
        gmp_printf("right answer %Zd\n", m_1_abs);
        printf("FABS error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(res, m_1_abs, NULL);
}

/**
 * 检查Trust_FTRN正确性
 * @param i
 * @param eres
 * @param em_1
 * @param em_2
 * @param em_3
 * @param pai
 */
void check_ftrn(int i, mpz_t eres, mpz_t m_1, mpz_t m_2, mpz_t m_3, Paillier pai)
{
    mpz_t half_N;
    mpz_init(half_N);
    mpz_div_ui(half_N, pai.pubkey.n, 2);
    if (mpz_cmp(m_1, half_N) > 0)
    {
        mpz_sub(m_1, m_1, pai.pubkey.n);
    }
    if (mpz_cmp(m_2, half_N) > 0)
    {
        mpz_sub(m_2, m_2, pai.pubkey.n);
    }


    mpz_t res, m_1_trn_m_2_or_m_3;
    mpz_inits(res, m_1_trn_m_2_or_m_3, NULL);
    if (mpz_cmp_si(m_1, 1) == 0)
    {
        mpz_set(m_1_trn_m_2_or_m_3, m_2);
    } else
    {
        mpz_set(m_1_trn_m_2_or_m_3, m_3);
    }

    pai.decrypt(res, eres);
    if (mpz_cmp(res, m_1_trn_m_2_or_m_3) != 0)
    {
        gmp_printf("m_1 = %Zd\n", m_1);
        gmp_printf("m_2 = %Zd\n", m_2);
        gmp_printf("m_3 = %Zd\n", m_3);
        gmp_printf("dec answer = %Zd\n", res);
        gmp_printf("right answer = %Zd\n", m_1_trn_m_2_or_m_3);
        printf("FTRN error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(res, m_1_trn_m_2_or_m_3, NULL);
}

#endif