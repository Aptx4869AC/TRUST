/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-31
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#ifndef PROTOCOL_H
#define PROTOCOL_H

#pragma once

#include "fastPai.h"


using namespace PHESPACE;
using namespace std;

namespace PROTOCOLSPACE
{
    class protocol
    {
    public:

        protocol();

        ~protocol();

        /** 密文乘法 **/
        void trust_fmul(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp);

        /** 密文比较 **/
        void trust_fcmp(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp);

        /** 密文相等判断 **/
        void trust_feql(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp);

        /** 密文绝对值 **/
        void trust_fabs(mpz_t res, mpz_t em_1, PaillierThdDec cp, PaillierThdDec csp);

        /** 密文三目运算，论文版本 **/
        void trust_ftrn_version1(mpz_t res, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, PaillierThdDec csp);

        /** 密文三目运算，更高效方案，非论文版本 **/
        void trust_ftrn_version2(mpz_t res, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, PaillierThdDec csp);

    };

    /**
     * 密文乘法
     * @param res  密文结果，即 [m1*m2]
     * @param em_1 密文输入m1，记[m1]
     * @param em_2 密文输入m2，记[m2]
     * @param cp   服务器cp
     * @param csp  服务器csp
     */
    void protocol::trust_fmul(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp)
    {

        Offline_val tuple_S0 = cp.pai.offlineVal;

        // Step 1, CP computes
        mpz_t r, neg_r, er;
        mpz_inits(r, neg_r, er, NULL);
        mpz_set(r, tuple_S0.trust_r);
        cp.pai.add(er, tuple_S0.e_trust_r, tuple_S0.e0); //refresh

        mpz_t X, Y, X_PDec_1;
        mpz_inits(X, Y, X_PDec_1, NULL);
        // X = [m_1 + r]
        cp.pai.add(X, em_1, er);

        // Y = [-r * m_2]
        mpz_neg(neg_r, r);
        cp.pai.scl_mul(Y, em_2, neg_r);

        // PDec(X)
        cp.PDec(X_PDec_1, X);

        /* send <[m_1 + r], (PDec)[m_1 + r], [m_2], [-r * m_2]> to csp */

        // Step 2, CSP computes
        mpz_t X_PDec_2, m_1_add_r;
        mpz_inits(X_PDec_2, m_1_add_r, NULL);
        // PDec + TDec
        csp.PDec(X_PDec_2, X);
        csp.TDec(m_1_add_r, X_PDec_1, X_PDec_2);

        csp.pai.scl_mul(res, em_2, m_1_add_r);
        csp.pai.add(res, res, Y);
        csp.pai.add(res, res, tuple_S0.e0); //refresh

        /* send [R] to cp */

        mpz_clears(r, neg_r, er, X, Y, X_PDec_1, X_PDec_2, m_1_add_r, NULL);

    }

    /**
     * 密文比较
     * @param res  密文结果，能判定[m1]与[m2]大小，当m1>=m2，结果为0，否则结果为1
     * @param em_1 密文输入m1，记[m1]
     * @param em_2 密文输入m2，记[m2]
     * @param cp   服务器cp
     * @param csp  服务器csp
     */
    void protocol::trust_fcmp(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;
        Offline_val tuple_S1 = csp.pai.offlineVal;

        mpz_t e_trust_r2, e_trust_r1_add_r2;
        mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);

        mpz_t ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1;
        mpz_inits(ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);
        mpz_set(r1, tuple_S0.trust_r1);
        mpz_neg(neg_r1, r1);

        int pi = (int) random() % 2;
        // Step 1
        if (pi == 0)
        {
            cp.pai.scl_mul(em_1r1, em_1, r1);
            cp.pai.scl_mul(neg_em_2r1, em_2, neg_r1);
            cp.pai.add(ed, em_1r1, neg_em_2r1);
            cp.pai.add(e_trust_r1_add_r2, tuple_S0.e_trust_r1_add_r2, tuple_S0.e0);  //refresh
            cp.pai.add(ed, ed, e_trust_r1_add_r2);
        } else
        {
            cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
            cp.pai.scl_mul(em_2r1, em_2, r1);
            cp.pai.add(ed, neg_em_1r1, em_2r1);
            cp.pai.add(e_trust_r2, tuple_S0.e_trust_r2, tuple_S0.e0);  //refresh
            cp.pai.add(ed, ed, e_trust_r2);
        }
        cp.PDec(ed_1, ed);

        mpz_t epi;
        mpz_inits(epi, NULL);
        if (pi == 0)
        {
            cp.pai.add(epi, tuple_S0.e0, tuple_S0.e0);   //refresh

        } else
        {
            cp.pai.add(epi, tuple_S0.e1, tuple_S0.e0);   //refresh
        }

        /* send <[d], (PDec)[d], [pai]> to csp */

        // Step 2
        mpz_t ed_2, d, miu0, emiu0, one_sub_2miu0, half_N;
        mpz_inits(ed_2, d, miu0, emiu0, one_sub_2miu0, half_N, NULL);
        mpz_div_ui(half_N, csp.pk.N, 2);
        csp.PDec(ed_2, ed);
        csp.TDec(d, ed_1, ed_2);
        if (mpz_cmp(d, half_N) > 0)
        {
            mpz_set_si(miu0, 0);
            csp.pai.add(emiu0, tuple_S1.e0, tuple_S1.e0);   //refresh
        } else if (mpz_cmp(d, half_N) < 0)
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
        csp.pai.scl_mul(res, epi, one_sub_2miu0);
        csp.pai.add(res, emiu0, res);

        /* send [R] to cp */

        mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1,
                   epi, ed_2, d, miu0, emiu0, one_sub_2miu0, half_N, NULL);
    }

    /**
     * 密文判断相等
     * @param res  密文结果，能判定[m1]与[m2]是否相等，当m1=m2，结果为0，否则结果为1
     * @param em_1 密文输入m1，记[m1]
     * @param em_2 密文输入m2，记[m2]
     * @param cp   服务器cp
     * @param csp  服务器csp
     */
    void protocol::trust_feql(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;
        Offline_val tuple_S1 = csp.pai.offlineVal;

        mpz_t e_trust_r2, e_trust_r1_add_r2;
        mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);
        mpz_t ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1;
        mpz_inits(ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);
        mpz_set(r1, tuple_S0.trust_r1);
        mpz_neg(neg_r1, r1);

        // Step 1
        int pi_1 = (int) random() % 2;
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

        int pi_2 = (int) random() % 2;
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


        /* send <[d_1], (PDec)[d_1], [pi_1], [d_2], (PDec)[d_2], [pi_2]> to csp */

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

        mpz_t left, right;
        mpz_inits(left, right, NULL);
        csp.pai.scl_mul(left, epi_1, one_sub_2miu0);
        csp.pai.add(left, emiu0, left);
        csp.pai.scl_mul(right, epi_2, one_sub_2miu0_prime);
        csp.pai.add(right, emiu0_prime, right);
        csp.pai.add(res, left, right);

        /* send [R] to cp */
        mpz_clears(e_trust_r2, e_trust_r1_add_r2, e_trust_r2_prime, e_trust_r1_prime_add_r2_prime, NULL);
        mpz_clears(ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);
        mpz_clears(ed_2, ed_2_P1, em_1r1_prime, em_2r1_prime, r1_prime, neg_r1_prime, neg_em_2r1_prime, neg_em_1r1_prime, NULL);
        mpz_clears(epi_1, epi_2, NULL);
        mpz_clears(ed_1_P2, ed_2_P2, d1, d2, miu0, emiu0, miu0_prime, emiu0_prime, one_sub_2miu0, one_sub_2miu0_prime, half_N, NULL);
        mpz_clears(left, right, NULL);
    }

    /**
     * 密文绝对值
     * @param res  密文结果，输出绝对值[|m1|]
     * @param em_1 密文输入m1，记[m1]
     * @param cp   服务器cp
     * @param csp  服务器csp
     */
    void protocol::trust_fabs(mpz_t res, mpz_t em_1, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;
        Offline_val tuple_S1 = csp.pai.offlineVal;

        mpz_t e_trust_r2, e_trust_r1_add_r2;
        mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);

        mpz_t ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1;
        mpz_inits(ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
        mpz_set(r1, tuple_S0.trust_r1);
        mpz_neg(neg_r1, r1);

        int pi = (int) random() % 2;
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

        /* send <[d], (PDec)[d], [pai * m_1], [m_1]> to csp */

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


        mpz_t left, right;
        mpz_inits(left, right, NULL);
        csp.pai.scl_mul(left, em_1, one_sub_2miu0);
        csp.pai.scl_mul(right, em_1pi, neg_two_sub_4miu0);
        csp.pai.add(res, left, right);
        csp.pai.add(res, res, tuple_S1.e0);   //refresh

        /* send [R] to cp */

        mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
        mpz_clears(em_1pi, ed_2, d, miu0, one_sub_2miu0, neg_two_sub_4miu0, half_N, left, right, NULL);
    }

    /**
     * 密文三目运算
     * @param res  密文结果，m1任意取值，当m1=1时，输出[m2]否则输出[m3]
     * @param em_1 密文输入m1，记[m1]
     * @param em_2 密文输入m2，记[m2]
     * @param em_3 密文输入m3，记[m3]
     * @param cp   服务器cp
     * @param csp  服务器csp
     */
    void protocol::trust_ftrn_version1(mpz_t res, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;
        Offline_val tuple_S1 = csp.pai.offlineVal;

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

        /* send <[m3-m2], [m2], [d_1], (PDec)[d_1], [pi_1 * (m3-m2)], [d_2], (PDec)[d_2], [pi_2 * (m3-m2)]> to csp */

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
        mpz_t left, mid, right;
        mpz_inits(left, mid, right, NULL);
        csp.pai.scl_mul(left, em_3_sub_m_2, miu0_add_miu0_prime); // left
        csp.pai.add(left, em_2, left); // left
        csp.pai.scl_mul(mid, epi_1_m_3_sub_m_2, one_sub_2miu0); // mid
        csp.pai.scl_mul(right, epi_2_m_3_sub_m_2, one_sub_2miu0_prime); // right

        csp.pai.add(res, left, mid);
        csp.pai.add(res, res, right);
        csp.pai.add(res, res, tuple_S1.e0);   //refresh

        /* send [R] to cp */
        mpz_clears(ed_1_P2, ed_2_P2, d1, d2, miu0, miu0_prime, one_sub_2miu0, one_sub_2miu0_prime, miu0_add_miu0_prime, half_N, NULL);
        mpz_clears(left, mid, right, NULL);
    }

    /**
     * 密文三目运算
     * @param res  密文结果，m1只在1与0之间取值，当m1=1时，输出[m2]否则输出[m3]
     * @param em_1 密文输入m1，记[m1]
     * @param em_2 密文输入m2，记[m2]
     * @param em_3 密文输入m3，记[m3]
     * @param cp   服务器cp
     * @param csp  服务器csp
     */
    void protocol::trust_ftrn_version2(mpz_t res, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, PaillierThdDec csp)
    {

        Offline_val tuple_S0 = cp.pai.offlineVal;
        Offline_val tuple_S1 = csp.pai.offlineVal;

        mpz_t e_trust_r2, e_trust_r1_add_r2;
        mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);

        mpz_t ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1;
        mpz_inits(ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
        mpz_set(r1, tuple_S0.trust_r1);
        mpz_neg(neg_r1, r1);

        int pi = rand() % 2;
        if (pi == 0)
        {
            // [r1*m1 + r2]
            cp.pai.scl_mul(em_1r1, em_1, r1);
            cp.pai.add(e_trust_r2, tuple_S0.e_trust_r2, tuple_S0.e0);  //refresh
            cp.pai.add(ed, em_1r1, e_trust_r2);
        } else
        {
            // [-r1*m1 + r1+r2]
            cp.pai.scl_mul(neg_em_1r1, em_1, neg_r1);
            cp.pai.add(e_trust_r1_add_r2, tuple_S0.e_trust_r1_add_r2, tuple_S0.e0);  //refresh
            cp.pai.add(ed, neg_em_1r1, e_trust_r1_add_r2);
        }
        cp.PDec(ed_1, ed);

        mpz_t em_3_sub_m_2, epi_m_3_sub_m_2;
        mpz_inits(em_3_sub_m_2, epi_m_3_sub_m_2, NULL);

        // [m3-m2]
        cp.pai.scl_mul(em_3_sub_m_2, em_2, -1);
        cp.pai.add(em_3_sub_m_2, em_3, em_3_sub_m_2);
        // [pi * (m3-m2)]
        cp.pai.scl_mul(epi_m_3_sub_m_2, em_3_sub_m_2, pi);

        /* send <[m3], [m2], [d_1], (PDec)[d_1], [pi * (m3-m2)] > to csp */

        // Step 2
        mpz_t ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N;
        mpz_inits(ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N, NULL);
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

        mpz_neg(one_sub_miu0, miu0);
        mpz_add_ui(one_sub_miu0, one_sub_miu0, 1);

        mpz_mul_si(one_sub_2miu0, miu0, 2);
        mpz_neg(one_sub_2miu0, one_sub_2miu0);
        mpz_add_ui(one_sub_2miu0, one_sub_2miu0, 1);

        // [R] compute
        mpz_t left, mid, right;
        mpz_inits(left, mid, right, NULL);
        csp.pai.scl_mul(left, em_2, one_sub_miu0); // left
        csp.pai.scl_mul(mid, em_3, miu0); // mid
        csp.pai.scl_mul(right, epi_m_3_sub_m_2, one_sub_2miu0); // right

        csp.pai.add(res, left, mid);
        csp.pai.add(res, res, right);
        csp.pai.add(res, res, tuple_S1.e0);   //refresh

        /* send [R] to cp */
        mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, NULL);
        mpz_clears(em_3_sub_m_2, epi_m_3_sub_m_2, NULL);
        mpz_clears(ed_2, d, miu0, one_sub_miu0, one_sub_2miu0, half_N, NULL);
        mpz_clears(left, mid, right, NULL);
    }


    protocol::protocol()
    {}

    protocol::~protocol()
    {}

    /**
     * 检查Trust_FMUL正确性
     * @param i
     * @param eres
     * @param em_1
     * @param em_2
     * @param pai
     */
    void check_fmul(int i, mpz_t eres, mpz_t m_1, mpz_t m_2, Paillier pai)
    {
        mpz_t half_N;
        mpz_init(half_N);
        mpz_div_ui(half_N, pai.pk.N, 2);
        if (mpz_cmp(m_1, half_N) > 0)
        {
            mpz_sub(m_1, m_1, pai.pk.N);
        }
        if (mpz_cmp(m_2, half_N) > 0)
        {
            mpz_sub(m_2, m_2, pai.pk.N);
        }

        mpz_t res, m_1_mul_m_2;
        mpz_inits(res, m_1_mul_m_2, NULL);
        pai.decrypt(res, eres);
        mpz_mul(m_1_mul_m_2, m_1, m_2);
        if (mpz_cmp(res, m_1_mul_m_2) != 0)
        {
            gmp_printf("m_1 = %Zd\n", m_1);
            gmp_printf("m_2 = %Zd\n", m_2);
            gmp_printf("dec answer = %Zd\n", res);
            gmp_printf("right answer = %Zd\n", m_1_mul_m_2);
            printf("FMUL error!\n");
            printf("i = %d is error \n", i);
            exit(-1);
        }
        mpz_clears(res, m_1_mul_m_2, NULL);
    }

    /**
     * 检查Trust_FCMP正确性
     * @param i
     * @param eres
     * @param em_1
     * @param em_2
     * @param pai
     */
    void check_fcmp(int i, mpz_t eres, mpz_t m_1, mpz_t m_2, Paillier pai)
    {
        mpz_t half_N;
        mpz_init(half_N);
        mpz_div_ui(half_N, pai.pk.N, 2);
        if (mpz_cmp(m_1, half_N) > 0)
        {
            mpz_sub(m_1, m_1, pai.pk.N);
        }
        if (mpz_cmp(m_2, half_N) > 0)
        {
            mpz_sub(m_2, m_2, pai.pk.N);
        }


        mpz_t res, m_1_cmp_m_2;
        mpz_inits(res, m_1_cmp_m_2, NULL);
        if (mpz_cmp(m_1, m_2) >= 0)
        {
            mpz_set_si(m_1_cmp_m_2, 0);
        } else mpz_set_si(m_1_cmp_m_2, 1);

        pai.decrypt(res, eres);
        if (mpz_cmp(res, m_1_cmp_m_2) != 0)
        {
            gmp_printf("m_1 = %Zd\n", m_1);
            gmp_printf("m_2 = %Zd\n", m_2);
            gmp_printf("dec answer = %Zd\n", res);
            gmp_printf("right answer = %Zd\n", m_1_cmp_m_2);
            printf("FCMP error!\n");
            printf("i = %d is error \n", i);
            exit(-1);
        }
        mpz_clears(res, m_1_cmp_m_2, NULL);
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
        mpz_div_ui(half_N, pai.pk.N, 2);
        if (mpz_cmp(m_1, half_N) > 0)
        {
            mpz_sub(m_1, m_1, pai.pk.N);
        }
        if (mpz_cmp(m_2, half_N) > 0)
        {
            mpz_sub(m_2, m_2, pai.pk.N);
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
        mpz_div_ui(half_N, pai.pk.N, 2);
        if (mpz_cmp(m_1, half_N) > 0)
        {
            mpz_sub(m_1, m_1, pai.pk.N);
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
        mpz_div_ui(half_N, pai.pk.N, 2);
        if (mpz_cmp(m_1, half_N) > 0)
        {
            mpz_sub(m_1, m_1, pai.pk.N);
        }
        if (mpz_cmp(m_2, half_N) > 0)
        {
            mpz_sub(m_2, m_2, pai.pk.N);
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

}


#endif //PROTOCOL_H

