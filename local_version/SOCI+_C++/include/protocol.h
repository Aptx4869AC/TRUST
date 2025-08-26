/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-01
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

        /**
         * Weiquan SOCI+ verison 密文乘法
         */
        void smul(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * Weiquan SOCI+ verison 密文比较
         */
        void scmp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * Weiquan SOCI+ verison 密文符号位获取
         */
        void ssba(mpz_t s_x, mpz_t u_x, mpz_t ex, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * Weiquan SOCI+ verison 密文除法
         */
        void sdiv(mpz_t eq, mpz_t ee, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp);


        /**
         * SOCI-TEE verison 密文乘法
         */
        void tee_smul(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * SOCI-TEE verison 密文比较
         */
        void tee_scmp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * SOCI-TEE verison 密文绝对值
         */
        void tee_sabs(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * SOCI-TEE verison 密文除法
         */
        void tee_sdiv(mpz_t eq, mpz_t ee, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * TRUST verison 密文乘法
         */
        void trust_fmul(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * TRUST verison 密文比较
         */
        void trust_fcmp(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * TRUST verison 密文相等判断
         */
        void trust_feql(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * TRUST verison 密文绝对值
         */
        void trust_fabs(mpz_t res, mpz_t em_1, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * TRUST verison 密文三目运算，论文版本
         */
        void trust_ftrn_version1(mpz_t res, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * TRUST verison 密文三目运算，本人自己开发更高效方案，非论文版本
         */
        void trust_ftrn_version2(mpz_t res, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * Zhao SOCI verison 密文除法
         */
        void online_smul(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * Zhao SOCI verison 密文比较
         */
        void online_scmp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * Zhao SOCI verison 密文符号位获取
         */
        void online_ssba(mpz_t s_x, mpz_t u_x, mpz_t ex, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * Zhao SOCI verison 密文除法
         */
        void online_sdiv(mpz_t eq, mpz_t ee, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp);


        /**
         * 改进 POCF verison，安全
         */
        void online_seqc(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * 改进 POCF verison，安全
         */
        void online_smms(mpz_t A, mpz_t I, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * 改进 POCF verison，安全
         */
        void online_sgcd(mpz_t res, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * 改进 POCF verison，安全
         */
        void online_sexp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * 改进 POCF verison，此版本不完善，不安全
         */
        void online_sinv(mpz_t res, mpz_t ex, mpz_t em, int ell, PaillierThdDec cp, PaillierThdDec csp, Paillier pai);

        void getRandom(mpz_t r, int len, gmp_randstate_t state);

        void getRandomwithUpper(mpz_t r, int len, mpz_t up, gmp_randstate_t state);

    };

    void protocol::smul(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
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
        /* send <C, C1> to csp */

        // Step 2, CSP computes
        mpz_t C2, L_mul_xr1_yr2, xr1, yr2, e_xr1_mul_yr2;
        mpz_inits(C2, L_mul_xr1_yr2, xr1, yr2, e_xr1_mul_yr2, NULL);
        csp.PDec(C2, C);
        csp.TDec(L_mul_xr1_yr2, C1, C2);
        mpz_divmod(xr1, yr2, L_mul_xr1_yr2, csp.pk.L);
        mpz_mul(e_xr1_mul_yr2, xr1, yr2);
        mpz_mod(e_xr1_mul_yr2, e_xr1_mul_yr2, csp.pk.N);
        csp.pai.encrypt(e_xr1_mul_yr2, e_xr1_mul_yr2);
        /* send <e_xr1_mul_yr2> to cp */

        // Step 3, CP computes
        mpz_t neg_er2x, neg_er1y, exy, neg_r1, neg_r2;
        mpz_inits(neg_er2x, neg_er1y, exy, neg_r1, neg_r2, NULL);
        mpz_neg(neg_r1, tuple_S0.r1);
        mpz_neg(neg_r2, tuple_S0.r2);
        cp.pai.scl_mul(neg_er2x, ex, neg_r2);
        cp.pai.scl_mul(neg_er1y, ey, neg_r1);
        cp.pai.add(exy, e_xr1_mul_yr2, neg_er2x);
        cp.pai.add(exy, exy, neg_er1y);
        cp.pai.add(exy, exy, neg_er1r2);
        mpz_set(res, exy);

        mpz_clears(X, Y, C, C1, C2, L_mul_xr1_yr2, xr1, yr2, e_xr1_mul_yr2, neg_er2x, neg_er1y, neg_r1, neg_r2, exy,
                   er1, er2, neg_er1r2, NULL);
    }

    void protocol::tee_smul(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;

        // Step 1, CP computes
        mpz_t r1, er1;
        mpz_inits(r1, er1, NULL);
        mpz_set(r1, tuple_S0.r1);
        cp.pai.add(er1, tuple_S0.er1, tuple_S0.e0);

        mpz_t X, Y, C, C1;
        mpz_inits(X, Y, C, C1, NULL);
        cp.pai.add(X, ex, er1);
        mpz_set(C, X);

        mpz_neg(r1, r1);
        cp.pai.scl_mul(Y, ey, r1);
        cp.PDec(C1, C);
        /* send <C, C1, [y], [-ry]> to csp */

        // Step 2, CSP computes
        mpz_t C2, x_add_r1;
        mpz_inits(C2, x_add_r1, NULL);
        csp.PDec(C2, C);
        csp.TDec(x_add_r1, C1, C2);
        csp.pai.scl_mul(res, ey, x_add_r1);
        csp.pai.add(res, res, Y);
        /* send [xy] to cp */

        mpz_clears(r1, er1, X, Y, C, C1, C2, x_add_r1, NULL);
    }

    void protocol::scmp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;
        Offline_val tuple_S1 = csp.pai.offlineVal;

        mpz_t er4, er3r4;
        mpz_inits(er4, er3r4, NULL);

        int pi = (int) random() % 2;

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

        // Step 2
        mpz_t D2, d, eu0, mid;
        mpz_inits(D2, d, eu0, mid, NULL);
        mpz_set(mid, csp.pk.N);
        mpz_div_ui(mid, mid, 2);
        csp.PDec(D2, D);
        csp.TDec(d, D1, D2);
        mpz_t e0_csp, e1_csp;
        mpz_inits(e0_csp, e1_csp, NULL);
        csp.pai.add(e0_csp, tuple_S1.e0, tuple_S1.e0);
        csp.pai.add(e1_csp, tuple_S1.e1, tuple_S1.e0);
        mpz_cmp(d, mid) > 0 ? mpz_set(eu0, e0_csp) : mpz_set(eu0, e1_csp);

        // Step 3
        mpz_set(res, eu0);
        mpz_t e1_cp;
        mpz_init(e1_cp);
        if (pi == 1)
        {
            cp.pai.add(e1_cp, tuple_S0.e1, tuple_S0.e0);
            cp.pai.sub(res, tuple_S0.e1, eu0);
        }

        mpz_clears(D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1, D2, d, eu0, mid, e0_csp, e1_csp, e1_cp, er4, er3r4,
                   NULL);
    }


    void protocol::tee_scmp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;

        mpz_t er4, er3r4;
        mpz_inits(er4, er3r4, NULL);

        int pi = (int) random() % 2;
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
        mpz_t epi, e_one_sub_pi;
        mpz_inits(epi, e_one_sub_pi, NULL);
        if (pi == 0)
        {
            cp.pai.add(epi, tuple_S0.e0, tuple_S0.e0);//refresh
            cp.pai.add(e_one_sub_pi, tuple_S0.e1, tuple_S0.e0);//refresh
        } else
        {
            mpz_set(epi, tuple_S0.e1);
            mpz_set(e_one_sub_pi, tuple_S0.e0);
            cp.pai.add(epi, tuple_S0.e1, tuple_S0.e0);//refresh
            cp.pai.add(e_one_sub_pi, tuple_S0.e0, tuple_S0.e0);//refresh
        }


        // Step 2
        mpz_t D2, d, u0, one_sub_u0, mid;
        mpz_inits(D2, d, u0, mid, one_sub_u0, NULL);
        mpz_set(mid, csp.pk.N);
        mpz_div_ui(mid, mid, 2);
        csp.PDec(D2, D);
        csp.TDec(d, D1, D2);
        if (mpz_cmp(d, mid) > 0)
        {
            mpz_set_si(u0, 0);
            mpz_set_si(one_sub_u0, 1);
        } else
        {
            mpz_set_si(u0, 1);
            mpz_set_si(one_sub_u0, 0);
        }

        mpz_t e1, e2;
        mpz_inits(e1, e2, NULL);
        csp.pai.scl_mul(e1, e_one_sub_pi, u0);
        csp.pai.scl_mul(e2, epi, one_sub_u0);
        csp.pai.add(res, e1, e2);

        mpz_clears(er4, er3r4, D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1, epi, e_one_sub_pi, D2, d, u0, one_sub_u0,
                   mid, e1, e2, NULL);
    }

    void protocol::tee_sabs(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;

        mpz_t er4, er3r4;
        mpz_inits(er4, er3r4, NULL);

        mpz_t e_x_sub_y;
        mpz_init(e_x_sub_y);
        cp.pai.sub(e_x_sub_y, ex, ey);

        int pi = (int) random() % 2;
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

        mpz_t one_sub_twopi;
        mpz_init(one_sub_twopi);
        if (pi == 0)
        {
            mpz_set_si(one_sub_twopi, 1);
        } else
        {
            mpz_set_si(one_sub_twopi, -1);
        }
        cp.pai.scl_mul(e_x_sub_y, e_x_sub_y, one_sub_twopi);

        // Step 2
        mpz_t D2, d, u0, mid, one_sub_twou0;
        mpz_inits(D2, d, u0, mid, one_sub_twou0, NULL);
        mpz_set(mid, csp.pk.N);
        mpz_div_ui(mid, mid, 2);
        csp.PDec(D2, D);
        csp.TDec(d, D1, D2);
        if (mpz_cmp(d, mid) > 0)
        {
            mpz_set_si(u0, 0);
            mpz_set_si(one_sub_twou0, 1);
        } else
        {
            mpz_set_si(u0, 1);
            mpz_set_si(one_sub_twou0, -1);
        }

        csp.pai.scl_mul(res, e_x_sub_y, one_sub_twou0);
        mpz_clears(er4, er3r4, D, D1, exr1, neg_eyr1, neg_r1, eyr1, neg_exr1, one_sub_twopi, D2, d, u0, mid,
                   one_sub_twou0, NULL);

    }

    void protocol::ssba(mpz_t s_x, mpz_t u_x, mpz_t ex, PaillierThdDec cp, PaillierThdDec csp)
    {
        Offline_val tuple_S0 = cp.pai.offlineVal;

        mpz_t e0, e1;
        mpz_inits(e0, e1, NULL);
        cp.pai.add(e0, tuple_S0.e0, tuple_S0.e0);//refresh
        cp.pai.add(e1, tuple_S0.e1, tuple_S0.e0);//refresh

        // Step-1
        scmp(s_x, ex, e0, cp, csp);

        // Step-2
        mpz_t n_sub_2;
        mpz_init(n_sub_2);
        mpz_sub_ui(n_sub_2, cp.pai.pk.N, 2);

        mpz_t sign;
        mpz_init(sign);
        cp.pai.scl_mul(sign, s_x, n_sub_2);   // [s_x]^(N-2)
        cp.pai.add(sign, e1, sign);    // [1]*[s_x^(N-2)]

        // Step-3
        smul(u_x, sign, ex, cp, csp);

        mpz_clears(e0, e1, n_sub_2, sign, NULL);

    }

    void protocol::tee_sdiv(mpz_t eq, mpz_t ee, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp)
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
            tee_scmp(u, ee, c, cp, csp);

            // Step-4
            cp.pai.scl_mul(u, u, neg_one);
            cp.pai.add(u, e1, u);
            cp.pai.scl_mul(ue, u, e);
            cp.pai.add(eq, eq, ue);

            // Step-5
            tee_smul(m, u, c, cp, csp);

            // Step-6
            cp.pai.scl_mul(m, m, neg_one);
            cp.pai.add(ee, ee, m);
        }
        mpz_clears(c, u, e, ue, m, two, neg_one, NULL);

    }

    void protocol::sdiv(mpz_t eq, mpz_t ee, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp)
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
            scmp(u, ee, c, cp, csp);

            // Step-4
            cp.pai.scl_mul(u, u, neg_one);
            cp.pai.add(u, e1, u);
            cp.pai.scl_mul(ue, u, e);
            cp.pai.add(eq, eq, ue);

            // Step-5
            smul(m, u, c, cp, csp);

            // Step-6
            cp.pai.scl_mul(m, m, neg_one);
            cp.pai.add(ee, ee, m);
        }
        mpz_clears(c, u, e, ue, m, two, neg_one, NULL);

    }


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


    void protocol::online_seqc(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {
        // jointly calculate
        mpz_t u1, u2, f1, f2, neg_one, eone;
        mpz_inits(u1, u2, f1, f2, neg_one, eone, NULL);

        online_scmp(u1, ex, ey, cp, csp);
        online_scmp(u2, ey, ex, cp, csp);


        mpz_set_si(neg_one, -1);
        mpz_set_si(eone, 1);
        cp.pai.encrypt(eone, eone);


        mpz_t newu1, newu2;
        mpz_inits(newu1, newu2, NULL);
        mpz_powm(newu1, u1, neg_one, cp.pai.pk.N_Square);
        mpz_powm(newu2, u2, neg_one, cp.pai.pk.N_Square);
        cp.pai.add(newu1, eone, newu1);
        cp.pai.add(newu2, eone, newu2);

        online_smul(f1, newu1, u2, cp, csp);
        online_smul(f2, u1, newu2, cp, csp);

        // CP calculates
        cp.pai.add(res, f1, f2);
    }

    void protocol::online_smms(mpz_t A, mpz_t I, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {
        // Step-1
        mpz_t eone, neg_one;
        mpz_inits(eone, neg_one, NULL);

        mpz_set_si(neg_one, -1);
        mpz_set_si(eone, 1);
        cp.pai.encrypt(eone, eone);

        mpz_t mid_N;
        mpz_init(mid_N);
        mpz_div_ui(mid_N, csp.pai.pk.N, 2);

        mpz_t r, r1, r2, r1r2, er, eyr, exr, er1, er2, er1r2, l1, l2, l3;
        mpz_inits(r, r1, r2, r1r2, er, eyr, exr, er1, er2, er1r2, l1, l2, l3, NULL);
        mpz_rrandomb(r, randstate, sigma);
        mpz_rrandomb(r1, randstate, sigma);
        getRandomwithUpper(r2, sigma, r1, randstate); // r2 = random([0, r1-1])
        mpz_sub(r2, mid_N, r2); // r2 = mid - r2
        mpz_add(r1r2, r1, r2);        // r1r2 = r1 + r2
        csp.pai.encrypt(er1, r1);
        csp.pai.encrypt(er2, r2);
        cp.pai.encrypt(er1r2, r1r2);    // er1r2 = [r1 + r2]

        if (mpz_cmp(r2, mid_N) > 0)
        {
            printf("error r2 > N/2\n");
            exit(-1);
        }
        if (mpz_cmp(r1r2, mid_N) <= 0)
        {
            printf("error r1+r2 <= N/2\n");
            exit(-1);
        }

        int pi = rand() % 2;
        if (pi == 1) // D = [r_1*(x-y+1)+r2]
        {
            // l1
            cp.pai.scl_mul(exr, ex, r1); // exr = [r1 * x]
            cp.pai.scl_mul(eyr, ey, neg_one); // eyr = [-y]
            cp.pai.scl_mul(eyr, eyr, r1); // eyr = [-r1 * y]

            cp.pai.add(l1, exr, eyr);

            cp.pai.add(l1, l1, er1r2);

            // l2
            mpz_powm(l2, ex, neg_one, cp.pai.N_Square);
            cp.pai.add(l2, er1, l2);
            cp.pai.add(l2, ey, l2);

            // l3
            mpz_powm(l3, ey, neg_one, cp.pai.N_Square);
            cp.pai.add(l3, er2, l3);
            cp.pai.add(l3, ex, l3);

        } else  // D = [r_1*(y-x)+r2]
        {
            // l1
            cp.pai.encrypt(er2, r2);
            cp.pai.scl_mul(eyr, ey, r1);
            cp.pai.scl_mul(exr, ex, neg_one);
            cp.pai.scl_mul(exr, exr, r1);

            cp.pai.add(l1, eyr, exr);
            cp.pai.add(l1, l1, er2);


            // l2
            mpz_powm(l2, ey, neg_one, cp.pai.N_Square);
            cp.pai.add(l2, er1, l2);
            cp.pai.add(l2, ex, l2);


            // l3
            mpz_powm(l3, ex, neg_one, cp.pai.N_Square);
            cp.pai.add(l3, er2, l3);
            cp.pai.add(l3, ey, l3);

        }
        mpz_t K1;
        mpz_init(K1);
        cp.PDec(K1, l1);

        // send<K1,l1,l2,l3>

        // Step-2
        mpz_t K, K2, u, eu, D1, D2;
        mpz_inits(K, K2, u, eu, D1, D2, NULL);
        csp.PDec(K2, l1);
        csp.TDec(K, K1, K2);

        if (mpz_cmp(K, mid_N) > 0)
        {
            // printf("K > N/2\n");
            mpz_set_si(u, 0);

            mpz_t zero, ezero;
            mpz_inits(zero, ezero, NULL);
            mpz_set_si(zero, 0);
            csp.pai.encrypt(ezero, zero);
            mpz_set(D1, ezero);
            mpz_set(D2, ezero);

            // Offline_val tuple_S0 = cp.pai.offlineVal;
            // cp.pai.add(D1, D1, tuple_S0.e0);//refresh
            // cp.pai.add(D2, D2, tuple_S0.e0);//refresh


        } else
        {

            // printf("K <= N/2\n");
            mpz_set_si(u, 1);
            // ciphertext refresh(CR)
            Offline_val tuple_S0 = cp.pai.offlineVal;
            cp.pai.add(D1, l2, tuple_S0.e0);//refresh
            cp.pai.add(D2, l3, tuple_S0.e0);//refresh


        }

        csp.pai.encrypt(eu, u);
        // send<u,D1, D2>

        // Step-3
        mpz_t neg_r1, neg_r2;
        mpz_inits(neg_r1, neg_r2, NULL);

        mpz_neg(neg_r1, r1);
        mpz_neg(neg_r2, r2);
        if (pi == 1)
        {
            // A
            mpz_powm(A, eu, neg_r1, cp.pai.N_Square);
            cp.pai.add(A, D1, A);
            cp.pai.add(A, ex, A);

            // I
            mpz_powm(I, eu, neg_r2, cp.pai.N_Square);
            cp.pai.add(I, D2, I);
            cp.pai.add(I, ey, I);

        } else
        {
            // A
            mpz_powm(A, eu, neg_r1, cp.pai.N_Square);
            cp.pai.add(A, D1, A);
            cp.pai.add(A, ey, A);

            // I
            mpz_powm(I, eu, neg_r2, cp.pai.N_Square);
            cp.pai.add(I, D2, I);
            cp.pai.add(I, ex, I);
        }


    }

    void protocol::online_sgcd(mpz_t res, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp)
    {
        mpz_t neg_one, eone, ezero;
        mpz_inits(neg_one, eone, ezero, NULL);
        mpz_set_si(eone, 1);
        mpz_set_si(ezero, 0);
        cp.pai.encrypt(eone, eone);
        cp.pai.encrypt(ezero, ezero);
        mpz_set_si(neg_one, -1);


        mpz_t A, I;
        mpz_inits(A, I, NULL);
        online_smms(A, I, ex, ey, cp, csp);


        vector <mpz_t> q_vec(ell + 1);
        vector <mpz_t> r_vec(ell + 1);
        vector <mpz_t> u_vec(ell + 1);
        vector <mpz_t> x_vec(ell + 1);
        vector <mpz_t> s_vec(ell + 1);
        for (int i = 0; i <= ell; i++)
        {
            mpz_inits(q_vec[i], r_vec[i], u_vec[i], x_vec[i], s_vec[i], NULL);
        }
        mpz_set(r_vec[0], I);

        for (int i = 1; i < ell; i++) //Euclidean Algorithm
        {
            // 确保I不等于0
            online_sdiv(q_vec[i], r_vec[i], A, I, ell, cp, csp);
            mpz_set(A, I); //
            mpz_set(I, r_vec[i]);// 余数

            // mpz_t temp1, temp2;
            // mpz_inits(temp1, temp2, NULL);
            // pai.decrypt(temp1, q_vec[i]);
            // pai.decrypt(temp2, r_vec[i]);
            // gmp_printf("q_vec[%d] = %Zd\n", i, temp1);
            // gmp_printf("r_vec[%d] = %Zd\n", i, temp2);
        }


        // printf("-------------------------------\n");
        for (int i = 1; i < ell; i++) // 非零余数表示为 1，将零余数表示为 0
        {

            online_seqc(u_vec[i], r_vec[i], ezero, cp, csp);
            // mpz_t temp1, temp2;
            // mpz_inits(temp1, temp2, NULL);
            // pai.decrypt(temp1, u_vec[i]);
            // gmp_printf("u_vec[%d] = %Zd\n", i, temp1);

        }
        // printf("-------------------------------\n");

        mpz_set(u_vec[0], eone);

        // 前缀异或
        for (int i = 1; i < ell; i++)
        {
            mpz_t f1, f2, f3, left1, right1, left2, right2;
            mpz_inits(f1, f2, f3, left1, right1, left2, right2, NULL);

            // f1
            mpz_powm(left1, u_vec[i - 1], neg_one, cp.pai.N_Square);
            cp.pai.add(left1, eone, left1);
            mpz_set(right1, u_vec[i]);
            online_smul(f1, left1, right1, cp, csp);

            // f2
            mpz_powm(right2, u_vec[i], neg_one, cp.pai.N_Square);
            cp.pai.add(right2, eone, right2);
            mpz_set(left2, u_vec[i - 1]);
            online_smul(f2, left2, right2, cp, csp);

            // f3
            cp.pai.add(f3, f1, f2);

            mpz_set(x_vec[i], f3);

            // mpz_t temp;
            // pai.decrypt(temp, x_vec[i]);
            // gmp_printf("x_vec[%d] = %Zd\n", i, temp);

        }
        // printf("-------------------------------\n");
        mpz_set(x_vec[0], ezero);
        mpz_set_si(res, 1);
        // 前缀减法
        for (int i = 1; i < ell; i++)
        {

            mpz_t sub, ans;
            mpz_inits(sub, ans, NULL);

            mpz_powm(sub, x_vec[i - 1], neg_one, cp.pai.N_Square);
            cp.pai.add(sub, x_vec[i], sub);

            // mpz_t temp;
            // pai.decrypt(temp, sub);
            // gmp_printf("s_vec[%d] = %Zd\n", i, temp);

            online_smul(ans, sub, r_vec[i - 1], cp, csp);
            cp.pai.add(res, res, ans);
        }

    }

    void protocol::online_sexp(mpz_t res, mpz_t x, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {

        // Step-1
        mpz_t r, er, y1, Y1, S, xr, neg_one;
        mpz_inits(r, er, y1, Y1, S, xr, neg_one, NULL);
        mpz_rrandomb(r, randstate, sigma);

        // printf("the bits of r is %ld\n", mpz_sizeinbase(r, 2));
        cp.pai.encrypt(er, r);
        cp.pai.add(y1, ey, er);

        mpz_set_si(neg_one, -1);
        mpz_powm(xr, x, r, cp.pai.pk.N);
        mpz_powm(S, xr, neg_one, cp.pai.pk.N);

        cp.PDec(Y1, y1);

        // send<Y1,y1>
        // Step-2
        mpz_t y2, Y2, h, H;
        mpz_inits(y2, Y2, h, H, NULL);
        csp.PDec(Y2, y1);
        csp.TDec(y1, Y1, Y2);
        mpz_powm(h, x, y1, csp.pai.pk.N);
        csp.pai.encrypt(H, h);

        // Step-3
        mpz_powm(res, H, S, cp.pai.pk.N_Square);

        mpz_clears(r, er, y1, Y1, S, xr, neg_one, NULL);
        mpz_clears(y2, Y2, h, H, NULL);

    }

    void protocol::online_sinv(mpz_t res, mpz_t ex, mpz_t em, int ell, PaillierThdDec cp, PaillierThdDec csp,
                               Paillier pai)
    {
        /*
        // Step-1
        mpz_t r, ex1, X1;
        mpz_inits(r, ex1, X1, NULL);
        mpz_rrandomb(r, randstate, sigma);
        mpz_powm(ex1, ex, r, cp.pai.pk.N_Square);
        cp.PDec(X1, ex1);
        // send<X1,ex1>

        // Step-2
        mpz_t neg_one, x1, X2, h, H;
        mpz_inits(neg_one, x1, X2, h, H,NULL);
        mpz_set_si(neg_one, -1);

        csp.PDec(X2, ex1);
        csp.TDec(x1, X1, X2);
        mpz_powm(h, x1, neg_one, csp.pai.pk.N);
        csp.pai.encrypt(H, h);

        // Step-3
        mpz_powm(res, H, r, cp.pai.pk.N_Square);


        mpz_clears(r, ex1, X1, NULL);
        mpz_clears(neg_one, x1, X2, h, H,NULL);*/

        mpz_t neg_one, eone, ezero;
        mpz_inits(neg_one, eone, ezero, NULL);
        mpz_set_si(eone, 1);
        mpz_set_si(ezero, 0);
        cp.pai.encrypt(eone, eone);
        cp.pai.encrypt(ezero, ezero);
        mpz_set_si(neg_one, -1);


        mpz_t A, I;
        mpz_inits(A, I, NULL);
        mpz_set(I, ex);
        mpz_set(A, em);

        vector <mpz_t> q_vec(ell + 1);
        vector <mpz_t> r_vec(ell + 1);
        vector <mpz_t> u_vec(ell + 1);
        vector <mpz_t> x_vec(ell + 1);
        vector <mpz_t> s_vec(ell + 1);
        vector <mpz_t> pres_vec(ell + 1);
        for (int i = 0; i <= ell; i++)
        {
            mpz_inits(q_vec[i], r_vec[i], u_vec[i], x_vec[i], s_vec[i], pres_vec[i], NULL);
        }
        mpz_set(r_vec[0], I);

        for (int i = 1; i < ell; i++) //Euclidean Algorithm
        {
            // 确保I不等于0
            online_sdiv(q_vec[i], r_vec[i], A, I, ell, cp, csp);
            mpz_set(A, I); //
            mpz_set(I, r_vec[i]);// 余数

            mpz_t temp1, temp2;
            mpz_inits(temp1, temp2, NULL);
            pai.decrypt(temp1, q_vec[i]);
            pai.decrypt(temp2, r_vec[i]);
//	       gmp_printf("q_vec[%d] = %Zd\n", i, temp1);
//	       gmp_printf("r_vec[%d] = %Zd\n", i, temp2);
        }


//        printf("-------------------------------\n");
        for (int i = 1; i < ell; i++) // 非零余数表示为 1，将零余数表示为 0
        {

            online_seqc(u_vec[i], r_vec[i], ezero, cp, csp);
            mpz_t temp1, temp2;
            mpz_inits(temp1, temp2, NULL);
            pai.decrypt(temp1, u_vec[i]);
            //gmp_printf("u_vec[%d] = %Zd\n", i, temp1);

        }
        //printf("-------------------------------\n");

        mpz_set(u_vec[0], eone);

        // 前缀异或
        for (int i = 1; i < ell; i++)
        {
            mpz_t f1, f2, f3, left1, right1, left2, right2;
            mpz_inits(f1, f2, f3, left1, right1, left2, right2, NULL);

            // f1
            mpz_powm(left1, u_vec[i - 1], neg_one, cp.pai.N_Square);
            cp.pai.add(left1, eone, left1);
            mpz_set(right1, u_vec[i]);
            online_smul(f1, left1, right1, cp, csp);

            // f2
            mpz_powm(right2, u_vec[i], neg_one, cp.pai.N_Square);
            cp.pai.add(right2, eone, right2);
            mpz_set(left2, u_vec[i - 1]);
            online_smul(f2, left2, right2, cp, csp);

            // f3
            cp.pai.add(f3, f1, f2);

            mpz_set(x_vec[i], f3);

            mpz_t temp;
            pai.decrypt(temp, x_vec[i]);
//            gmp_printf("x_vec[%d] = %Zd\n", i, temp);

        }
//         printf("-------------------------------\n");
        mpz_set(x_vec[0], ezero);
        // 前缀减操作
        for (int i = 1; i < ell; i++)
        {
            mpz_t sub;
            mpz_inits(sub, NULL);
            mpz_powm(sub, x_vec[i - 1], neg_one, cp.pai.N_Square);
            cp.pai.add(sub, x_vec[i], sub);
            mpz_set(pres_vec[i], sub);

            mpz_t temp;
            pai.decrypt(temp, pres_vec[i]);
//	           gmp_printf("pres_vec[%d] = %Zd\n", i, temp);
        }
//         printf("-------------------------------\n");


        for (int i = 1; i < ell; i++)
        {
            // 非操作
            mpz_t sub, ans;
            mpz_inits(sub, ans, NULL);
            mpz_powm(sub, x_vec[i], neg_one, cp.pai.N_Square);
            cp.pai.add(s_vec[i], eone, sub);

            // 密文加
            cp.pai.add(s_vec[i], s_vec[i], pres_vec[i]);

            mpz_t temp;
            pai.decrypt(temp, s_vec[i]);
//	           gmp_printf("s_vec[%d] = %Zd\n", i, temp);



        }
        mpz_set(s_vec[ell], ezero);
//         printf("-------------------------------\n");

        mpz_t t, x, y;
        mpz_inits(t, x, y, NULL);
        mpz_set(x, eone);
        mpz_set(y, ezero);

        for (int i = 3; i >= 1; i--)
        {

            mpz_t temp0, temp1, temp2, temp3;
            mpz_inits(temp0, temp1, temp2, temp3, NULL);
            online_smul(q_vec[i], q_vec[i], s_vec[i], cp, csp);


            // 这部分还差加个变与不变的，应该配合CSP

            pai.decrypt(temp0, q_vec[i]);

            online_smul(t, q_vec[i], y, cp, csp);
            mpz_t sub;
            mpz_init(sub);
            mpz_powm(sub, t, neg_one, cp.pai.N_Square);
            cp.pai.add(t, x, sub); // [t] = [x-t]
            online_smul(t, s_vec[i], t, cp, csp);
            // mpz_powm(t, t, eone, em);

            mpz_set(x, y);
            mpz_set(y, t);
            pai.decrypt(temp3, t);
            pai.decrypt(temp2, y);
            pai.decrypt(temp1, x);
            mpz_set(res, t);

        }


    }


    void protocol::online_smul(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {

        mpz_t neg_one, eone, ezero;
        mpz_inits(neg_one, eone, ezero, NULL);
        mpz_set_si(eone, 1);
        mpz_set_si(ezero, 0);
        cp.pai.encrypt(eone, eone);
        cp.pai.encrypt(ezero, ezero);
        mpz_set_si(neg_one, -1);


        // Step-1 CP
        mpz_t r1, r2, er1, er2, X, Y, X1, Y1, r1r2, er1r2;
        mpz_inits(r1, r2, er1, er2, X, Y, X1, Y1, r1r2, er1r2, NULL);

        mpz_rrandomb(r1, randstate, sigma);
        mpz_rrandomb(r2, randstate, sigma);
        cp.pai.encrypt(er1, r1);
        cp.pai.encrypt(er2, r2);
        cp.pai.add(X, ex, er1);
        cp.pai.add(Y, ey, er2);
        cp.PDec(X1, X);
        cp.PDec(Y1, Y);


        // Step-2 CSP
        mpz_t X2, Y2, maskedX, maskedY, maskedXY, exy;
        mpz_inits(X2, Y2, maskedX, maskedY, maskedXY, exy, NULL);
        csp.PDec(X2, X);
        csp.PDec(Y2, Y);
        csp.TDec(maskedX, X1, X2);
        csp.TDec(maskedY, Y1, Y2);

        mpz_mul(maskedXY, maskedX, maskedY);
        mpz_mod(maskedXY, maskedXY, csp.pai.pk.N);
        csp.pai.encrypt(exy, maskedXY);

        // Step-3 CP
        mpz_t exr2, eyr1;
        mpz_inits(exr2, eyr1, NULL);

        mpz_powm(exr2, ex, r2, cp.pai.pk.N_Square);
        mpz_powm(eyr1, ey, r1, cp.pai.pk.N_Square);
        mpz_mul(r1r2, r1, r2);
        cp.pai.encrypt(er1r2, r1r2);


        cp.pai.add(res, exr2, eyr1);
        cp.pai.add(res, res, er1r2);

        cp.pai.scl_mul(res, res, neg_one);
        cp.pai.add(res, exy, res);

        mpz_clears(exr2, eyr1, NULL);
        mpz_clears(r1, r2, er1, er2, X, Y, X1, Y1, r1r2, er1r2, NULL);

    }

    void protocol::online_scmp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp)
    {
        mpz_t mid_N;
        mpz_init(mid_N);
        mpz_div_ui(mid_N, csp.pai.pk.N, 2);

        mpz_t neg_one, eone, ezero;
        mpz_inits(neg_one, eone, ezero, NULL);
        mpz_set_si(eone, 1);
        mpz_set_si(ezero, 0);
        cp.pai.encrypt(eone, eone);
        cp.pai.encrypt(ezero, ezero);
        mpz_set_si(neg_one, -1);

        mpz_t r1, r2, r0, r1r2, er2, D, D1, exr, eyr;
        mpz_inits(r1, r2, r0, r1r2, er2, D, D1, exr, eyr, NULL);

        mpz_t mid;
        mpz_init(mid);
        mpz_set(mid, mid_N); // mid = N/2
        getRandom(r1, sigma, randstate); // r1 = random(sigma)
        getRandomwithUpper(r2, sigma, r1, randstate); // r2 = random([0, r1-1])
        mpz_sub(r2, mid, r2); // r2 = mid - r2

        mpz_add(r1r2, r1, r2);        // r1r2 = r1 + r2
        if (mpz_cmp(r2, mid_N) > 0)
        {
            printf("error r2 > N/2\n");
            exit(-1);
        }
        if (mpz_cmp(r1r2, mid_N) <= 0)
        {
            printf("error r1+r2 <= N/2\n");
            exit(-1);
        }
        int pi = rand() % 2;
        if (pi == 0)         // D = [r_1*(x-y+1)+r2]
        {

            //printf("pi = 0\n");
            mpz_add(r2, r1, r2);        // r2 = r1 + r2
            cp.pai.encrypt(er2, r2);    // er2 = [r1+r2]
            cp.pai.scl_mul(exr, ex, r1); // ex = [r1 * x]

            cp.pai.scl_mul(eyr, ey, neg_one); // ey = [-r1 * y]
            cp.pai.scl_mul(eyr, eyr, r1); // ey = [-r1 * y]
            cp.pai.add(D, exr, eyr);
            cp.pai.add(D, D, er2);


        } else                              // D = [r_1*(y-x)+r2]
        {
            //printf(" pi = 1 \n");
            cp.pai.encrypt(er2, r2);
            cp.pai.scl_mul(eyr, ey, r1);
            cp.pai.scl_mul(exr, ex, neg_one);
            cp.pai.scl_mul(exr, exr, r1);

            cp.pai.add(D, eyr, exr);
            cp.pai.add(D, D, er2);
        }
        cp.PDec(D1, D);

        //Step-2
        mpz_t d, D2;
        mpz_inits(d, D2, NULL);
        csp.PDec(D2, D);
        csp.TDec(d, D1, D2);

        mpz_cmp(d, mid_N) > 0 ? mpz_set(res, ezero) : mpz_set(res, eone);

        //Step-3
        if (pi == 0)
        {
            mpz_set(res, res);
        } else
        {
            cp.pai.scl_mul(res, res, neg_one);
            cp.pai.add(res, eone, res);
        }
        mpz_clears(r1, r2, r0, r1r2, er2, D, D1, exr, eyr, mid, d, D2, NULL);

    }

    void protocol::online_ssba(mpz_t s_x, mpz_t u_x, mpz_t ex, PaillierThdDec cp, PaillierThdDec csp)
    {

        mpz_t neg_one, eone, ezero;
        mpz_inits(neg_one, eone, ezero, NULL);
        mpz_set_si(eone, 1);
        mpz_set_si(ezero, 0);
        cp.pai.encrypt(eone, eone);
        cp.pai.encrypt(ezero, ezero);
        mpz_set_si(neg_one, -1);

        // Step-1
        online_scmp(s_x, ex, ezero, cp, csp);

        // Step-2
        mpz_t sign;
        mpz_init(sign);

        mpz_t n_sub_2;
        mpz_init(n_sub_2);
        mpz_sub_ui(n_sub_2, cp.pai.pk.N, 2);

        cp.pai.scl_mul(sign, s_x, n_sub_2);   // [s_x]^(N-2)
        cp.pai.add(sign, eone, sign);    // [1]*[s_x^(N-2)]

        // Step-3
        online_smul(u_x, sign, ex, cp, csp);
        mpz_clear(sign);


    }

    void
    protocol::online_sdiv(mpz_t eq, mpz_t ee, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp)
    {
        // Step-1

        mpz_t c, u, e, ue, m, two;
        mpz_inits(c, u, e, ue, m, two, NULL);
        mpz_set_si(two, 2);
        mpz_set(ee, ex);

        mpz_t neg_one, eone, ezero;
        mpz_inits(neg_one, eone, ezero, NULL);
        mpz_set_si(eone, 1);
        mpz_set_si(ezero, 0);
        cp.pai.encrypt(eone, eone);
        cp.pai.encrypt(ezero, ezero);
        mpz_set_si(neg_one, -1);
        mpz_set(eq, ezero);

        for (int i = ell; i >= 0; i--)
        {
            // Step-2
            mpz_pow_ui(e, two, i);      // e=2^i
            cp.pai.scl_mul(c, ey, e);   // [y]^{2^i}

            // Step-3
            online_scmp(u, ee, c, cp, csp);

            // Step-4
            cp.pai.scl_mul(u, u, neg_one);
            cp.pai.add(u, eone, u);
            cp.pai.scl_mul(ue, u, e);
            cp.pai.add(eq, eq, ue);

            // Step-5
            online_smul(m, u, c, cp, csp);

            // Step-6
            cp.pai.scl_mul(m, m, neg_one);
            cp.pai.add(ee, ee, m);
        }
        mpz_clears(c, u, e, ue, m, two, NULL);

    }

    void protocol::getRandom(mpz_t r, int len, gmp_randstate_t state)
    {

        do
        {
            mpz_urandomb(r, state, len);
        } while (mpz_cmp_ui(r, 0) == 0);
    }

    void protocol::getRandomwithUpper(mpz_t r, int len, mpz_t up, gmp_randstate_t state)
    {
        do
        {
            getRandom(r, len, state);
        } while (mpz_cmp_ui(r, 0) == 0 || mpz_cmp(r, up) >= 0);
    }

    protocol::protocol()
    {}

    protocol::~protocol()
    {}

    void check_sabs(int i, mpz_t result_cmp_x_y, mpz_t x, mpz_t y, Paillier pai)
    {
        mpz_t z, answer;
        mpz_inits(z, answer, NULL);
        if (mpz_cmp(x, y) >= 0)
        {
            mpz_sub(answer, x, y);
        } else mpz_sub(answer, y, x);

        pai.decrypt(z, result_cmp_x_y);
        if (mpz_cmp(z, answer) != 0)
        {
            gmp_printf("dec answer x>=y? = %Zd\n", z);
            gmp_printf("right answer %Zd\n", answer);
            printf("SABS error!\n");
            printf("i = %d is error \n", i);
            exit(-1);
        }
        mpz_clears(z, answer, NULL);
    }

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

