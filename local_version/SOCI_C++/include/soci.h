/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#ifndef SOCI_H
#define SOCI_H
#pragma once

#include <gmp.h>
#include "paillier.h"

using namespace std;
using namespace PHESPACE;

const int sigma = 128;

namespace SOCISPACE {

    class soci {

    public:
        mpz_t zero, one, neg_one, neg_two;

        soci() {
            mpz_inits(this->zero, this->one, this->neg_one, this->neg_two, NULL);
            mpz_set_si(this->zero, 0);
            mpz_set_si(this->one, 1);
            mpz_set_si(this->neg_one, -1);
            mpz_set_si(this->neg_two, -2);
        }

        /**
         * SOCI乘法
         * @param res
         * @param ex
         * @param ey
         * @param cp
         * @param csp
         */
        void smul(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * SOCI比较
         * @param res
         * @param ex
         * @param ey
         * @param cp
         * @param csp
         */
        void scmp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * SOCI符号位获取
         * @param s_x
         * @param u_x
         * @param c
         * @param cp
         * @param csp
         */
        void ssba(mpz_t s_x, mpz_t u_x, mpz_t c, PaillierThdDec cp, PaillierThdDec csp);

        /**
         * SOCI除法
         * @param eq
         * @param er
         * @param ex
         * @param ey
         * @param ell
         * @param cp
         * @param csp
         */
        void sdiv(mpz_t eq, mpz_t er, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp);

        void getRandom(mpz_t r, int len, gmp_randstate_t state);

        void getRandomwithUpper(mpz_t r, int len, mpz_t up, gmp_randstate_t state);


    };


    void get_secRandNum(mpz_t r, int SIGMA) {
        setrandom(&randstate);
        mpz_rrandomb(r, randstate, SIGMA); // 生成 SIGMA 位数的随机数
    }

    void get_secRandNum(mpz_t r) {
        get_secRandNum(r, sigma);
    }

    /**
     * SOCI乘法
     * @param res
     * @param ex
     * @param ey
     * @param cp
     * @param csp
     */
    void soci::smul(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp) {
        // Step-1 CP
        mpz_t r1, r2, er1, er2, X, Y, X1, Y1, r1r2, er1r2;
        mpz_inits(r1, r2, er1, er2, X, Y, X1, Y1, r1r2, er1r2, NULL);

        mpz_rrandomb(r1, randstate, sigma); // 生成 SIGMA 位数的随机数
        mpz_rrandomb(r2, randstate, sigma); // 生成 SIGMA 位数的随机数
        cp.pai.encrypt(er1, r1);
        cp.pai.encrypt(er2, r2);
        cp.pai.add(X, ex, er1);
        cp.pai.add(Y, ey, er2);
        cp.pdec(X1, X);
        cp.pdec(Y1, Y);

        // Step-2 CSP
        mpz_t X2, Y2, maskedX, maskedY, maskedXY, exy;
        mpz_inits(X2, Y2, maskedX, maskedY, maskedXY, exy, NULL);
        csp.pdec(X2, X);
        csp.pdec(Y2, Y);
        csp.fdec(maskedX, X1, X2);
        csp.fdec(maskedY, Y1, Y2);

        mpz_mul(maskedXY, maskedX, maskedY);
        mpz_mod(maskedXY, maskedXY, csp.pai.pubkey.n);
        csp.pai.encrypt(exy, maskedXY);

        // Step-3 CP
        mpz_t exr2, eyr1;
        mpz_inits(exr2, eyr1, NULL);
        mpz_powm(exr2, ex, r2, cp.pai.pubkey.nsquare);
        mpz_powm(eyr1, ey, r1, cp.pai.pubkey.nsquare);
        mpz_mul(r1r2, r1, r2);
        cp.pai.encrypt(er1r2, r1r2);

        cp.pai.add(res, exr2, eyr1);
        cp.pai.add(res, res, er1r2);

        cp.pai.scl_mul(res, res, neg_one);
        cp.pai.add(res, exy, res);

        mpz_clears(r1, r2, er1, er2, X, Y, X1, Y1, r1r2, er1r2, exr2, eyr1, NULL);
    }

    /**
     * SOCI比较
     * @param res
     * @param ex
     * @param ey
     * @param cp
     * @param csp
     */
    void soci::scmp(mpz_t res, mpz_t ex, mpz_t ey, PaillierThdDec cp, PaillierThdDec csp) {
        //Step-1
        mpz_t r1, r2, r0, r1r2, er2, D, D1, exr, eyr;
        mpz_inits(r1, r2, r0, r1r2, er2, D, D1, exr, eyr, NULL);
        setrandom(&randstate);

        mpz_t mid;
        mpz_init(mid);
        mpz_set(mid, cp.pai.pubkey.half_n); // mid = n
        getRandom(r1, sigma, randstate); // r1 = random(sigma)
        getRandomwithUpper(r2, sigma, r1, randstate); // r2 = random([0, r1-1])
        mpz_sub(r2, mid, r2); // r2 = mid - r2

        mpz_add(r1r2, r1, r2);        // r1r2 = r1 + r2
        if (mpz_cmp(r2, cp.pai.pubkey.half_n) > 0) {
            printf("r2 <= N/2\n");
            exit(-1);
        }
        if (mpz_cmp(r1r2, cp.pai.pubkey.half_n) <= 0) {
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

        //Step-2
        mpz_t d, D2;
        mpz_inits(d, D2, NULL);
        csp.pdec(D2, D);
        csp.fdec(d, D1, D2);

        if (mpz_cmp(d, csp.pai.pubkey.half_n) > 0) {
            csp.pai.encrypt(res, zero);
        } else {
            csp.pai.encrypt(res, one);
        }

        //Step-3
        if (pi == 0) {
            mpz_set(res, res);
        } else {
            cp.pai.scl_mul(res, res, neg_one);
            mpz_t eone;
            mpz_init(eone);
            cp.pai.encrypt(eone, one);
            cp.pai.add(res, eone, res);
            mpz_clear(eone);
        }
        mpz_clears(r1, r2, r0, r1r2, er2, D, D1, exr, eyr, mid, d, D2, NULL);
    }

    /**
     * SOCI符号位获取
     * @param s_x
     * @param u_x
     * @param cipher
     * @param cp
     * @param csp
     */
    void soci::ssba(mpz_t s_x, mpz_t u_x, mpz_t cipher, PaillierThdDec cp, PaillierThdDec csp) {
        // Step-1
        mpz_t ezero;
        mpz_init(ezero);
        cp.pai.encrypt(ezero, zero);
        scmp(s_x, cipher, ezero, cp, csp);

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
        smul(u_x, sign, cipher, cp, csp);
        mpz_clear(sign);
    }

    /**
     * SOCI除法
     * @param eq
     * @param er
     * @param ex
     * @param ey
     * @param ell
     * @param cp
     * @param csp
     */
    void soci::sdiv(mpz_t eq, mpz_t er, mpz_t ex, mpz_t ey, int ell, PaillierThdDec cp, PaillierThdDec csp) {
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

        for (int i = ell; i >= 0; i--) {
            // Step-2
            mpz_pow_ui(e, two, i);      // e=2^i
            cp.pai.scl_mul(c, ey, e);   // [y]^{2^i}

            // Step-3
            scmp(u, er, c, cp, csp);

            // Step-4
            cp.pai.scl_mul(u, u, neg_one);
            cp.pai.add(u, eone, u);
            cp.pai.scl_mul(ue, u, e);
            cp.pai.add(eq, eq, ue);

            // Step-5
            smul(m, u, c, cp, csp);

            // Step-6
            cp.pai.scl_mul(m, m, neg_one);
            cp.pai.add(er, er, m);
        }
        mpz_clears(c, u, e, ue, m, two, NULL);
    }

    void soci::getRandom(mpz_t r, int len, gmp_randstate_t state) {

        do {
            mpz_rrandomb(r, state, len);
        } while (mpz_cmp_ui(r, 0) == 0);
    }

    void soci::getRandomwithUpper(mpz_t r, int len, mpz_t up, gmp_randstate_t state) {
        do {
            getRandom(r, len, state);
        } while (mpz_cmp_ui(r, 0) == 0 || mpz_cmp(r, up) >= 0);
    }

}
#endif

