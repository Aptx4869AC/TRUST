/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#ifndef PAILLIER_H
#define PAILLIER_H
#pragma once

#include <iostream>
#include <gmp.h>
#include <omp.h>
#include <ctime>
#include <cstring>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <netinet/tcp.h>
#include <arpa/inet.h>
#include <fstream>
#include <vector>
#include <chrono>
#include <thread>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <fcntl.h>
#include <sys/uio.h>
#include <errno.h>


static gmp_randstate_t randstate;
const int sigma = 128;

using namespace std;

namespace PHESPACE {

    class PaillierKey {

    public:
        mpz_t g, n, nsquare, half_n, mid;

        PaillierKey() {
            mpz_inits(this->g, this->n, this->nsquare, this->half_n, this->mid, NULL);
        }

        PaillierKey(mpz_t p, mpz_t q) {
            mpz_t two;
            mpz_inits(this->g, this->n, this->nsquare, this->half_n, two, this->mid, NULL);
            mpz_mul(this->n, p, q);        // n = p * q
            mpz_add_ui(this->g, this->n, 1);    // g = n + 1
            mpz_mul(this->nsquare, this->n, this->n); // nsqaure = n * n
            mpz_set_si(two, 2);
            mpz_fdiv_q(this->half_n, this->n, two);      // half_n = n / 2
        }

        PaillierKey(mpz_t n) {
            mpz_t two;
            mpz_inits(this->g, this->n, this->nsquare, this->half_n, two, this->mid, NULL);
            mpz_set(this->n, n);
            mpz_add_ui(this->g, this->n, 1);    // g = n + 1;
            mpz_mul(this->nsquare, this->n, this->n); // nsqaure = n * n;
            mpz_set_si(two, 2);
            mpz_fdiv_q(this->half_n, this->n, two);      // half_n = n / 2
            mpz_set(this->mid, this->half_n); //mid = n/2
        }

        PaillierKey(mpz_t g, mpz_t n, mpz_t nsqaure) {
            mpz_t two;
            mpz_inits(this->g, this->n, this->nsquare, this->half_n, two, this->mid, NULL);

            mpz_set(this->g, g);
            mpz_set(this->n, n);
            mpz_set(this->nsquare, nsqaure);
            mpz_set_si(two, 2);
            mpz_fdiv_q(this->half_n, this->n, two);      // half_n = n / 2
        }


    };

    class PaillierPrivateKey : public PaillierKey {

    public:
        mpz_t lambda, lmdInv;

        PaillierPrivateKey() : PaillierKey() {
            mpz_inits(this->lambda, this->lmdInv, NULL);
        }


        PaillierPrivateKey(mpz_t p, mpz_t q, mpz_t lambda) : PaillierKey(p, q) {
            mpz_inits(this->lambda, this->lmdInv, NULL);

            mpz_set(this->lambda, lambda);
            mpz_invert(this->lmdInv, this->lambda, this->n);
        }

        PaillierPrivateKey(mpz_t n, mpz_t lambda) : PaillierKey(n) {
            mpz_inits(this->lambda, this->lmdInv, NULL);

            mpz_set(this->lambda, lambda);
            mpz_invert(this->lmdInv, this->lambda, this->n);
        }

    };

    class PaillierThdPrivateKey {

    public:
        mpz_t sk, n, nsqaure;

        PaillierThdPrivateKey() {
            mpz_inits(this->sk, this->n, this->nsqaure, NULL);
        }

        PaillierThdPrivateKey(mpz_t sk, mpz_t n, mpz_t nsqaure) {
            mpz_inits(this->sk, this->n, this->nsqaure, NULL);

            mpz_set(this->sk, sk);
            mpz_set(this->n, n);
            mpz_set(this->nsqaure, nsqaure);
        }
    };

    class Paillier {

    public:
        PaillierKey pubkey;
        PaillierPrivateKey prikey;

        Paillier() {}

        Paillier(PaillierKey pubkey) {
            this->pubkey = pubkey;
        }

        Paillier(PaillierPrivateKey prikey) {
            this->pubkey = PaillierKey(prikey.g, prikey.n, prikey.nsquare);
            this->prikey = prikey;
        }

        Paillier(PaillierKey pubkey, PaillierPrivateKey prikey) {
            this->pubkey = pubkey;
            this->prikey = prikey;
        }

        void keygen(mpz_t p, mpz_t q);

        void keygen(unsigned long bitLen);

        void encrypt(mpz_t c, mpz_t m);

        void encrypt(mpz_t c, mpz_t m, mpz_t r);

        void decrypt(mpz_t m, mpz_t c);

        void add(mpz_t res, mpz_t c1, mpz_t c2);

        void scl_mul(mpz_t resc, mpz_t c, mpz_t e);

        void scl_mul(mpz_t res, mpz_t c, int e);
    };

    class PaillierThdDec {

    public:
        PaillierThdPrivateKey psk;
        Paillier pai;

        PaillierThdDec(PaillierThdPrivateKey psk) : psk(psk) {};

        PaillierThdDec(PaillierThdPrivateKey psk, PaillierKey pubkey) : psk(psk), pai(pubkey) {};

        void pdec(mpz_t pc, mpz_t c);

        void fdec(mpz_t m, mpz_t c1, mpz_t c2);
    };

    void setrandom(gmp_randstate_t *randstate) {
        time_t seed = time(NULL);
        gmp_randinit_default(*randstate);
        gmp_randseed_ui(*randstate, (unsigned long int) seed);
    }

    void Paillier::keygen(unsigned long bitLen) {

        mpz_t p, q;
        mpz_inits(p, q, NULL);
        setrandom(&randstate);
        do {
            mpz_rrandomb(p, randstate, bitLen);
        } while (mpz_probab_prime_p(p, 80) == 0 || mpz_sizeinbase(p, 2) != bitLen);
        setrandom(&randstate);
        do {
            mpz_rrandomb(q, randstate, bitLen);
        } while (mpz_probab_prime_p(q, 80) == 0 || mpz_cmp(p, q) == 0 || mpz_sizeinbase(q, 2) != bitLen);

        //gmp_printf("p = %Zd\n", p);
        //gmp_printf("q = %Zd\n", q);
        // 获取 p 和 q 的位数
        int pBits = mpz_sizeinbase(p, 2);
        int qBits = mpz_sizeinbase(q, 2);

        keygen(p, q);
    }

    void Paillier::keygen(mpz_t p, mpz_t q) {

        mpz_t n, lambda;
        mpz_inits(n, lambda, NULL);
        mpz_mul(n, p, q);

        int nBits = mpz_sizeinbase(n, 2);
        //printf("N 的位数：%d\n", nBits);

        pubkey = PaillierKey(n);
        mpz_sub_ui(p, p, 1);
        mpz_sub_ui(q, q, 1);
        mpz_mul(lambda, p, q);

        prikey = PaillierPrivateKey(n, lambda);
    }

    PaillierThdPrivateKey *thdkeygen(PaillierPrivateKey prikey, int sigma) {

        mpz_t sk1, sk2;
        mpz_inits(sk1, sk2, NULL);
        setrandom(&randstate);
        mpz_rrandomb(sk1, randstate, sigma);    // sk1 is a ranodm number with sigma bits
        mpz_mul(sk2, prikey.lambda, prikey.lmdInv);
        mpz_sub(sk2, sk2, sk1);                // sk2 = lambda ?mu - sk1

        static PaillierThdPrivateKey res[2] = {PaillierThdPrivateKey(sk1, prikey.n, prikey.nsquare),
                                               PaillierThdPrivateKey(sk2, prikey.n, prikey.nsquare)
        };

        return res;
    }
    

    void Paillier::encrypt(mpz_t c, mpz_t m) {
        mpz_t r;
        mpz_init(r);
        mpz_urandomm(r, randstate, pubkey.n);
        encrypt(c, m, r);
    }

    void Paillier::encrypt(mpz_t c, mpz_t m, mpz_t r) {
//        double start_time, end_time;
//        start_time = omp_get_wtime();
        if (mpz_cmp_ui(m, 0) < 0) {
            mpz_add(m, m, pubkey.n);
        }
        mpz_mul(c, m, pubkey.n);             // m·n
        mpz_add_ui(c, c, 1);             // 1 + m·n
        mpz_powm(r, r, pubkey.n, pubkey.nsquare);  // r^n mod n^2
        mpz_mul(c, c, r);                // (1+m·N)·r^n
        mpz_mod(c, c, pubkey.nsquare);             // (1+m·N)·r^n mod n^2
//        end_time = omp_get_wtime();
//        printf("encrypted time is  ------  %f ms\n", (end_time - start_time) * 1000);
    }

    void Paillier::decrypt(mpz_t m, mpz_t c) {
        if (mpz_cmp(c, prikey.nsquare) >= 0) {
            throw ("ciphertext must be less than n^2");
            return;
        }

        // c=c^lambda mod n^2

        mpz_powm(m, c, prikey.lambda, prikey.nsquare);

        // (c - 1) / n * lambda^(-1) mod n
        mpz_sub_ui(m, m, 1);            // c=c-1
        mpz_fdiv_q(m, m, prikey.n);        // c=(c-1)/n
        mpz_mul(m, m, prikey.lmdInv);    // c=c*lambda^(-1)
        mpz_mod(m, m, prikey.n);        // m=c mod n
        if (mpz_cmp(m, prikey.mid) > 0)
            mpz_sub(m, m, prikey.n);
    }

    void Paillier::add(mpz_t res, mpz_t c1, mpz_t c2) {
        mpz_mul(res, c1, c2);
        mpz_mod(res, res, pubkey.nsquare);
    }

    void Paillier::scl_mul(mpz_t res, mpz_t c, int e) {
        mpz_t mp_e;
        mpz_init(mp_e);
        mpz_set_si(mp_e, e);
        scl_mul(res, c, mp_e);
    }

    void Paillier::scl_mul(mpz_t res, mpz_t c, mpz_t e) {
        mpz_powm(res, c, e, pubkey.nsquare);
    }


    void PaillierThdDec::pdec(mpz_t pc, mpz_t c) {
        // c^sk % n^2
        mpz_powm(pc, c, psk.sk, psk.nsqaure);
    }


    void PaillierThdDec::fdec(mpz_t m, mpz_t c1, mpz_t c2) {

        // (c1 * c2 % n^2 - 1)/n
        mpz_mul(m, c1, c2);
        mpz_mod(m, m, psk.nsqaure);
        mpz_sub_ui(m, m, 1);
        mpz_fdiv_q(m, m, psk.n);
    }
}

#endif

