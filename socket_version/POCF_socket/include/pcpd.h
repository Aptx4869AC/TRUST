/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#ifndef POCD_H
#define POCD_H
#pragma once

#include <gmp.h>
#include <time.h>
#include <stdlib.h>
#include <cmath>
#include <vector>
#include <algorithm>
#include <omp.h>
#include <iostream>
#include <ctime>
#include <cstring>
#include <unistd.h>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <fstream>
#include <chrono>
#include <thread>
#include <iomanip>
#include <netinet/tcp.h>
#include <fcntl.h>
#include <sys/uio.h>
#include <errno.h>

using namespace std;

static gmp_randstate_t randstate;
const static int sigma = 128;

namespace PHESPACE
{

    void setrandom(gmp_randstate_t *randstate)
    {
//        time_t seed = time(NULL);
        gmp_randinit_default(*randstate);
//        gmp_randseed_ui(*randstate, (unsigned long int) seed);
    }


    // 公钥
    class PublicKey
    {
    public:
        mpz_t N;
        mpz_t N_square;
        mpz_t N_half;
        int sigma;

        PublicKey()
        {
            mpz_inits(this->N, this->N_square, this->N_half, NULL);
        }

        PublicKey(mpz_t N, int sigma)
        {
            mpz_inits(this->N, this->N_square, this->N_half, NULL);
            mpz_set(this->N, N);
            mpz_mul(this->N_square, N, N);
            mpz_fdiv_q_ui(this->N_half, N, 2);
            this->sigma = sigma;
        }

        ~PublicKey()
        {}
    };

    // 私钥
    class PrivateKey
    {
    public:
        mpz_t lambda;
        mpz_t lambda_reverse;
        mpz_t N;
        mpz_t N_half;
        mpz_t N_square;

        PrivateKey()
        {
            mpz_inits(this->lambda, this->lambda_reverse, this->N, this->N_square, this->N_half, NULL);
        }

        PrivateKey(mpz_t lambda, mpz_t lambda_reverse, mpz_t N)
        {
            mpz_inits(this->lambda, this->lambda_reverse, this->N, this->N_square, this->N_half, NULL);
            mpz_set(this->lambda, lambda);
            mpz_set(this->lambda_reverse, lambda_reverse);
            mpz_set(this->N, N);
            mpz_mul(this->N_square, N, N);
            mpz_fdiv_q_ui(this->N_half, N, 2);
        }

        ~PrivateKey()
        {}
    };

    class Paillier
    {
    public:

        PrivateKey private_key;
        PublicKey public_key;
        PrivateKey private_key_1;
        PrivateKey private_key_2;

        Paillier()
        { setrandom(&randstate); }

        ~Paillier()
        {}

        void KeyGen(int bitLen, int sigma, int yita);
    };

    class Paillier_Third_Party
    {
    public:
        PublicKey public_key;
        PrivateKey private_key;
        int sigma;

        Paillier_Third_Party()
        {
            setrandom(&randstate);
        }

        Paillier_Third_Party(PublicKey public_key, PrivateKey private_key)
        {
            setrandom(&randstate);
            this->public_key = public_key;
            this->private_key = private_key;
        }

        ~Paillier_Third_Party()
        {}

    };


    void Paillier::KeyGen(int bitLen, int sigma, int yita)
    {
        mpz_t p, q, N, lambda, u, lambda_1, lambda_2;
        mpz_inits(p, q, N, lambda, u, lambda_1, lambda_2, NULL);

        setrandom(&randstate);
        do
        {
            mpz_rrandomb(p, randstate, bitLen);
        } while (mpz_probab_prime_p(p, 80) == 0 || mpz_sizeinbase(p, 2) != bitLen);
        setrandom(&randstate);
        do
        {
            mpz_rrandomb(q, randstate, bitLen);
        } while (mpz_probab_prime_p(q, 80) == 0 || mpz_cmp(p, q) == 0 || mpz_sizeinbase(q, 2) != bitLen);

        mpz_mul(N, p, q);
        mpz_sub_ui(p, p, 1);
        mpz_sub_ui(q, q, 1);
        mpz_mul(lambda, p, q);
        mpz_fdiv_q_ui(lambda, lambda, 2);
        mpz_invert(u, lambda, N);

        PrivateKey private_key(lambda, u, N);
        PublicKey public_key(N, sigma);
        this->private_key = private_key;
        this->public_key = public_key;

        mpz_t lambda_1_generated;
        mpz_init(lambda_1_generated);
        mpz_rrandomb(lambda_1_generated, randstate, sigma);

        mpz_t temp;
        mpz_init(temp);
        mpz_mul(lambda_2, lambda, u);
        mpz_mul(temp, lambda, N);
        mpz_mul_ui(temp, temp, yita);
        mpz_add(lambda_2, lambda_2, temp);
        mpz_sub(lambda_2, lambda_2, lambda_1_generated);

        PrivateKey private_key_1(lambda_1_generated, u, N);
        PrivateKey private_key_2(lambda_2, u, N);

        this->private_key_1 = private_key_1;
        this->private_key_2 = private_key_2;
    }

    void L_function(mpz_t result, const mpz_t N)
    {
        mpz_sub_ui(result, result, 1);
        mpz_fdiv_q(result, result, N);
    }

    void add(mpz_t res, mpz_t c1, mpz_t c2, mpz_t N_square)
    {
        mpz_mul(res, c1, c2);
        mpz_mod(res, res, N_square);
    }

    void scl_mul(mpz_t res, mpz_t c, mpz_t e, mpz_t N_square)
    {
        mpz_powm(res, c, e, N_square);
    }

    void scl_mul(mpz_t res, mpz_t c, long long e, mpz_t N_square)
    {
        mpz_t mp_e;
        mpz_init(mp_e);
        mpz_set_si(mp_e, e);
        scl_mul(res, c, mp_e, N_square);
    }


    void Enc(mpz_t result, const PublicKey &public_key, mpz_t m)
    {
        if (mpz_cmp_si(m, 0) < 0)
        {
            mpz_add(m, m, public_key.N);
        }

        // 选择r
        mpz_t r;
        mpz_init(r);
        mpz_t flag;
        mpz_init(flag);
        while (true)
        {
            mpz_urandomm(r, randstate, public_key.N);
            mpz_gcd(flag, r, public_key.N);
            if (mpz_cmp_si(flag, 1) == 0)
            {
                break;
            }
        }
        mpz_urandomm(r, randstate, public_key.N);

        mpz_t temp1, temp2;
        mpz_inits(temp1, temp2, NULL);

//        mpz_t g;
//        mpz_init(g);
//        mpz_add_ui(g, public_key.N, 1);
//        mpz_powm(temp1, g, m, public_key.N_square); // g^m mod n^2


        // 计算 temp1 = (1 + m * N) % N_square
        mpz_mul(temp1, m, public_key.N);
        mpz_add_ui(temp1, temp1, 1);
        mpz_mod(temp1, temp1, public_key.N_square);


        // 计算 temp2 = r^N % N_square
        mpz_powm(temp2, r, public_key.N, public_key.N_square);

        // 计算 result = temp1 * temp2 % N_square
        mpz_mul(result, temp1, temp2);
        mpz_mod(result, result, public_key.N_square);

        mpz_clears(r, temp1, temp2, NULL);
    }

    void Dec(mpz_t result, const PrivateKey &private_key, const mpz_t cipher_text)
    {

        mpz_powm(result, cipher_text, private_key.lambda, private_key.N_square);
        L_function(result, private_key.N);
        mpz_mul(result, result, private_key.lambda_reverse);
        mpz_mod(result, result, private_key.N);

        if (mpz_cmp(result, private_key.N_half) > 0)
        {
            mpz_sub(result, result, private_key.N);
        }
    }


    void PDec(mpz_t result, const PrivateKey &private_key, const mpz_t ciphertext)
    {
        mpz_powm(result, ciphertext, private_key.lambda, private_key.N_square);
    }

    void TDec(mpz_t result, const mpz_t ciphertext1, const mpz_t ciphertext2,
              const mpz_t N_square, const mpz_t N, const mpz_t N_half)
    {

        // 计算 result = (ciphertext1 * ciphertext2) % (N * N)
        mpz_mul(result, ciphertext1, ciphertext2);
        mpz_mod(result, result, N_square);

        // 计算 result = L_function(result, N)，即m
        L_function(result, N);
        if (mpz_cmp(result, N_half) > 0)
        {
            mpz_sub(result, result, N);
        }
    }

    void CR(mpz_t ciphertext, const mpz_t N, const mpz_t N_square)
    {
        mpz_t r, flag;
        mpz_inits(r, flag, NULL);

        // 选择 r
        while (true)
        {
            mpz_urandomm(r, randstate, N);
            mpz_gcd(flag, r, N);
            if (mpz_cmp_ui(flag, 1) == 0)
            {
                break;
            }
        }

        // 计算 new_ciphertext = (old_ciphertext * gmpy2.powmod(r, N, N_square)) % N_square
        mpz_powm(r, r, N, N_square);
        mpz_mul(ciphertext, ciphertext, r);
        mpz_mod(ciphertext, ciphertext, N_square);

        mpz_clears(r, flag, NULL);
    }

}
#endif

