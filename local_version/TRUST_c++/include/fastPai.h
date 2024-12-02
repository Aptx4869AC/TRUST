/******************************************************************************
 * Author: APTX4869
 * Created: 2024-12-01
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#ifndef FASTPAI_H
#define FASTPAI_H

#pragma once

#include <gmp.h>
#include <omp.h>
#include <vector>
#include <random>
#include <ctime>
#include <iostream>
#include <cstring>
#include <random>
#include <fcntl.h>
#include <unistd.h>
#include <sys/uio.h>
#include "offline.h"
using namespace OFFLINESPACE;
using namespace std;

static gmp_randstate_t randstate;
static int sigma;
static int k; //BIT SECURITY , 64、80、104、112



ssize_t save_to_file(const char *filename, const vector <mpz_ptr> &src) {
    // Open the file for writing
    int fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd == -1) {
        perror("open");
        exit(EXIT_FAILURE);
    }
    vector<struct iovec> iov(src.size());
    for (size_t i = 0; i < src.size(); ++i) {
        size_t n = mpz_sizeinbase(src[i], 10);
        unsigned char *msg = new unsigned char[n];
        mpz_export(msg, &n, 1, 1, 0, 0, src[i]);

        iov[i].iov_base = msg;
        iov[i].iov_len = n;

    }
    ssize_t bytesWritten = writev(fd, iov.data(), src.size());
    if (bytesWritten == -1) {
        perror("writev");
        exit(EXIT_FAILURE);
    }
//        printf("%ld bytes written to file '%s'.\n", bytesWritten, filename);

    for (size_t i = 0; i < src.size(); ++i) {
        delete[] static_cast<char *>(iov[i].iov_base);
    }

    // Close the file descriptor
    close(fd);

    return bytesWritten;
}

namespace PHESPACE
{

    class PrivateKey
    {
    public:
        mpz_t alpha, N;

        PrivateKey()
        {
            mpz_inits(this->alpha, this->N, NULL);
        }

        PrivateKey(mpz_t alpha, mpz_t N)
        {
            mpz_inits(this->alpha, this->N, NULL);
            mpz_set(this->alpha, alpha);
            mpz_set(this->N, N);
        }

        ~PrivateKey()
        {}
    };

    class PublicKey
    {
    public:
        mpz_t N, h, h_N, L, N_Square;
        int L_k;

        PublicKey()
        {
            mpz_inits(this->N, this->h, this->h_N, this->L, this->N_Square, NULL);
        }

        PublicKey(mpz_t N, mpz_t h, mpz_t h_N, mpz_t L, int L_k)
        {
            mpz_inits(this->N, this->h, this->h_N, this->L, this->N_Square, NULL);
            mpz_set(this->N, N);
            mpz_set(this->h, h);
            mpz_set(this->h_N, h_N);
            mpz_set(this->L, L);
            mpz_mul(this->N_Square, N, N);
            this->L_k = L_k;
        }

        ~PublicKey()
        {}
    };

    class PaillierThirdParty
    {
    public:
        mpz_t N, partial_key;

        PaillierThirdParty(mpz_t N, mpz_t partial_key)
        {
            mpz_inits(this->N, this->partial_key, NULL);
            mpz_set(this->N, N);
            mpz_set(this->partial_key, partial_key);
        }

        ~PaillierThirdParty()
        {}
    };

    void setrandom(gmp_randstate_t *randstate)
    {
        time_t seed = time(NULL);
        gmp_randinit_default(*randstate);
        gmp_randseed_ui(*randstate, (unsigned long int) seed);
    }

    class Paillier
    {
    public:
        PublicKey pk;
        PrivateKey sk;
        mpz_t N_Square;
        mpz_t **table;
        Offline_val offlineVal;

        Paillier()
        {
            mpz_init(this->N_Square);
        }

        Paillier(PublicKey pk, int sigma)
        {
            mpz_init(this->N_Square);
            mpz_mul(this->N_Square, pk.N, pk.N);
            this->pk = pk;
            double start_time = omp_get_wtime();
            this->table = Pre_Compute::constructTable(pk.h_N, this->N_Square, pk.L_k);
            Offline_Phase(sigma);
            double end_time = omp_get_wtime();
            printf("Offline_Phase time is  ------  %f ms\n", (end_time - start_time) * 1000);
        }

        Paillier(PublicKey pk, PrivateKey sk, int sigma)
        {
            mpz_init(this->N_Square);
            mpz_mul(this->N_Square, pk.N, pk.N);
            this->pk = pk;
            this->sk = sk;
            double start_time = omp_get_wtime();
            this->table = Pre_Compute::constructTable(pk.h_N, this->N_Square, pk.L_k);
            Offline_Phase(sigma);
            double end_time = omp_get_wtime();
            printf("Offline_Phase time is  ------  %f ms\n", (end_time - start_time) * 1000);
        }

        ~Paillier()
        {}

        void Offline_Phase(int sigma);

        void encrypt(mpz_t dst, mpz_t m);

        void decrypt(mpz_t dst, mpz_t c);

        void add(mpz_t dst, mpz_t ex, mpz_t ey);

        void sub(mpz_t dst, mpz_t c1, mpz_t c2);

        void scl_mul(mpz_t dst, mpz_t c, int cons);

        void scl_mul(mpz_t dst, mpz_t c, mpz_t cons);
    };

    void Paillier::add(mpz_t dst, mpz_t ex, mpz_t ey)
    {
        mpz_mul(dst, ex, ey);
        mpz_mod(dst, dst, N_Square);
    }

    void Paillier::sub(mpz_t dst, mpz_t c1, mpz_t c2)
    {
        mpz_t neg_one;
        mpz_inits(neg_one, NULL);
        mpz_set_si(neg_one, -1);
        scl_mul(dst, c2, neg_one);
        add(dst, dst, c1);
        mpz_clears(neg_one, NULL);
    }

    void Paillier::scl_mul(mpz_t dst, mpz_t c, int cons)
    {
        mpz_t cons_mpz;
        mpz_inits(cons_mpz, NULL);
        mpz_set_si(cons_mpz, cons);
        scl_mul(dst, c, cons_mpz);
        mpz_clears(cons_mpz, NULL);
    }

    void Paillier::scl_mul(mpz_t dst, mpz_t c, mpz_t cons)
    {
        mpz_powm(dst, c, cons, N_Square);
    }

    static void L_function(mpz_t dst, mpz_t x, mpz_t N)
    {
        mpz_sub_ui(dst, x, 1);
        mpz_div(dst, dst, N);
    }


    void Paillier::Offline_Phase(int sigma)
    {
        mpz_inits(offlineVal.r1, offlineVal.r2, offlineVal.er1, offlineVal.er2, offlineVal.neg_er1r2, offlineVal.r3,
                  offlineVal.r4, offlineVal.er4, offlineVal.er3r4, offlineVal.e0, offlineVal.e1, NULL);
        mpz_urandomb(offlineVal.r1, randstate, sigma);
        mpz_urandomb(offlineVal.r2, randstate, sigma);
        mpz_urandomb(offlineVal.r3, randstate, sigma);
        encrypt(offlineVal.er1, offlineVal.r1);
        encrypt(offlineVal.er2, offlineVal.r2);
        mpz_mul(offlineVal.neg_er1r2, offlineVal.r1, offlineVal.r2);
        mpz_neg(offlineVal.neg_er1r2, offlineVal.neg_er1r2);
        encrypt(offlineVal.neg_er1r2, offlineVal.neg_er1r2);
        mpz_t half_N;
        mpz_init(half_N);
        mpz_div_ui(half_N, pk.N, 2);

        mpz_t tmp;
        mpz_init(tmp);
        mpz_urandomb(tmp, randstate, sigma);
        while (true)
        {
            if (mpz_cmp(tmp, offlineVal.r3) < 0)
            {
                break;
            }
            mpz_urandomb(tmp, randstate, sigma);
        }
        mpz_sub(offlineVal.r4, half_N, tmp);
        encrypt(offlineVal.er4, offlineVal.r4);
        mpz_add(offlineVal.er3r4, offlineVal.r3, offlineVal.r4);
        encrypt(offlineVal.er3r4, offlineVal.er3r4);
        mpz_t zero, one;
        mpz_inits(zero, one, NULL);
        mpz_set_si(zero, 0);
        mpz_set_si(one, 1);
        encrypt(offlineVal.e0, zero);
        encrypt(offlineVal.e1, one);

        mpz_inits(offlineVal.trust_r, offlineVal.trust_r1, offlineVal.trust_r2,
                  offlineVal.trust_r1_prime, offlineVal.trust_r2_prime, NULL);

        mpz_urandomb(offlineVal.trust_r, randstate, sigma);
        mpz_urandomb(offlineVal.trust_r1, randstate, sigma);
        mpz_urandomb(offlineVal.trust_r1_prime, randstate, sigma);


        mpz_urandomb(tmp, randstate, sigma);
        while (true)
        {
            if (mpz_cmp(tmp, offlineVal.trust_r1) < 0)
            {
                break;
            }
            mpz_urandomb(tmp, randstate, sigma);
        }
        mpz_sub(offlineVal.trust_r2, half_N, tmp);

        mpz_urandomb(tmp, randstate, sigma);
        while (true)
        {
            if (mpz_cmp(tmp, offlineVal.trust_r1_prime) < 0)
            {
                break;
            }
            mpz_urandomb(tmp, randstate, sigma);
        }
        mpz_sub(offlineVal.trust_r2_prime, half_N, tmp);

        mpz_inits(offlineVal.e_trust_r, offlineVal.e_trust_r1, offlineVal.e_trust_r2,
                  offlineVal.e_trust_r1_prime, offlineVal.e_trust_r2_prime, NULL);
        encrypt(offlineVal.e_trust_r, offlineVal.trust_r);
        encrypt(offlineVal.e_trust_r1, offlineVal.trust_r1);
        encrypt(offlineVal.e_trust_r2, offlineVal.trust_r2);
        encrypt(offlineVal.e_trust_r1_prime, offlineVal.trust_r1_prime);
        encrypt(offlineVal.e_trust_r2_prime, offlineVal.trust_r2_prime);


        mpz_inits(offlineVal.e_trust_r1_add_r2, offlineVal.e_trust_r1_prime_add_r2_prime, NULL);
        mpz_add(offlineVal.e_trust_r1_add_r2, offlineVal.trust_r1, offlineVal.trust_r2);
        encrypt(offlineVal.e_trust_r1_add_r2, offlineVal.e_trust_r1_add_r2);
        mpz_add(offlineVal.e_trust_r1_prime_add_r2_prime, offlineVal.trust_r1_prime, offlineVal.trust_r2_prime);
        encrypt(offlineVal.e_trust_r1_prime_add_r2_prime, offlineVal.e_trust_r1_prime_add_r2_prime);

        mpz_inits(offlineVal.e_trust_2r1_prime_add_r2_prime, offlineVal.e_trust_r2_prime_sub_r1_prime, NULL);
        mpz_add(offlineVal.e_trust_2r1_prime_add_r2_prime, offlineVal.trust_r1_prime, offlineVal.trust_r1_prime);
        mpz_add(offlineVal.e_trust_2r1_prime_add_r2_prime, offlineVal.e_trust_2r1_prime_add_r2_prime, offlineVal.trust_r2_prime);
        encrypt(offlineVal.e_trust_2r1_prime_add_r2_prime, offlineVal.e_trust_2r1_prime_add_r2_prime);
        mpz_sub(offlineVal.e_trust_r2_prime_sub_r1_prime, offlineVal.trust_r2_prime, offlineVal.trust_r1_prime);
        encrypt(offlineVal.e_trust_r2_prime_sub_r1_prime, offlineVal.e_trust_r2_prime_sub_r1_prime);

        mpz_clears(half_N, tmp, zero, one, NULL);

    }

    void Paillier::encrypt(mpz_t dst, mpz_t m)
    {
        if (mpz_cmp_ui(m, 0) < 0)
        {
            mpz_add(m, m, pk.N);

        }//如果m是负数
        mpz_t r;
        mpz_init(r);
        mpz_urandomb(r, randstate, pk.L_k);
        if (table != nullptr)
        {
            mpz_t mask;
            mpz_init(mask);
            Pre_Compute::compute(mask, r, table, N_Square);//(h^N mod N^2)^r
            mpz_mod(mask, mask, N_Square);//(h^N mod N^2)^r mod N^2
            mpz_mul(dst, m, pk.N);//m*N
            mpz_add_ui(dst, dst, 1);//1+m*N
            mpz_mod(dst, dst, N_Square);//1+m*N mod N^2
            mpz_mul(dst, dst, mask);
            mpz_mod(dst, dst, N_Square);
            mpz_clears(mask, NULL);
        } else
        {
            mpz_t tmp;
            mpz_init(tmp);
            mpz_powm(tmp, pk.h_N, r, N_Square);//正常操作

            mpz_mul(dst, m, pk.N);
            mpz_add_ui(dst, dst, 1);
            mpz_mod(dst, dst, N_Square);
            mpz_mul(dst, dst, tmp);
            mpz_mod(dst, dst, N_Square);
            mpz_clears(tmp, NULL);
        }
    }

    void Paillier::decrypt(mpz_t dst, mpz_t c)
    {
        mpz_t parameter, inverse, m, alpha2;
        mpz_inits(parameter, inverse, m, alpha2, NULL);
        mpz_mul_ui(alpha2, sk.alpha, 2);

        mpz_powm(parameter, c, alpha2, N_Square);
        mpz_invert(inverse, alpha2, sk.N);

        L_function(parameter, parameter, sk.N);
        mpz_mul(m, inverse, parameter);
        mpz_mod(m, m, sk.N);

        mpz_t mid;
        mpz_init(mid);
        mpz_div_ui(mid, sk.N, 2);
        if (mpz_cmp(m, mid) > 0)
        {
            mpz_sub(m, m, sk.N);
        }
        mpz_set(dst, m);
        mpz_clears(parameter, inverse, m, alpha2, NULL);
    }

    class PaillierThdDec
    {
    public:
        mpz_t N, partial_key, N_Square;
        PublicKey pk;
        Paillier pai;

        PaillierThdDec()
        {
            mpz_inits(this->N, this->partial_key, this->N_Square, NULL);
        }

        PaillierThdDec(mpz_t N, mpz_t partial_key, PublicKey pk, int sigma)
        {
            mpz_inits(this->N, this->partial_key, this->N_Square, NULL);
            mpz_set(this->N, N);
            mpz_set(this->partial_key, partial_key);
            mpz_mul(this->N_Square, N, N);
            this->pk = pk;
            this->pai = Paillier(pk, sigma);
        }

        ~PaillierThdDec()
        {}

        void PDec(mpz_t dst, mpz_t c);

        void TDec(mpz_t dst, mpz_t c1, mpz_t c2);
    };

    void PaillierThdDec::PDec(mpz_t dst, mpz_t c)
    {
        mpz_powm(dst, c, partial_key, N_Square);
    }

    void PaillierThdDec::TDec(mpz_t dst, mpz_t c1, mpz_t c2)
    {
        mpz_mul(dst, c1, c2);
        mpz_mod(dst, dst, N_Square);
        L_function(dst, dst, N);
    }

    class KeyGen
    {
    public:
        PublicKey pk;
        PrivateKey sk;
        PaillierThirdParty *pai_third_parties;

        KeyGen(int k, int sigma)
        {
            pai_third_parties = KGen(k, sigma);
        }

        ~KeyGen()
        {}

        PaillierThirdParty *KGen(int k, int sigma)
        {
            int yita = 0;
            mpz_t N, P, Q, p, q;
            mpz_inits(N, P, Q, p, q, NULL);

            mpz_t *params = Ngen(k);
            mpz_set(N, params[0]);
            mpz_set(P, params[1]);
            mpz_set(Q, params[2]);
            mpz_set(p, params[3]);
            mpz_set(q, params[4]);

            mpz_t alpha;
            mpz_init(alpha);
            mpz_mul(alpha, p, q);//私钥

            mpz_t beta;
            mpz_init(beta);
            mpz_sub_ui(P, P, 1); // P-1
            mpz_sub_ui(Q, Q, 1); // Q-1
            mpz_mul(beta, P, Q); // PQ
            mpz_mul_si(p, p, 4);
            mpz_mul(p, p, q); // 4PQ
            mpz_div(beta, beta, p);

            mpz_t y;
            mpz_init(y);
            mpz_t y_gcd_N;
            mpz_init(y_gcd_N);
            while (true)
            {
                mpz_urandomm(y, randstate, N);
                mpz_gcd(y_gcd_N, y, N);
                if (mpz_cmp_ui(y_gcd_N, 1) == 0)
                {
                    break;
                }
            }

            mpz_t h;
            mpz_init(h);
            mpz_mul_ui(h, beta, 2);
            mpz_neg(y, y);
            mpz_powm(h, y, h, N);

            mpz_t h_N;
            mpz_init(h_N);
            mpz_t N_square;
            mpz_init(N_square);
            mpz_mul(N_square, N, N);
            mpz_powm(h_N, h, N, N_square);

            mpz_t L;
            mpz_init(L);
            mpz_urandomb(L, randstate, sigma + 2);
            pk = PublicKey(N, h, h_N, L, getL_k(k));
            sk = PrivateKey(alpha, N);

            mpz_t alpha1;
            mpz_init(alpha1);
            mpz_urandomb(alpha1, randstate, sigma);

            mpz_t alpha2, tmp, two_alpha;
            mpz_inits(alpha2, tmp, two_alpha, NULL);
            mpz_mul_ui(two_alpha, alpha, 2);
            mpz_invert(tmp, two_alpha, N);

            mpz_mul(alpha2, two_alpha, tmp);
            mpz_sub(alpha2, alpha2, alpha1);
            mpz_mul_ui(tmp, N, yita);
            mpz_mul(tmp, tmp, two_alpha);
            mpz_add(alpha2, alpha2, tmp);

//            double total_bytes = 0.0;
//            vector<mpz_ptr> src_vec;
//            src_vec.push_back(N);
//            src_vec.push_back(h);
//            total_bytes += save_to_file("CP_output.txt", src_vec);
//            total_bytes /= 1024.0;
//            printf("bytes Costs is  ---  %f KB\n", total_bytes);
//            exit(-1);
            // two_alpha bytes Costs is  ---  0.055664 KB

            PaillierThirdParty *res = new PaillierThirdParty[2]{
                    PaillierThirdParty(N, alpha1), // 使用已有的构造函数来初始化对象
                    PaillierThirdParty(N, alpha2)  // 使用相同的参数再次初始化对象
            };
            mpz_clears(N, P, Q, p, q, alpha, beta, y, h, h_N, L, alpha1, alpha2, tmp, two_alpha, NULL);
            return res;
        }

        static int getN_k(int k)
        {
            switch (k)
            {
                case 64:
                    return 512;
                case 80:
                    return 1024;
                case 104:
                    return 1536;
                case 112:
                    return 2048;
                case 128:
                    return 3072;
                case 192:
                    return 7680;
                default:
                    return 2048;
            }
        }

        static int getL_k(int k)
        {
            return 4 * k;
        }

        bool areCoprime(const mpz_t p, const mpz_t q)
        {
            mpz_t gcd;
            mpz_init(gcd);

            // 计算 p 和 q 的最大公约数
            mpz_gcd(gcd, p, q);

            // 如果最大公约数是1，则 p 和 q 是互素的
            bool coprime = (mpz_cmp_ui(gcd, 1) == 0);

            // 释放内存
            mpz_clear(gcd);

            return coprime;
        }

        bool isCoPrime(mpz_t p, mpz_t q, mpz_t p_another, mpz_t q_another)
        {
            mpz_t p_another_gcd_q_another, p_another_mod_p, p_another_mod_q, q_another_mod_p, q_another_mod_q;
            mpz_inits(p_another_gcd_q_another, p_another_mod_p, p_another_mod_q, q_another_mod_p, q_another_mod_q,
                      NULL);
            mpz_gcd(p_another_gcd_q_another, p_another, q_another);
            mpz_mod(p_another_mod_p, p_another, p);
            mpz_mod(p_another_mod_q, p_another, q);
            mpz_mod(q_another_mod_p, q_another, p);
            mpz_mod(q_another_mod_q, q_another, q);

            int res = mpz_cmp_ui(p_another_gcd_q_another, 1) == 0 &&
                      mpz_cmp_ui(p_another_mod_p, 0) != 0 &&
                      mpz_cmp_ui(p_another_mod_q, 0) != 0 &&
                      mpz_cmp_ui(q_another_mod_p, 0) != 0 &&
                      mpz_cmp_ui(q_another_mod_q, 0) != 0;
//            int p_q = areCoprime(p, q);
//            int p_p_another = areCoprime(p, p_another);
//            int p_q_another = areCoprime(p, q_another);
//            int q_p_another = areCoprime(q, p_another);
//            int q_q_another = areCoprime(q, q_another);
//            int p_another_q_another = areCoprime(p_another, q_another);
//            printf("p_q = %d\n", p_q);
//            printf("p_p_another = %d\n", p_p_another);
//            printf("p_q_another = %d\n", p_q_another);
//            printf("q_p_another = %d\n", q_p_another);
//            printf("q_q_another = %d\n", q_q_another);
//            printf("p_another_q_another = %d\n", p_another_q_another);
//            int res = p_q && p_p_another && p_q_another && q_p_another && q_q_another && p_another_q_another;
            mpz_clears(p_another_gcd_q_another, p_another_mod_p, p_another_mod_q, q_another_mod_p, q_another_mod_q,
                       NULL);
            return res;
        }

    private:
        void get_odd_len_integer(mpz_t dst, int len)
        {
            mpz_rrandomb(dst, randstate, len);
            mpz_setbit(dst, 0); // 确保最低位为1，以保证生成的数是奇数
        }

        mpz_t *Ngen(int k)
        {
            int n_k = getN_k(k);
            int l_k = getL_k(k);
//            printf("n_k = %d\n", n_k);
//            printf("l_k = %d\n", l_k);
            mpz_t p, q, P, Q;
            mpz_inits(p, q, P, Q, NULL);

//            mpz_rrandomb(p, randstate, l_k / 2);//P: l_k / 2 bit
//            get_odd_len_integer(p, l_k / 2);
//            get_odd_len_integer(q, l_k / 2);
//            while (mpz_cmp(p, q) == 0) {
//                get_odd_len_integer(q, l_k / 2);
//            }
//            printf("p 的比特长度为 %ld bits\n", mpz_sizeinbase(p, 2));
//            printf("q 的比特长度为 %ld bits\n", mpz_sizeinbase(q, 2));

            do
            {
                mpz_rrandomb(p, randstate, l_k / 2);
            } while (mpz_probab_prime_p(p, 80) == 0 || mpz_sizeinbase(p, 2) != l_k / 2);

            do
            {
                mpz_rrandomb(q, randstate, l_k / 2);
            } while (mpz_probab_prime_p(q, 80) == 0 || mpz_cmp(p, q) == 0 || mpz_sizeinbase(q, 2) != l_k / 2);
//            printf("p 的比特长度为 %ld bits\n", mpz_sizeinbase(p, 2));
//            printf("q 的比特长度为 %ld bits\n", mpz_sizeinbase(q, 2));
            do
            {
                mpz_t p_another, q_another;
                mpz_inits(p_another, q_another, NULL);
                do
                {
                    get_odd_len_integer(p_another, (n_k - l_k) / 2 - 1);
                } while (mpz_sizeinbase(p_another, 2) != (n_k - l_k) / 2 - 1);
                do
                {
                    get_odd_len_integer(q_another, (n_k - l_k) / 2 - 1);
                } while (mpz_cmp(p_another, q_another) == 0 || mpz_sizeinbase(q_another, 2) != (n_k - l_k) / 2 - 1);
//                printf("p_another 的比特长度为 %ld bits\n", mpz_sizeinbase(p_another, 2));
//                printf("q_another 的比特长度为 %ld bits\n", mpz_sizeinbase(q_another, 2));

                if (isCoPrime(p, q, p_another, q_another))
                {
                    mpz_mul_ui(P, p, 2);
                    mpz_mul(P, P, p_another);
                    mpz_add_ui(P, P, 1);
                    mpz_mul_ui(Q, q, 2);
                    mpz_mul(Q, Q, q_another);
                    mpz_add_ui(Q, Q, 1);
                    bool P_is_prime = mpz_probab_prime_p(P, 80); // 判断是不是素数
                    bool Q_is_prime = mpz_probab_prime_p(Q, 80); // 判断是不是素数

                    // case1
                    if (P_is_prime && !Q_is_prime)
                    {
                        while (true)
                        {
                            get_odd_len_integer(q_another, (n_k - l_k) / 2 - 1);
                            mpz_t p_another_gcd_q_another, q_another_mod_p, q_another_mod_q;
                            mpz_inits(p_another_gcd_q_another, q_another_mod_p, q_another_mod_q, NULL);
                            mpz_gcd(p_another_gcd_q_another, p_another, q_another);
                            mpz_mod(q_another_mod_p, q_another, p);
                            mpz_mod(q_another_mod_q, q_another, q);
                            if (mpz_cmp_ui(p_another_gcd_q_another, 1) == 0 &&
                                mpz_cmp_ui(q_another_mod_p, 0) != 0 &&
                                mpz_cmp_ui(q_another_mod_q, 0) != 0)
                            {
                                mpz_mul_ui(Q, q, 2);
                                mpz_mul(Q, Q, q_another);
                                mpz_add_ui(Q, Q, 1);
                                Q_is_prime = mpz_probab_prime_p(Q, 80);
                                if (Q_is_prime) break;
                            }
                        }
                    }
                    // case2
                    if (Q_is_prime && !P_is_prime)
                    {
                        while (true)
                        {
                            get_odd_len_integer(p_another, (n_k - l_k) / 2 - 1);
                            mpz_t p_another_gcd_q_another, p_another_mod_p, p_another_mod_q;
                            mpz_inits(p_another_gcd_q_another, p_another_mod_p, p_another_mod_q, NULL);
                            mpz_gcd(p_another_gcd_q_another, p_another, q_another);
                            mpz_mod(p_another_mod_p, p_another, p);
                            mpz_mod(p_another_mod_q, p_another, q);
                            if (mpz_cmp_ui(p_another_gcd_q_another, 1) == 0 &&
                                mpz_cmp_ui(p_another_mod_p, 0) != 0 &&
                                mpz_cmp_ui(p_another_mod_q, 0) != 0)
                            {
                                mpz_mul_ui(P, p, 2);
                                mpz_mul(P, P, p_another);
                                mpz_add_ui(P, P, 1);
                                P_is_prime = mpz_probab_prime_p(P, 80);
                                if (P_is_prime) break;
                            }
                        }
                    }
                    // case3
                    if (P_is_prime && Q_is_prime)
                    {
                        break;
                    }

                }
            } while (true);

            mpz_t N;
            mpz_init(N);
            mpz_mul(N, P, Q);
//            printf("P 的比特长度为 %ld bits\n", mpz_sizeinbase(P, 2));
//            printf("Q 的比特长度为 %ld bits\n", mpz_sizeinbase(Q, 2));
//            printf("N 的比特长度为 %ld bits\n", mpz_sizeinbase(N, 2));

            auto *res = new mpz_t[5];
            mpz_inits(res[0], res[1], res[2], res[3], res[4], NULL);
            mpz_set(res[0], N);
            mpz_set(res[1], P);
            mpz_set(res[2], Q);
            mpz_set(res[3], p);
            mpz_set(res[4], q);

            mpz_clears(p, q, P, Q, N, NULL);
            return res;
        }
    };
}

#endif //FASTPAI_H
