/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include <iostream>
#include <gmp.h>
#include <omp.h>
#include <ctime>
#include <unistd.h>
#include "include/paillier.h"
#include "include/soci.h"
#include <fstream>
#include <vector>
#include <fcntl.h>
#include <unistd.h>
#include <sys/uio.h>

using namespace std;
using namespace PHESPACE;
using namespace SOCISPACE;



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

int main() {

    double start_time, end_time;
    double average_time_keygen = 0;
    double average_time_enc = 0;
    double average_time_dec = 0;

    double average_time_HE_Addition = 0;
    double average_time_HE_ScalarMul = 0;
    double average_time_HE_Subtraction = 0;


    double average_time_PDec_TDec = 0;
    double average_time_SMUL = 0;
    double average_time_SCMP = 0;
    double average_time_SSBA = 0;
    double average_time_SDIV = 0;
    double average_bytes_PK = 0;
    double average_bytes_SK = 0;

    int epoch = 20;
    printf("epoch = %d\n", epoch);

    mpz_t x, y, z, c1, c2, cx, cy, cz, temp;
    mpz_inits(x, y, z, c1, c2, cx, cy, cz, temp, NULL);
    gmp_randstate_t randstate;
    setrandom(&randstate);
    mpz_rrandomb(x, randstate, 8); // 8-bit random number
    mpz_rrandomb(y, randstate, 8); // 8-bit random number
    gmp_printf("x = %Zd\n", x);
    gmp_printf("y = %Zd\n", y);

    soci sc;
    Paillier pai;
    int key_len = 1024; // 256, 512, 768, 1024
    pai.keygen(key_len);

    printf("key_len = %d\n", key_len);
    printf("N 的比特长度为 %d bits\n", 2 * key_len);

    printf("----------------------------------------------------------\n");
    start_time = omp_get_wtime();
    for (int i = 0; i < 10; i++) {
        Paillier temp;
        temp.keygen(key_len);
    }
    end_time = omp_get_wtime();
    average_time_keygen += (end_time - start_time);
    printf("Keygen (%d average) time is  ------  %f ms\n", epoch, average_time_keygen / 10 * 1000);

    PaillierThdPrivateKey *psk = thdkeygen(pai.prikey, sigma);
    PaillierThdDec cp(psk[0], pai.pubkey);
    PaillierThdDec csp(psk[1], pai.pubkey);
    
//    start_time = omp_get_wtime();
//    for (int i = 0; i < 10; i++) {
//        PaillierThdPrivateKey *temp_psk = thdkeygen(pai.prikey, sigma);
//        PaillierThdDec temp_cp(temp_psk[0], pai.pubkey);
//        PaillierThdDec temp_csp(temp_psk[1], pai.pubkey);
//    }
//    end_time = omp_get_wtime();
//    printf("thdkeygen (%d average) time is  ------  %f ms\n", epoch, (end_time - start_time) / 10 * 1000);

    pai.encrypt(cx, x);
    pai.encrypt(cy, y);

    printf("----------------------------------------------------------\n");
    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        start_time = omp_get_wtime();
        pai.encrypt(cz, cz);
        end_time = omp_get_wtime();
        average_time_enc += (end_time - start_time);
        mpz_clear(cz);

    }
    printf("Enc (%d average) time is  ------  %f ms\n", epoch, average_time_enc / epoch * 1000);


    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        start_time = omp_get_wtime();
        pai.decrypt(cz, cx);
        end_time = omp_get_wtime();
        average_time_dec += (end_time - start_time);
        if (mpz_cmp(cz, x) != 0) {
            gmp_printf("x = %Zd\n", y);
            gmp_printf("dec(x) = %Zd\n", cz);
            printf("Dec is error!\n");
            exit(-1);
        }
        mpz_clear(cz);
    }
    printf("Dec (%d average) time is  ------  %f ms\n", epoch, average_time_dec / epoch * 1000);


    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        //pai.decrypt(cz, cx);
        start_time = omp_get_wtime();
        cp.pdec(c1, cy);
        csp.pdec(c2, cy);
        cp.fdec(cz, c1, c2);
        end_time = omp_get_wtime();
        average_time_PDec_TDec += (end_time - start_time);
        if (mpz_cmp(cz, y) != 0) {
            gmp_printf("y = %Zd\n", y);
            gmp_printf("dec(y) = %Zd\n", cz);
            printf("PDec + TDec is error!\n");
            exit(-1);
        }
        mpz_clear(cz);

    }


    printf("PDec + TDec (%d average) time is  ------  %f ms\n", epoch, average_time_PDec_TDec / epoch * 1000);

    if (mpz_cmp_si(y, 0) < 0) {
        gmp_printf("y = %Zd\n", y);
        mpz_add(y, y, pai.pubkey.n);
        gmp_printf("y = %Zd (mod n)\n", y);
    }


    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        // x+y
        start_time = omp_get_wtime();
        mpz_mul(cz, cx, cy);
        mpz_mod(cz, cz, pai.pubkey.nsquare);
        end_time = omp_get_wtime();
        average_time_HE_Addition += (end_time - start_time);
        mpz_clear(cz);
    }

    printf("HE_Addition (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Addition / epoch * 1000);

    mpz_t ten;
    mpz_init(ten);
    mpz_set_si(ten, 10);
    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        // 10x
        start_time = omp_get_wtime();
        mpz_powm(cz, cx, ten, pai.pubkey.nsquare);
        end_time = omp_get_wtime();
        average_time_HE_ScalarMul += (end_time - start_time);
        mpz_clear(cz);
    }

    printf("HE_ScalarMul (%d average) time is  ------  %f ms\n", epoch, average_time_HE_ScalarMul / epoch * 1000);

    mpz_t neg_one;
    mpz_init(neg_one);
    mpz_set_si(neg_one, -1);
    for (int i = 0; i < epoch; i++) {
        mpz_t cz;
        mpz_init(cz);
        // x-y
        start_time = omp_get_wtime();
        mpz_powm(cz, cy, neg_one, pai.pubkey.nsquare);
        mpz_mul(cz, cx, cy);
        mpz_mod(cz, cz, pai.pubkey.nsquare);
        end_time = omp_get_wtime();
        average_time_HE_Subtraction += (end_time - start_time);
        mpz_clear(cz);
    }

    printf("HE_Subtraction (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Subtraction / epoch * 1000);
    
    for (int i = 0; i < epoch; i++) {
        mpz_t cz, temp, z;
        mpz_inits(z, cz, temp, NULL);
        start_time = omp_get_wtime();
        sc.smul(cz, cx, cy, cp, csp);
        end_time = omp_get_wtime();
        average_time_SMUL += (end_time - start_time);
        mpz_mul(temp, x, y);
        pai.decrypt(z, cz);
        if (mpz_cmp(z, temp) != 0) {
            gmp_printf("z = %Zd\n", z);
            gmp_printf("%Zd * %Zd = %Zd\n", x, y, temp);
            printf("SMUL is error!\n");
            exit(-1);
        }
        mpz_clears(cz, temp, z, NULL);
    }

    for (int i = 0; i < epoch; i++) {
        mpz_t z, cz, temp;
        mpz_inits(z, cz, temp, NULL);
        start_time = omp_get_wtime();
        sc.scmp(cz, cx, cy, cp, csp);
        end_time = omp_get_wtime();
        average_time_SCMP += (end_time - start_time);

        pai.decrypt(z, cz);
        if (mpz_cmp(x, y) >= 0) {
            mpz_set_si(temp, 0);
        } else mpz_set_si(temp, 1);
        if (mpz_cmp(z, temp) != 0) {

            gmp_printf("dec answer x>=y? = %Zd\n", z);
            gmp_printf("right answer %Zd\n", temp);
            printf("SCMP is error!\n");
            exit(-1);
        }
        mpz_clears(z, cz, temp, NULL);
    }


    for (int i = 0; i < epoch; i++) {
        mpz_t s_x, u_x, s, u, temp_s, temp_u;
        mpz_inits(s_x, u_x, s, u, temp_s, temp_u, NULL);
        start_time = omp_get_wtime();
        sc.ssba(s_x, u_x, cx, cp, csp);
        end_time = omp_get_wtime();
        average_time_SSBA += (end_time - start_time);

        pai.decrypt(s, s_x);
        pai.decrypt(u, u_x);
        if (mpz_cmp_si(cx, 0) >= 0) {
            mpz_set_si(temp_s, 0);
            mpz_set(temp_u, x);
        } else {
            mpz_set_si(temp_s, 1);
            mpz_set(temp_u, x);
            mpz_neg(temp_u, temp_u);// temp_u = -temp_u
        }

        if (mpz_cmp(temp_u, u) != 0 && mpz_cmp(s, s) != 0) {
            gmp_printf("dec answer s = %Zd u = %Zd\n", s, u);
            gmp_printf("right answer s = %Zd u = %Zd\n", temp_s, temp_u);
            printf("SSBA is error!\n");
            exit(-1);
        }
        mpz_clears(s_x, u_x, s, u, temp_s, temp_u, NULL);
    }

    printf("----------------------------------------------------------\n");

    int bit_length = 10;
    for (int i = 0; i < epoch; i++) {
        mpz_t eq, er, q, r, right_q, right_r;
        mpz_inits(eq, er, q, r, right_q, right_r, NULL);
        start_time = omp_get_wtime();
        sc.sdiv(eq, er, cx, cy, bit_length, cp, csp);
        end_time = omp_get_wtime();
        average_time_SDIV += (end_time - start_time);

        pai.decrypt(q, eq);
        pai.decrypt(r, er);
        mpz_div(right_q, x, y);
        mpz_mod(right_r, x, y);
        if (mpz_cmp(q, right_q) != 0 && mpz_cmp(r, right_r) != 0) {
            gmp_printf("q = %Zd r = %Zd\n", q, r);
            printf("SDIV is error!\n");
            exit(-1);
        }
        mpz_clears(eq, er, q, r, right_q, right_r, NULL);
    }
    printf("SMUL (%d average) time is  ------  %f ms\n", epoch, average_time_SMUL / epoch * 1000);
    printf("SCMP (%d average) time is  ------  %f ms\n", epoch, average_time_SCMP / epoch * 1000);
    printf("SSBA (%d average) time is  ------  %f ms\n", epoch, average_time_SSBA / epoch * 1000);
    printf("SDIV (%d-bits) (%d average) time is  ------  %f ms\n", epoch, bit_length, average_time_SDIV / epoch * 1000);


    return 0;
}

