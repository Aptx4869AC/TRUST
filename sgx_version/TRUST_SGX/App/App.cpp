/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-5-27
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#include <iostream>
#include <stdio.h>
#include <string.h>
#include <assert.h>
#include <unistd.h>
#include <pwd.h>
#include <omp.h>
#include <gmp.h>
#include <vector>
#include <random>
#include <ctime>
#include <cstring>
#include <sys/socket.h>
#include <sys/types.h>
#include <netinet/in.h>
#include <arpa/inet.h>
#include <iomanip>
#include <cmath>
#include <algorithm>
#include <netinet/tcp.h>
#include <fcntl.h>
#include <sys/uio.h>
#include <errno.h>
#include <fstream>
#include <sstream>
#include "sgx_urts.h"
#include "App.h"
#include "Enclave_u.h"
#include "fastPai.h"
#include "offline.h"

using namespace std;
using namespace PHESPACE;
using namespace OFFLINESPACE;

sgx_enclave_id_t global_eid = 0;
int sigma;
static ssize_t total_bytes;

typedef struct _sgx_errlist_t
{
    sgx_status_t err;
    const char *msg;
    const char *sug; /* Suggestion */
} sgx_errlist_t;

/* Error code returned by sgx_create_enclave */
static sgx_errlist_t sgx_errlist[] =
        {
                {
                        SGX_ERROR_UNEXPECTED, "Unexpected error occurred.", NULL
                }, {
                SGX_ERROR_INVALID_PARAMETER, "Invalid parameter.", NULL
        }, {
                SGX_ERROR_OUT_OF_MEMORY, "Out of memory.", NULL
        }, {
                SGX_ERROR_ENCLAVE_LOST, "Power transition occurred.", "Please refer to the sample \"PowerTransition\" for details."
        }, {
                SGX_ERROR_INVALID_ENCLAVE, "Invalid enclave image.", NULL
        }, {
                SGX_ERROR_INVALID_ENCLAVE_ID, "Invalid enclave identification.", NULL
        }, {
                SGX_ERROR_INVALID_SIGNATURE, "Invalid enclave signature.", NULL
        }, {
                SGX_ERROR_OUT_OF_EPC, "Out of EPC memory.", NULL
        }, {
                SGX_ERROR_NO_DEVICE, "Invalid SGX device.", "Please make sure SGX module is enabled in the BIOS, and install SGX driver afterwards."
        }, {
                SGX_ERROR_MEMORY_MAP_CONFLICT, "Memory map conflicted.", NULL
        }, {
                SGX_ERROR_INVALID_METADATA, "Invalid enclave metadata.", NULL
        }, {
                SGX_ERROR_DEVICE_BUSY, "SGX device was busy.", NULL
        }, {
                SGX_ERROR_INVALID_VERSION, "Enclave version was invalid.", NULL
        }, {
                SGX_ERROR_INVALID_ATTRIBUTE, "Enclave was not authorized.", NULL
        }, {
                SGX_ERROR_ENCLAVE_FILE_ACCESS, "Can't open enclave file.", NULL
        }, {
                SGX_ERROR_MEMORY_MAP_FAILURE, "Failed to reserve memory for the enclave.", NULL
        },};


/* Check error conditions for loading enclave */
void print_error_message(sgx_status_t ret)
{
    size_t idx = 0;
    size_t ttl = sizeof sgx_errlist / sizeof sgx_errlist[0];

    for (idx = 0; idx < ttl; idx++)
    {
        if (ret == sgx_errlist[idx].err)
        {
            if (NULL != sgx_errlist[idx].sug)
                printf("Info: %s\n", sgx_errlist[idx].sug);
            printf("Error: %s\n", sgx_errlist[idx].msg);
            break;
        }
    }

    if (idx == ttl)
        printf("Error code is 0x%X. Please refer to the \"Intel SGX SDK Developer Reference\" for more details.\n",
               ret);
}

/* Initialize the enclave:
 *   Call sgx_create_enclave to initialize an enclave instance
 */
int initialize_enclave(void)
{
    sgx_status_t ret = SGX_ERROR_UNEXPECTED;

    /* Call sgx_create_enclave to initialize an enclave instance */
    /* Debug Support: set 2nd parameter to 1 */
    ret = sgx_create_enclave(ENCLAVE_FILENAME, SGX_DEBUG_FLAG, NULL, NULL, &global_eid, NULL);
    if (ret != SGX_SUCCESS)
    {
        print_error_message(ret);
        return -1;
    }

    return 0;
}


/* Other functions */

ssize_t save_to_file(const char *filename, const vector<mpz_ptr> &src)
{
    // Open the file for writing
    int fd = open(filename, O_WRONLY | O_CREAT | O_TRUNC, 0644);
    if (fd == -1)
    {
        perror("open");
        exit(EXIT_FAILURE);
    }
    vector<struct iovec> iov(src.size());
    for (size_t i = 0; i < src.size(); ++i)
    {
        size_t n = mpz_sizeinbase(src[i], 10);
        unsigned char *msg = new unsigned char[n];
        mpz_export(msg, &n, 1, 1, 0, 0, src[i]);

        iov[i].iov_base = msg;
        iov[i].iov_len = n;

    }
    ssize_t bytesWritten = writev(fd, iov.data(), src.size());
    if (bytesWritten == -1)
    {
        perror("writev");
        exit(EXIT_FAILURE);
    }
//        printf("%ld bytes written to file '%s'.\n", bytesWritten, filename);

    for (size_t i = 0; i < src.size(); ++i)
    {
        delete[] static_cast<char *>(iov[i].iov_base);
    }

    // Close the file descriptor
    close(fd);

    return bytesWritten;
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


/**
 * 序列化
 * @param num
 * @param buffer
 * @param size
 */
void serialize_mpz(mpz_t num, unsigned char **buffer, size_t *size)
{
    *size = (mpz_sizeinbase(num, 2) + 7) / 8; // 计算所需字节数
    *buffer = new unsigned char[*size];
    mpz_export(*buffer, size, 1, 1, 0, 0, num); // 导出为字节数组
}

/**
 * 反序列化
 * @param res
 * @param buffer
 * @param size
 */
void deserialize_mpz(mpz_t res, unsigned char *buffer, size_t size)
{
    mpz_import(res, size, 1, 1, 0, 0, buffer); // 从字节数组导入
}

void test_serialize_deserialize()
{
    mpz_t original_var, deserialized_var;
    mpz_inits(original_var, deserialized_var, NULL);

    // Assign a value to original_var (replace with your value)
    mpz_set_ui(original_var, 123456789);

    // Serialize original_var
    unsigned char *buffer;
    size_t size;
    serialize_mpz(original_var, &buffer, &size);

    // Deserialize into deserialized_var
    deserialize_mpz(deserialized_var, buffer, size);

    // Print original_var and deserialized_var
    gmp_printf("Original variable: %Zd\n", original_var);
    gmp_printf("Deserialized variable: %Zd\n", deserialized_var);

    // Free memory allocated for buffer
    free(buffer);

    // Clear mpz variables
    mpz_clears(original_var, deserialized_var, NULL);
}


/* OCall functions */
void ocall_print_string(const char *str)
{
    printf("%s", str);
}

/* ECall functions */
/**
 * 密文乘
 * @param res
 * @param em_1
 * @param em_2
 * @param cp
 * @param record_bytes_flag
 */
void REE_Trust_FMUL(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, bool record_bytes_flag)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

    // Step 1, REE computes
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
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(X);
        src_vec.push_back(X_PDec_1);
        src_vec.push_back(Y);
        src_vec.push_back(em_2);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    /* send <[m_1 + r], (PDec)[m_1 + r], [m_2], [-r * m_2]> to csp */
    unsigned char *seal_res, *seal_X, *seal_X_PDec_1, *seal_Y, *seal_em_2;
    size_t res_size, X_size, X_PDec_1_size, Y_size, em_2_size;
    serialize_mpz(X, &seal_X, &X_size);
    serialize_mpz(X_PDec_1, &seal_X_PDec_1, &X_PDec_1_size);
    serialize_mpz(Y, &seal_Y, &Y_size);
    serialize_mpz(em_2, &seal_em_2, &em_2_size);
    res_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_res = new unsigned char[res_size];
    Enclave_Trust_FMUL(global_eid, &seal_res, &res_size, &seal_X, &X_size, &seal_X_PDec_1, &X_PDec_1_size,
                       &seal_Y, &Y_size, &seal_em_2, &em_2_size);


    deserialize_mpz(res, seal_res, res_size);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(res);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    mpz_clears(r, neg_r, er, X, Y, X_PDec_1, NULL);
    delete[] seal_res;
    delete[] seal_X;
    delete[] seal_X_PDec_1;
    delete[] seal_Y;
    delete[] seal_em_2;


}

/**
 * 密文比较
 * @param res
 * @param em_1
 * @param em_2
 * @param cp
 * @param record_bytes_flag
 */
void REE_Trust_FCMP(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, bool record_bytes_flag)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;
    mpz_t e_trust_r2, e_trust_r1_add_r2;
    mpz_inits(e_trust_r2, e_trust_r1_add_r2, NULL);

    mpz_t ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1;
    mpz_inits(ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, NULL);
    mpz_set(r1, tuple_S0.trust_r1);
    mpz_neg(neg_r1, r1);

    // Step 1
    int pi = (int) random() % 2;
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

    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(ed);
        src_vec.push_back(ed_1);
        src_vec.push_back(epi);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    /* send <[d], (PDec)[d], [pai]> to csp */
    unsigned char *seal_res, *seal_ed, *seal_ed_1, *seal_epi;
    size_t res_size, ed_size, ed_1_size, epi_size;
    serialize_mpz(ed, &seal_ed, &ed_size);
    serialize_mpz(ed_1, &seal_ed_1, &ed_1_size);
    serialize_mpz(epi, &seal_epi, &epi_size);
    res_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_res = new unsigned char[res_size];
    Enclave_Trust_FCMP(global_eid, &seal_res, &res_size, &seal_ed, &ed_size, &seal_ed_1, &ed_1_size,
                       &seal_epi, &epi_size);


    deserialize_mpz(res, seal_res, res_size);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(res);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed, ed_1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, epi, NULL);
    delete[] seal_res;
    delete[] seal_ed;
    delete[] seal_ed_1;
    delete[] seal_epi;
}

/**
 * 密文相等比较
 * @param res
 * @param em_1
 * @param em_2
 * @param cp
 * @param record_bytes_flag
 */
void REE_Trust_FEQL(mpz_t res, mpz_t em_1, mpz_t em_2, PaillierThdDec cp, bool record_bytes_flag)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;
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
    /* send <[d_1], (PDec)[d_1], [pi_1], [d_2], (PDec)[d_2], [pi_2]> to csp */
    unsigned char *seal_res, *seal_ed_1, *seal_ed_1_P1, *seal_epi_1;
    unsigned char *seal_ed_2, *seal_ed_2_P1, *seal_epi_2;
    size_t res_size, ed_1_size, ed_1_P1_size, epi_1_size;
    size_t ed_2_size, ed_2_P1_size, epi_2_size;
    serialize_mpz(ed_1, &seal_ed_1, &ed_1_size);
    serialize_mpz(ed_1_P1, &seal_ed_1_P1, &ed_1_P1_size);
    serialize_mpz(epi_1, &seal_epi_1, &epi_1_size);
    serialize_mpz(ed_2, &seal_ed_2, &ed_2_size);
    serialize_mpz(ed_2_P1, &seal_ed_2_P1, &ed_2_P1_size);
    serialize_mpz(epi_2, &seal_epi_2, &epi_2_size);
    res_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_res = new unsigned char[res_size];
    Enclave_Trust_FEQL(global_eid, &seal_res, &res_size,
                       &seal_ed_1, &ed_1_size, &seal_ed_1_P1, &ed_1_P1_size, &seal_epi_1, &epi_1_size,
                       &seal_ed_2, &ed_2_size, &seal_ed_2_P1, &ed_2_P1_size, &seal_epi_2, &epi_2_size);


    deserialize_mpz(res, seal_res, res_size);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(res);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed_1, ed_1_P1, em_1r1, em_2r1, r1, neg_r1, neg_em_2r1, neg_em_1r1, epi_1, NULL);
    mpz_clears(e_trust_r2_prime, e_trust_r1_prime_add_r2_prime, ed_2, ed_2_P1, em_1r1_prime, em_2r1_prime,
               r1_prime, neg_r1_prime, neg_em_2r1_prime, neg_em_1r1_prime, epi_2, NULL);
    delete[] seal_res;
    delete[] seal_ed_1;
    delete[] seal_ed_1_P1;
    delete[] seal_epi_1;
    delete[] seal_ed_2;
    delete[] seal_ed_2_P1;
    delete[] seal_epi_2;
}

/**
 * 密文绝对值
 * @param res
 * @param em_1
 * @param cp
 * @param record_bytes_flag
 */
void REE_Trust_FABS(mpz_t res, mpz_t em_1, PaillierThdDec cp, bool record_bytes_flag)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

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

    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(ed);
        src_vec.push_back(ed_1);
        src_vec.push_back(em_1pi);
        src_vec.push_back(em_1);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    /* send <[d], (PDec)[d], [pai * m_1], [m_1]> to csp */
    unsigned char *seal_res, *seal_ed, *seal_ed_1, *seal_em_1pi, *seal_em_1;
    size_t res_size, ed_size, ed_1_size, em_1pi_size, em_1_size;
    serialize_mpz(ed, &seal_ed, &ed_size);
    serialize_mpz(ed_1, &seal_ed_1, &ed_1_size);
    serialize_mpz(em_1pi, &seal_em_1pi, &em_1pi_size);
    serialize_mpz(em_1, &seal_em_1, &em_1_size);
    res_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_res = new unsigned char[res_size];
    Enclave_Trust_FABS(global_eid, &seal_res, &res_size,
                       &seal_ed, &ed_size, &seal_ed_1, &ed_1_size, &seal_em_1pi, &em_1pi_size,
                       &seal_em_1, &em_1_size);


    deserialize_mpz(res, seal_res, res_size);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(res);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, em_1pi, NULL);
    delete[] seal_res;
    delete[] seal_ed;
    delete[] seal_ed_1;
    delete[] seal_em_1pi;
    delete[] seal_em_1;

}

/**
 * 密文运算符 version1
 * @param res
 * @param em_1
 * @param em_2
 * @param em_3
 * @param cp
 * @param record_bytes_flag
 */
void REE_Trust_FTRN_v1(mpz_t res, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, bool record_bytes_flag)
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
    mpz_t epi_1;
    mpz_inits(epi_1, NULL);
    if (pi_1 == 0)
    {
        cp.pai.add(epi_1, tuple_S0.e0, tuple_S0.e0);   //refresh

    } else
    {
        cp.pai.add(epi_1, tuple_S0.e1, tuple_S0.e0);   //refresh
    }

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
    mpz_t epi_2;
    mpz_inits(epi_2, NULL);
    if (pi_2 == 0)
    {
        cp.pai.add(epi_2, tuple_S0.e0, tuple_S0.e0);   //refresh

    } else
    {
        cp.pai.add(epi_2, tuple_S0.e1, tuple_S0.e0);   //refresh
    }

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

    /* send <[m3-m2], [m2], [d_1], (PDec)[d_1], [pi_1 * (m3-m2)], [d_2], (PDec)[d_2], [pi_2 * (m3-m2)]> to csp */
    unsigned char *seal_res, *seal_em_3_sub_m_2, *seal_em_2;
    unsigned char *seal_ed_1, *seal_ed_1_P1, *seal_epi_1_m_3_sub_m_2;
    unsigned char *seal_ed_2, *seal_ed_2_P1, *seal_epi_2_m_3_sub_m_2;

    size_t res_size, em_3_sub_m_2_size, em_2_size;
    size_t ed_1_size, ed_1_P1_size, epi_1_m_3_sub_m_2_size;
    size_t ed_2_size, ed_2_P1_size, epi_2_m_3_sub_m_2_size;

    serialize_mpz(em_3_sub_m_2, &seal_em_3_sub_m_2, &em_3_sub_m_2_size);
    serialize_mpz(em_2, &seal_em_2, &em_2_size);

    serialize_mpz(ed_1, &seal_ed_1, &ed_1_size);
    serialize_mpz(ed_1_P1, &seal_ed_1_P1, &ed_1_P1_size);
    serialize_mpz(epi_1_m_3_sub_m_2, &seal_epi_1_m_3_sub_m_2, &epi_1_m_3_sub_m_2_size);

    serialize_mpz(ed_2, &seal_ed_2, &ed_2_size);
    serialize_mpz(ed_2_P1, &seal_ed_2_P1, &ed_2_P1_size);
    serialize_mpz(epi_2_m_3_sub_m_2, &seal_epi_2_m_3_sub_m_2, &epi_2_m_3_sub_m_2_size);

    res_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_res = new unsigned char[res_size];
    Enclave_Trust_FTRN_v1(global_eid, &seal_res, &res_size,
                          &seal_em_3_sub_m_2, &em_3_sub_m_2_size, &seal_em_2, &em_2_size,
                          &seal_ed_1, &ed_1_size, &seal_ed_1_P1, &ed_1_P1_size, &seal_epi_1_m_3_sub_m_2, &epi_1_m_3_sub_m_2_size,
                          &seal_ed_2, &ed_2_size, &seal_ed_2_P1, &ed_2_P1_size, &seal_epi_2_m_3_sub_m_2, &epi_2_m_3_sub_m_2_size);

    deserialize_mpz(res, seal_res, res_size);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(res);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed_1, ed_1_P1, em_1r1, r1, neg_r1, neg_em_1r1, epi_1, NULL);
    mpz_clears(e_trust_2r1_prime_add_r2_prime, e_trust_r2_prime_sub_r1_prime, ed_2, ed_2_P1, em_1r1_prime,
               r1_prime, neg_r1_prime, neg_em_1r1_prime, epi_2, NULL);
    mpz_clears(em_3_sub_m_2, epi_1_m_3_sub_m_2, epi_2_m_3_sub_m_2, NULL);
    delete[] seal_res;
    delete[] seal_em_3_sub_m_2;
    delete[] seal_em_2;
    delete[] seal_ed_1;
    delete[] seal_ed_1_P1;
    delete[] seal_epi_1_m_3_sub_m_2;
    delete[] seal_ed_2;
    delete[] seal_ed_2_P1;
    delete[] seal_epi_2_m_3_sub_m_2;

}


/**
 * 密文运算符 version2
 * @param res
 * @param em_1
 * @param em_2
 * @param em_3
 * @param cp
 * @param record_bytes_flag
 */
void REE_Trust_FTRN_v2(mpz_t res, mpz_t em_1, mpz_t em_2, mpz_t em_3, PaillierThdDec cp, bool record_bytes_flag)
{
    Offline_val tuple_S0 = cp.pai.offlineVal;

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

    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(em_3);
        src_vec.push_back(em_2);
        src_vec.push_back(ed);
        src_vec.push_back(ed_1);
        src_vec.push_back(epi_m_3_sub_m_2);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    /* send <[m3], [m2], [d_1], (PDec)[d_1], [pi * (m3-m2)] > to csp */
    unsigned char *seal_res, *seal_em_3, *seal_em_2, *seal_ed, *seal_ed_1, *seal_epi_m_3_sub_m_2;
    size_t res_size, em_3_size, em_2_size, ed_size, ed_1_size, epi_m_3_sub_m_2_size;
    serialize_mpz(em_3, &seal_em_3, &em_3_size);
    serialize_mpz(em_2, &seal_em_2, &em_2_size);
    serialize_mpz(ed, &seal_ed, &ed_size);
    serialize_mpz(ed_1, &seal_ed_1, &ed_1_size);
    serialize_mpz(epi_m_3_sub_m_2, &seal_epi_m_3_sub_m_2, &epi_m_3_sub_m_2_size);
    res_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
    seal_res = new unsigned char[res_size];
    Enclave_Trust_FTRN_v2(global_eid, &seal_res, &res_size,
                          &seal_em_3, &em_3_size, &seal_em_2, &em_2_size, &seal_ed, &ed_size,
                          &seal_ed_1, &ed_1_size, &seal_epi_m_3_sub_m_2, &epi_m_3_sub_m_2_size);


    deserialize_mpz(res, seal_res, res_size);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(res);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }

    mpz_clears(e_trust_r2, e_trust_r1_add_r2, ed, ed_1, em_1r1, r1, neg_r1, neg_em_1r1, em_3_sub_m_2, epi_m_3_sub_m_2, NULL);
    delete[] seal_res;
    delete[] seal_em_3;
    delete[] seal_em_2;
    delete[] seal_ed;
    delete[] seal_ed_1;
    delete[] seal_epi_m_3_sub_m_2;

}


/* Application entry */
int main(int argc, char *argv[])
{


    /* Initialize the enclave */
    if (initialize_enclave() < 0)
    {
        printf("Enter a character before exit ...\n");
        getchar();
        return -1;
    }

    double start_time, end_time;
    double average_time_keygen = 0;
    double average_time_enc = 0;
    double average_time_dec = 0;
    double average_time_PDec_TDec = 0;
    double average_time_HE_Addition = 0;
    double average_time_HE_ScalarMul = 0;
    double average_time_HE_Subtraction = 0;
    double average_time_Trust_FMUL = 0;
    double average_time_Trust_FCMP = 0;
    double average_time_Trust_FEQL = 0;
    double average_time_Trust_FABS = 0;
    double average_time_Trust_FTRN_v1 = 0;
    double average_time_Trust_FTRN_v2 = 0;

    double average_bytes_Trust_FMUL = 0;
    double average_bytes_Trust_FCMP = 0;
    double average_bytes_Trust_FEQL = 0;
    double average_bytes_Trust_FABS = 0;
    double average_bytes_Trust_FTRN_v1 = 0;
    double average_bytes_Trust_FTRN_v2 = 0;

    int k = 64; //BIT SECURITY, 64、80、104、112、124、128
    sigma = 128;
    switch (k)
    {
        case 64:
            printf("N 的比特长度为 512 bits\n");
            break;
        case 80:
            printf("N 的比特长度为 1024 bits\n");
            break;
        case 104:
            printf("N 的比特长度为 1536 bits\n");
            break;
        case 112:
            printf("N 的比特长度为 2048 bits\n");
            break;
        case 124:
            printf("N 的比特长度为 2560 bits\n");
            break;
        case 128:
            printf("N 的比特长度为 3072 bits\n");
            break;
        case 192:
            printf("N 的比特长度为 7680 bits\n");
            break;
    }
    int epoch = 10;
//    printf("循环次数 epoch = %d\n", epoch);
    setrandom(&randstate);
    KeyGen keyGen(k, sigma);

    start_time = omp_get_wtime();
    for (int i = 0; i < epoch; i++)
    {
        KeyGen temp_keyGen(k, sigma);
    }
    end_time = omp_get_wtime();
    average_time_keygen += (end_time - start_time);

    Paillier pai(keyGen.pk, keyGen.sk, sigma);
    PaillierThirdParty *psk = keyGen.pai_third_parties;
    PaillierThdDec cp = PaillierThdDec(psk[0].N, psk[0].partial_key, keyGen.pk, sigma);
    printf("REE 初始化成功\n");

    // 序列化TEE部分的私钥内容
    PaillierThdDec enclave = PaillierThdDec(psk[1].N, psk[1].partial_key, keyGen.pk, sigma);
    unsigned char *seal_partial_key, *seal_n, *seal_h_N;
    size_t partial_key_size, n_size, h_N_size;
    serialize_mpz(enclave.partial_key, &seal_partial_key, &partial_key_size);
    serialize_mpz(enclave.N, &seal_n, &n_size);
    serialize_mpz(enclave.pk.h_N, &seal_h_N, &h_N_size);
    Enclave_Init(global_eid, &seal_partial_key, &partial_key_size, &seal_n, &n_size, &seal_h_N, &h_N_size,
                 enclave.pk.L_k);

    // test Basically Cryptographic
    for (int i = 0; i < epoch; i++)
    {
        mpz_t a, ea;
        mpz_inits(a, ea, NULL);
        mpz_rrandomb(a, randstate, 16);
//        gmp_printf("a = %Zd\n", a);

        start_time = omp_get_wtime();
        pai.encrypt(ea, a);
        end_time = omp_get_wtime();
        average_time_enc += (end_time - start_time);
//        printf("%d [Enc] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        start_time = omp_get_wtime();
        pai.decrypt(a, ea);
        end_time = omp_get_wtime();
        average_time_dec += (end_time - start_time);
//        printf("%d [Dec] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t a1, a2, a_cp, a_csp;
        mpz_inits(a1, a2, a_cp, a_csp, NULL);
        start_time = omp_get_wtime();
        cp.PDec(a1, ea);
        unsigned char *seal_pa;
        size_t pa_size;
        unsigned char *seal_a;
        size_t a_size;
        serialize_mpz(ea, &seal_a, &a_size);
        pa_size = (mpz_sizeinbase(cp.pai.pk.N_Square, 2) + 7) / 8; // 计算所需字节数
        seal_pa = new unsigned char[pa_size];

        Enclave_PDec(global_eid, &seal_pa, &pa_size, &seal_a, &a_size);
        deserialize_mpz(a2, seal_pa, pa_size);
        cp.TDec(a_cp, a1, a2);
        end_time = omp_get_wtime();
        if (mpz_cmp(a_cp, a) != 0)
        {
            gmp_printf("a_cp = %Zd\n", a_cp);
            gmp_printf("a = %Zd\n", a);
            printf("(Trust) Dec error!\n");
            printf("i = %d is error \n", i);
            exit(-1);
        }
        free(seal_pa);
        free(seal_a);
        average_time_PDec_TDec += (end_time - start_time);
//        printf("%d [PDec and TDec] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t cz;
        mpz_init(cz);
        start_time = omp_get_wtime();
        mpz_mul(cz, ea, ea);
        mpz_mod(cz, cz, pai.pk.N_Square);
        end_time = omp_get_wtime();
        average_time_HE_Addition += (end_time - start_time);
//        printf("%d [HE_Addition] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t ten;
        mpz_init(ten);
        mpz_set_si(ten, 10);
        start_time = omp_get_wtime();
        mpz_powm(cz, ea, ten, pai.pk.N_Square);
        end_time = omp_get_wtime();
        average_time_HE_ScalarMul += (end_time - start_time);
//        printf("%d [HE_ScalarMul] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_t neg_one;
        mpz_init(neg_one);
        mpz_set_si(neg_one, -1);
        start_time = omp_get_wtime();
        mpz_powm(cz, ea, neg_one, pai.pk.N_Square);
        mpz_mul(cz, ea, ea);
        mpz_mod(cz, cz, pai.pk.N_Square);
        end_time = omp_get_wtime();
        average_time_HE_Subtraction += (end_time - start_time);
//        printf("%d [HE_Subtraction] time is  ------  %f ms\n", i, (end_time - start_time) * 1000);

        mpz_clears(a, ea, a1, a2, a_cp, a_csp, NULL);
//        printf("----------------------------------------------------------\n");
    }

    printf("----------------------------------------------------------\n");
    printf("KeyGen (%d average) time is  ------  %f ms\n", epoch, average_time_keygen / epoch * 1000);
    printf("Enc (%d average) time is  ------  %f ms\n", epoch, average_time_enc / epoch * 1000);
    printf("Dec (%d average) time is  ------  %f ms\n", epoch, average_time_dec / epoch * 1000);
    printf("PDec + TDec (%d average) time is  ------  %f ms\n", epoch, average_time_PDec_TDec / epoch * 1000);
    printf("HE_Addition (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Addition / epoch * 1000);
    printf("HE_ScalarMul (%d average) time is  ------  %f ms\n", epoch, average_time_HE_ScalarMul / epoch * 1000);
    printf("HE_Subtraction (%d average) time is  ------  %f ms\n", epoch, average_time_HE_Subtraction / epoch * 1000);
    printf("----------------------------------------------------------\n");


    ifstream file("Data/training_dataADD.csv");  // 打开CSV文件
    string line, value;

    if (!file.is_open())
    {
        cerr << "无法打开文件 training_dataADD.csv" << '\n';
        return 1;
    }

    // 跳过表头
    getline(file, line);

    int _epoch = 1;
    bool record_bytes_flag = true;
    // test Computation Costs / Communication Costs
    // 逐行读取CSV文件的内容
    while (getline(file, line))
    {
        stringstream ss(line);
        // 从每一行读取 Value1 和 Value2
        getline(ss, value, ',');
        int value1 = stoi(value);  // 转换为整数

        getline(ss, value, ',');
        int value2 = stoi(value);  // 转换为整数
//        cout << "epoch: " << _epoch << ", Value1: " << value1 << ", Value2: " << value2 << '\n';

        mpz_t m_1, m_2, m_3, em_1, em_2, em_3, eres, res;
        mpz_inits(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
        mpz_set_si(m_1, value1);
        mpz_set_si(m_2, value2);
        if (mpz_cmp_si(m_1, 0) < 0)
        {
            mpz_add(m_1, m_1, cp.pai.pk.N);
        }
        if (mpz_cmp_si(m_2, 0) < 0)
        {
            mpz_add(m_2, m_2, cp.pai.pk.N);
        }
        mpz_rrandomb(m_3, randstate, 8);
//        mpz_rrandomb(m_1, randstate, 8); // 8-bit random number
//        mpz_rrandomb(m_2, randstate, 8); // 8-bit random number
//        mpz_rrandomb(m_3, randstate, 8); // 8-bit random number
        cp.pai.encrypt(em_1, m_1);
        cp.pai.encrypt(em_2, m_2);
        cp.pai.encrypt(em_3, m_3);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        REE_Trust_FMUL(eres, em_1, em_2, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FMUL += end_time - start_time;
        check_fmul(_epoch, eres, m_1, m_2, pai);
        average_bytes_Trust_FMUL += total_bytes / 1024.0;
//        printf("%d [FMUL] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        REE_Trust_FCMP(eres, em_1, em_2, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FCMP += end_time - start_time;
        check_fcmp(_epoch, eres, m_1, m_2, pai);
        average_bytes_Trust_FCMP += total_bytes / 1024.0;
//        printf("%d [FCMP] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        REE_Trust_FEQL(eres, em_1, em_2, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FEQL += end_time - start_time;
        check_feql(_epoch, eres, m_1, m_2, pai);
        average_bytes_Trust_FEQL += total_bytes / 1024.0;
//        printf("%d [FEQL] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        total_bytes = 0.0;
        start_time = omp_get_wtime();
        REE_Trust_FABS(eres, em_1, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FABS += (end_time - start_time);
        check_fabs(_epoch, eres, m_1, pai);
        average_bytes_Trust_FABS += total_bytes / 1024.0;
//        printf("%d [FABS] time is  ------  %f ms\n", _epoch + 1, (end_time - start_time) * 1000);


        total_bytes = 0.0;
        start_time = omp_get_wtime();
        REE_Trust_FTRN_v1(eres, em_1, em_2, em_3, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FTRN_v1 += end_time - start_time;
        check_ftrn(_epoch, eres, m_1, m_2, m_3, pai);
        average_bytes_Trust_FTRN_v1 += total_bytes / 1024.0;
//        printf("%d [FTRN] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        int pi = rand() % 2;
        mpz_set_si(m_1, pi);
        pai.encrypt(em_1, m_1);
        total_bytes = 0.0;
        start_time = omp_get_wtime();
        REE_Trust_FTRN_v2(eres, em_1, em_2, em_3, cp, record_bytes_flag);
        end_time = omp_get_wtime();
        average_time_Trust_FTRN_v2 += end_time - start_time;
        check_ftrn(_epoch, eres, m_1, m_2, m_3, pai);
        average_bytes_Trust_FTRN_v2 += total_bytes / 1024.0;
//        printf("%d [FTRN] time is %f ms\n", _epoch + 1, (end_time - start_time) * 1000);

        _epoch++;
        if (_epoch > 200)
        {
            break;
        }
        mpz_clears(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
//        printf("----------------------------------------------------------\n");
    }


    printf("Trust_FMUL (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FMUL / _epoch * 1000);
    printf("Trust_FCMP (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FCMP / _epoch * 1000);
    printf("Trust_FEQL (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FEQL / _epoch * 1000);
    printf("Trust_FABS (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FABS / _epoch * 1000);
    printf("Trust_FTRN_v1 (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FTRN_v1 / _epoch * 1000);
    printf("Trust_FTRN_v2 (%d average) time is  ------  %f ms\n", _epoch, average_time_Trust_FTRN_v2 / _epoch * 1000);

    if (record_bytes_flag)
    {
        printf("Trust_FMUL (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FMUL / _epoch);
        printf("Trust_FCMP (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FCMP / _epoch);
        printf("Trust_FEQL (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FEQL / _epoch);
        printf("Trust_FABS (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FABS / _epoch);
        printf("Trust_FTRN_v1 (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FTRN_v1 / _epoch);
        printf("Trust_FTRN_v2 (%d average) Communication Costs is  ---  %f KB\n", _epoch, average_bytes_Trust_FTRN_v2 / _epoch);
    }

    Enclave_Release(global_eid);

    /* Destroy the enclave */
    sgx_destroy_enclave(global_eid);

    file.close();  // 关闭文件

    printf("Info: successfully returned.\n");
    return 0;
}

