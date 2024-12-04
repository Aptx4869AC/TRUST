/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#ifndef CP_PROTOCOL_H
#define CP_PROTOCOL_H
#pragma once

#include "pcpd.h"
#include "network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

ssize_t total_bytes = 0;
mpz_t zero, one, two, neg_one, N, N_square, N_half;

/**
 * POCF乘法
 * @param result
 * @param e_x
 * @param e_y
 * @param cp
 * @param record_bytes_flag
 */
void CP_SMUL(mpz_t result, const mpz_t &e_x, const mpz_t &e_y, const Paillier_Third_Party &cp, bool record_bytes_flag)
{
    // 正常步骤
    mpz_t flag, temp, e_r_x, e_r_y, X, Y, X1, Y1, h, H, S_1, S_2, S_3, e_x_y, r_x, r_y;
    mpz_inits(flag, temp, e_r_x, e_r_y, X, Y, X1, Y1, h, H, S_1, S_2, S_3, e_x_y, r_x, r_y, NULL);

    // step 1 (CP)
    // 选取两个随机数 r_x 和 r_y
    while (true)
    {
        mpz_urandomm(r_x, randstate, N);
        mpz_gcd(flag, r_x, N);
        if (mpz_cmp_ui(flag, 1) == 0)
        {
            break;
        }
    }
    while (true)
    {
        mpz_urandomm(r_y, randstate, N);
        mpz_gcd(flag, r_y, N);
        if (mpz_cmp_ui(flag, 1) == 0)
        {
            break;
        }
    }

    // 加密随机数 r_x 和 r_y
    Enc(e_r_x, cp.public_key, r_x);
    Enc(e_r_y, cp.public_key, r_y);

    // 计算 X = e_x * e_r_x mod N_square 和 Y = e_y * e_r_y mod N_square
    mpz_mul(X, e_x, e_r_x);
    mpz_mod(X, X, N_square);
    mpz_mul(Y, e_y, e_r_y);
    mpz_mod(Y, Y, N_square);
    // 部分解密 X 和 Y
    PDec(X1, cp.private_key, X);
    PDec(Y1, cp.private_key, Y);

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 1);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    // 发X,Y,X1,Y1给CSP
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

    // step 3 (CP)
    // 计算 S_1, S_2, S_3
    mpz_mul(temp, r_x, r_y);
    mpz_mod(temp, temp, N);
    Enc(S_1, cp.public_key, temp);
    mpz_powm(S_1, S_1, neg_one, N_square);

    mpz_neg(temp, r_y);
    mpz_powm(S_2, e_x, temp, N_square);

    recv_mpz(csp_sock, H);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(H);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }
    mpz_neg(temp, r_x);
    mpz_powm(S_3, e_y, temp, N_square);
    // 计算 e_x_y
    mpz_mul(temp, H, S_1);
    mpz_mod(temp, temp, N_square);
    mpz_mul(result, temp, S_2);
    mpz_mod(result, result, N_square);
    mpz_mul(result, result, S_3);
    mpz_mod(result, result, N_square);

    // 清理内存
    mpz_clears(flag, temp, e_r_x, e_r_y, X, Y, X1, Y1, h, H, S_1, S_2, S_3, r_x, r_y, NULL);

}

/**
 * POCF比较
 * @param e_u
 * @param e_x
 * @param e_y
 * @param cp
 * @param record_bytes_flag
 */
void CP_SCMP(mpz_t e_u, const mpz_t &e_x, const mpz_t &e_y, const Paillier_Third_Party &cp, bool record_bytes_flag)
{

    // 正常步骤
    mpz_t eone, temp, e_x_1, e_y_1, new_modulus, e_l, K, l, u_another;
    mpz_inits(eone, temp, e_x_1, e_y_1, new_modulus, e_l, K, l, u_another, NULL);

    // step1 cp
    // e_x_1 = 2x+1, e_y_1=2y
    Enc(eone, cp.public_key, one);
    mpz_powm(e_x_1, e_x, two, N_square);
    mpz_mul(e_x_1, e_x_1, eone);
    mpz_mod(e_x_1, e_x_1, N_square);
    mpz_powm(e_y_1, e_y, two, N_square);

    srand(time(NULL));
    int s = rand() % 2; // 生成0或1的随机数
    int N_bits = mpz_sizeinbase(N, 2);
    mpz_powm_ui(new_modulus, two, (N_bits / 4 - 1), N_square);
    mpz_t r;
    mpz_init(r);
    mpz_urandomm(r, randstate, new_modulus);
    if (mpz_cmp_si(r, 0) == 0)
    {
        mpz_add(r, r, one);
    }

    if (s == 1)
    {
        mpz_powm(temp, e_y_1, neg_one, N_square); // -y1
        mpz_mul(temp, e_x_1, temp); // x1-y1
        mpz_mod(temp, temp, N_square);
        mpz_powm(e_l, temp, r, N_square); // r(x1-y1)
    } else
    {
        mpz_powm(temp, e_x_1, neg_one, N_square); // -x1
        mpz_mul(temp, e_y_1, temp); // y1-x1
        mpz_mod(temp, temp, N_square);
        mpz_powm(e_l, temp, r, N_square); // r(y1-x1)
    }
    PDec(K, cp.private_key, e_l);
    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 2);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    // send K and e_l to CSP
    send_mpz(csp_sock, K);
    send_mpz(csp_sock, e_l);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(K);
        src_vec.push_back(e_l);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    // step3 CP
    recv_mpz(csp_sock, u_another);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(u_another);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }
    if (s == 1)
    {
        CR(u_another, N, N_square);
        mpz_set(e_u, u_another);
    } else if (s == 0)
    {
        mpz_set(temp, eone);
        mpz_powm(u_another, u_another, neg_one, N_square);
        mpz_mul(e_u, temp, u_another);
        mpz_mod(e_u, e_u, N_square);
    }
    mpz_clears(eone, temp, e_x_1, e_y_1, new_modulus, e_l, K, l, u_another, NULL);

}

/**
 * 符号位获取
 * @param result_e_s
 * @param result_x_another
 * @param e_x
 * @param cp
 * @param record_bytes_flag
 */
void CP_SSBA(mpz_t result_e_s, mpz_t result_x_another, const mpz_t &e_x, const Paillier_Third_Party &cp,
             bool record_bytes_flag)
{

    // 正常步骤
    mpz_t eone, new_modulus, r, temp, e_l, L, e_u, e_s, e_y, x_another, l;
    mpz_inits(eone, new_modulus, r, temp, e_l, L, e_u, e_s, e_y, x_another, l, NULL);

    Enc(eone, cp.public_key, one);

    // step 1 CP
    int s = rand() % 2; // 生成0或1的随机数
    int N_bits = mpz_sizeinbase(N, 2);
    mpz_powm_ui(new_modulus, two, N_bits / 4, N_square);
    mpz_urandomm(r, randstate, new_modulus);
    if (mpz_cmp_si(r, 0) == 0)
    {
        mpz_add(r, r, one);
    }
    if (s == 1)
    {
        mpz_powm(e_l, e_x, two, N_square);
        mpz_mul(e_l, e_l, eone);
        mpz_mod(e_l, e_l, N_square);
        mpz_powm(e_l, e_l, r, N_square);
    } else
    {
        mpz_powm(e_l, e_x, two, N_square);
        mpz_mul(e_l, e_l, eone);
        mpz_mod(e_l, e_l, N_square);
        mpz_t neg_r;
        mpz_init(neg_r);
        mpz_neg(neg_r, r);
        mpz_powm(e_l, e_l, neg_r, N_square);
    }
    PDec(L, cp.private_key, e_l);

    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 3);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);// 将这个 e_l 和 L 发给 csp
    send_mpz(csp_sock, e_l);
    send_mpz(csp_sock, L);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(e_l);
        src_vec.push_back(L);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    // 将这个 e_u 发给 cp
    recv_mpz(csp_sock, e_u);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(e_u);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }
    // step 3 cp
    if (s == 1)
    {
        CR(e_u, N, N_square);
        mpz_set(e_s, e_u);
    } else
    {
        mpz_powm(e_s, e_u, neg_one, N_square);
        mpz_mul(e_s, eone, e_s);
        mpz_mod(e_s, e_s, N_square);
        CR(e_s, N, N_square);
    }
    mpz_powm(e_y, e_s, two, N_square);
    mpz_powm(temp, eone, neg_one, N_square);
    mpz_mul(e_y, e_y, temp);
    mpz_mod(e_y, e_y, N_square);
    CP_SMUL(x_another, e_x, e_y, cp, record_bytes_flag);

    mpz_set(result_e_s, e_s);
    mpz_set(result_x_another, x_another);

    mpz_clears(eone, new_modulus, r, temp, e_l, L, e_u, e_s, e_y, x_another, l, NULL);

}


/**
 * 异或运算，限定范围 x, y ∈{0, 1}
 * @param result
 * @param e_x
 * @param e_y
 * @param cp
 * @param record_bytes_flag
 */
void SXOR(mpz_t result, const mpz_t &e_x, const mpz_t &e_y, const Paillier_Third_Party &cp, bool record_bytes_flag)
{
    mpz_t temp, e_1, f_1, f_2;
    mpz_inits(temp, e_1, f_1, f_2, NULL);

    // 计算 e_1
    Enc(e_1, cp.public_key, one);

    // 计算 f_1
    mpz_t inv_e_x;
    mpz_init(inv_e_x);
    mpz_powm(inv_e_x, e_x, neg_one, N_square);
    mpz_mul(inv_e_x, inv_e_x, e_1);
    mpz_mod(inv_e_x, inv_e_x, N_square);
    CP_SMUL(f_1, inv_e_x, e_y, cp, record_bytes_flag);
    mpz_clear(inv_e_x);


    // 计算 f_2
    mpz_t inv_e_y;
    mpz_init(inv_e_y);
    mpz_powm(inv_e_y, e_y, neg_one, N_square);
    mpz_mul(inv_e_y, inv_e_y, e_1);
    mpz_mod(inv_e_y, inv_e_y, N_square);
    CP_SMUL(f_2, e_x, inv_e_y, cp, record_bytes_flag);
    mpz_clear(inv_e_y);


    // 计算 f
    mpz_mul(result, f_1, f_2);
    mpz_mod(result, result, N_square);

    mpz_clears(temp, e_1, f_1, f_2, NULL);
}

/**
 * 相等判断
 * @param result
 * @param e_x
 * @param e_y
 * @param cp
 * @param record_bytes_flag
 */
void SEQ(mpz_t result, const mpz_t &e_x, const mpz_t &e_y, const Paillier_Third_Party &cp, bool record_bytes_flag)
{
    mpz_t u_1, u_2;
    mpz_inits(u_1, u_2, NULL);

    // 计算 u_1
    CP_SCMP(u_1, e_x, e_y, cp, record_bytes_flag);

    // 计算 u_2
    CP_SCMP(u_2, e_y, e_x, cp, record_bytes_flag);

    // 计算 f
    SXOR(result, u_1, u_2, cp, record_bytes_flag);

    mpz_clears(u_1, u_2, NULL);
}

/**
 * 在线版本部分解密调用
 * @param result
 * @param Y
 */
void CP_PDec(mpz_t result, mpz_t &Y)
{
    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);
    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 5);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    send_mpz(csp_sock, Y);
    recv_mpz(csp_sock, result);
    close(csp_sock);
}

void CP_Encrypted_LSB(mpz_t e_x_i, const mpz_t &T, const Paillier_Third_Party &cp, bool record_bytes_flag)
{

    mpz_t temp, eone, r, er, Y, Y1, Y2, y, alpha;
    mpz_inits(temp, eone, r, er, Y, Y1, Y2, y, alpha, NULL);

    Enc(eone, cp.public_key, one);

    // Bob:
    mpz_urandomm(r, randstate, N);
    Enc(er, cp.public_key, r);
    mpz_mul(Y, T, er);
    mpz_mod(Y, Y, N_square);
    PDec(Y1, cp.private_key, Y);
    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 6);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    // send Y1,Y to Alice
    send_mpz(csp_sock, Y1);
    send_mpz(csp_sock, Y);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(Y1);
        src_vec.push_back(Y);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    // Bob:
    recv_mpz(csp_sock, alpha);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(alpha);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }
    if (mpz_even_p(r))
    {
        mpz_set(e_x_i, alpha);
    } else
    {
        mpz_powm(Y1, alpha, neg_one, N_square);
        mpz_mul(e_x_i, eone, Y1);
        mpz_mod(e_x_i, e_x_i, N_square);
    }

    mpz_clears(temp, eone, r, er, Y, Y1, Y2, y, alpha, NULL);

}

int CP_SVR(const mpz_t &e_x, const vector<mpz_t *> &enc_bits, const Paillier_Third_Party &cp, bool record_bytes_flag)
{

    mpz_t temp, U, V, r, W, W1, W2, D_w, result, part_result1, part_result2;
    mpz_inits(temp, U, V, r, W, W1, W2, D_w, result, part_result1, part_result2, NULL);

    // Bob
    mpz_set(U, *enc_bits[0]);
    for (int i = 1; i < enc_bits.size(); ++i)
    {
        mpz_set_si(temp, i);
        mpz_powm(temp, two, temp, N_square);
        mpz_powm(temp, *enc_bits[i], temp, N_square);
        mpz_mul(U, U, temp);
        mpz_mod(U, U, N_square);
    }

    mpz_powm(V, e_x, neg_one, N_square);
    mpz_mul(V, V, U);
    mpz_mod(V, V, N_square);

    mpz_sub(temp, N, one);
    mpz_urandomm(r, randstate, temp);
    mpz_powm(W, V, r, N_square);
    PDec(W1, cp.private_key, W);
    int csp_sock = client_socket_initial(const_cast<char *>("127.0.0.1"), 8083);

    // 协议标志
    mpz_t tag;
    mpz_init(tag);
    mpz_set_si(tag, 7);
    send_mpz(csp_sock, tag);
    mpz_clear(tag);

    // send Y1,Y to Alice
    send_mpz(csp_sock, W1);
    send_mpz(csp_sock, W);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> src_vec;
        src_vec.push_back(W1);
        src_vec.push_back(W);
        total_bytes += save_to_file("CP_output.txt", src_vec);
    }

    // Bob
    recv_mpz(csp_sock, result);
    recv_mpz(csp_sock, part_result2);
    close(csp_sock);
    if (record_bytes_flag)
    {
        vector <mpz_ptr> dec_vec;
        dec_vec.push_back(result);
        dec_vec.push_back(part_result2);
        total_bytes += save_to_file("CP_output.txt", dec_vec);
    }
    PDec(part_result1, cp.private_key, result);
    TDec(D_w, part_result1, part_result2, N_square, N, N_half);
    if (mpz_cmp_si(D_w, 0) == 0)
    {
        return 0;
    } else
    {
        return 1;
    }
    mpz_clears(temp, U, V, r, W, W1, W2, D_w, result, part_result1, part_result2, NULL);

}
/**
 * 安全位分解
 * @param enc_bits
 * @param e_x
 * @param m
 * @param cp
 * @param record_bytes_flag
 */
void SBD(vector<mpz_t *> &enc_bits, const mpz_t &e_x, const int m, const Paillier_Third_Party &cp, bool record_bytes_flag)
{
    mpz_t temp, l, T, Z;
    mpz_inits(temp, l, T, Z, NULL);

    mpz_invert(l, two, N);

    while (true)
    {
        mpz_set(T, e_x);
        enc_bits.clear();
        for (int i = 0; i < m; ++i)
        {
            mpz_t e_x_i;
            mpz_init(e_x_i);
            CP_Encrypted_LSB(e_x_i, T, cp, record_bytes_flag);
            mpz_t *ptr = (mpz_t *) malloc(sizeof(mpz_t));
            mpz_init_set(*ptr, e_x_i);
            enc_bits.push_back(ptr);
            mpz_powm(Z, e_x_i, neg_one, N_square);
            mpz_mul(Z, Z, T);
            mpz_mod(Z, Z, N_square);
            mpz_powm(T, Z, l, N_square);
        }

        int gama = CP_SVR(e_x, enc_bits, cp, record_bytes_flag);
        if (gama == 1)
        {
            break;
        }
    }
    reverse(enc_bits.begin(), enc_bits.end());  // 高位在前面
}

/**
 * 安全除法
 * @param result_q_xing
 * @param result_r_xing
 * @param e_x
 * @param e_y
 * @param l
 * @param cp
 * @param record_bytes_flag
 */
void SDIV(mpz_t result_q_xing, mpz_t result_r_xing, const mpz_t &e_x, const mpz_t &e_y, const int l,
          const Paillier_Third_Party &cp, bool record_bytes_flag)
{
    mpz_t s_x, x_xing, s_y, y_xing, e_0, e_1, e_r, e_q, e_f, e_1_f, e_fx, e_K_1, e_y_another, e_x_another;
    // 初始化mpz_t变量
    mpz_inits(s_x, x_xing, s_y, y_xing, e_0, e_1, e_r, e_q, e_f, e_1_f, e_fx, e_K_1, e_y_another, e_x_another, NULL);

    // 计算e_0和e_1
    Enc(e_0, cp.public_key, zero);
    Enc(e_1, cp.public_key, one);

    // cp and csp
    SEQ(e_f, e_x, e_0, cp, record_bytes_flag);

    // cp
    mpz_powm(e_1_f, e_f, neg_one, N_square);
    mpz_mul(e_1_f, e_1, e_1_f);
    mpz_mod(e_1_f, e_1_f, N_square);

    // cp and csp
    CP_SMUL(e_fx, e_f, e_x, cp, record_bytes_flag);
    CP_SMUL(e_y_another, e_f, e_y, cp, record_bytes_flag);

    // cp
    mpz_mul(e_x_another, e_fx, e_1_f);
    mpz_mod(e_x_another, e_x_another, N_square);

    // cp and csp
    CP_SSBA(s_x, x_xing, e_x_another, cp, record_bytes_flag);
    CP_SSBA(s_y, y_xing, e_y_another, cp, record_bytes_flag);

    // 计算enc_q_bits
    vector < mpz_t * > enc_q_bits;
    SBD(enc_q_bits, y_xing, l, cp, record_bytes_flag);

    // 初始化enc_a_bits
    vector < mpz_t * > enc_a_bits;
    for (int i = 0; i < l; ++i)
    {
        mpz_t *ptr = (mpz_t *) malloc(sizeof(mpz_t));
        mpz_init_set(*ptr, e_0);
        enc_a_bits.push_back(ptr);
    }

    int exponent = l - 1;
    // 循环迭代
    for (int _ = 0; _ < l; ++_)
    {
        for (int i = 0; i < l - 1; i++)
        {
            mpz_set(*enc_a_bits[i], *enc_a_bits[i + 1]);
        }
        mpz_set(*enc_a_bits[l - 1], *enc_q_bits[0]);
        for (int i = 0; i < l - 1; i++)
        {
            mpz_set(*enc_q_bits[i], *enc_q_bits[i + 1]);
        }

        mpz_t e_A, e_Q, e_B;
        mpz_inits(e_A, e_Q, e_B, NULL);
        mpz_set(e_A, *enc_a_bits[l - 1]); // 最后一个e_a，也就是最低位的比特
        for (int i = 0; i < l - 1; ++i)
        {
            mpz_t temp_e_a;
            mpz_init(temp_e_a);
            mpz_powm_ui(temp_e_a, two, exponent - i, N_square);
            mpz_powm(temp_e_a, *enc_a_bits[i], temp_e_a, N_square);
            mpz_mul(e_A, temp_e_a, e_A);
            mpz_mod(e_A, e_A, N_square);
            mpz_clear(temp_e_a);
        }

        CP_SCMP(e_Q, e_A, x_xing, cp, record_bytes_flag);
        mpz_powm(*enc_q_bits[l - 1], e_Q, neg_one, N_square);
        mpz_mul(*enc_q_bits[l - 1], e_1, *enc_q_bits[l - 1]);
        mpz_mod(*enc_q_bits[l - 1], *enc_q_bits[l - 1], N_square);

        mpz_powm(e_B, x_xing, neg_one, N_square);
        CP_SMUL(e_B, e_B, *enc_q_bits[l - 1], cp, record_bytes_flag);

        mpz_mul(e_A, e_A, e_B);
        mpz_mod(e_A, e_A, N_square);
        SBD(enc_a_bits, e_A, l, cp, record_bytes_flag);
        mpz_clears(e_A, e_Q, e_B, NULL);
    }

    // 计算e_r和e_q
    mpz_set(e_r, *enc_a_bits[l - 1]);
    for (int i = 0; i < l - 1; ++i)
    {
        mpz_t temp_e_a;
        mpz_init(temp_e_a);
        mpz_powm_ui(temp_e_a, two, exponent - i, N_square);
        mpz_powm(temp_e_a, *enc_a_bits[i], temp_e_a, N_square);
        mpz_mul(e_r, temp_e_a, e_r);
        mpz_mod(e_r, e_r, N_square);
        mpz_clear(temp_e_a);
    }
    mpz_set(e_q, *enc_q_bits[l - 1]);
    for (int i = 0; i < l - 1; ++i)
    {
        mpz_t temp_e_q;
        mpz_init(temp_e_q);
        mpz_powm_ui(temp_e_q, two, exponent - i, N_square);
        mpz_powm(temp_e_q, *enc_q_bits[i], temp_e_q, N_square);
        mpz_mul(e_q, temp_e_q, e_q);
        mpz_mod(e_q, e_q, N_square);
        mpz_clear(temp_e_q);
    }

    // 计算e_K_1、r_xing和q_xing
    CP_SMUL(e_K_1, s_x, s_y, cp, record_bytes_flag);
    CP_SMUL(result_r_xing, e_r, s_y, cp, record_bytes_flag);
    CP_SMUL(result_q_xing, e_q, e_K_1, cp, record_bytes_flag);

    enc_a_bits.clear();
    enc_q_bits.clear();

    // 清理内存
    mpz_clears(s_x, x_xing, s_y, y_xing, e_0, e_1, e_r, e_q, e_f, e_1_f, e_fx, e_K_1, e_y_another, e_x_another, NULL);
}


/**
 * 乘法正确性检验
 * @param i
 * @param result_cxy
 * @param x
 * @param y
 * @param private_key
 */
void check_smul(int i, mpz_t result_cxy, mpz_t x, mpz_t y, PrivateKey private_key)
{
    mpz_t z, xy;
    mpz_inits(z, xy, NULL);
    Dec(z, private_key, result_cxy);
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
 * @param private_key
 */
void check_scmp(int i, mpz_t result_cmp_x_y, mpz_t x, mpz_t y, PrivateKey private_key)
{
    mpz_t z, answer;
    mpz_inits(z, answer, NULL);
    if (mpz_cmp(x, y) >= 0)
    {
        mpz_set_si(answer, 0);
    } else mpz_set_si(answer, 1);

    Dec(z, private_key, result_cmp_x_y);
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
 * 符号位获取正确性经验，结果与SOCI中SSBA的表示相反
 * @param i
 * @param result_cs
 * @param result_cu
 * @param x
 * @param pai
 */
void check_ssba(int i, mpz_t result_cs, mpz_t result_cu, mpz_t x, PrivateKey private_key)
{
    mpz_t s, u, answer_s, answer_u;
    mpz_inits(s, u, answer_s, answer_u, NULL);

    Dec(s, private_key, result_cs);
    Dec(u, private_key, result_cu);

    if (mpz_cmp_si(x, 0) >= 0)
    {
        mpz_set_si(answer_s, 1);
        mpz_set(answer_u, x);
    } else
    {
        mpz_set_si(answer_s, 0);
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
 * @param private_key
 */
void check_sdiv(int i, mpz_t result_cq, mpz_t result_ce, mpz_t x, mpz_t y, PrivateKey private_key)
{
    mpz_t q, e, right_q, right_e;
    mpz_inits(q, e, right_q, right_e, NULL);
    Dec(q, private_key, result_cq);
    Dec(e, private_key, result_ce);

    mpz_div(right_q, y, x);
    mpz_mod(right_e, y, x);
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
 * 相等判断正确性检验
 * @param i
 * @param result_eq_x_y
 * @param x
 * @param y
 * @param private_key
 */
void check_seq(int i, mpz_t result_eq_x_y, mpz_t x, mpz_t y, PrivateKey private_key)
{
    mpz_t z, answer;
    mpz_inits(z, answer, NULL);
    if (mpz_cmp(x, y) == 0)
    {
        mpz_set_si(answer, 0);
    } else mpz_set_si(answer, 1);

    Dec(z, private_key, result_eq_x_y);
    if (mpz_cmp(z, answer) != 0)
    {
        gmp_printf("dec answer x==y? = %Zd\n", z);
        gmp_printf("right answer %Zd\n", answer);
        printf("SEQ error!\n");
        printf("i = %d is error \n", i);
        exit(-1);
    }
    mpz_clears(z, answer, NULL);
}


#endif