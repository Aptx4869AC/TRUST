/******************************************************************************
 * Author: Aptx4869AC 
 * Created: 2024-12-01 
 * GitHub: https://github.com/Aptx4869AC  
 *****************************************************************************/

#include "include/fastPai.h"
#include "include/network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

int main() {

    #ifdef _WIN32 
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2,2), &wsaData) != 0) {
        std::cerr << "WSAStartup failed: " << WSAGetLastError() << std::endl;
        return EXIT_FAILURE;
    }
    #endif


    double start_time, end_time;
    double average_time_keygen = 0;
    double average_time_enc = 0;
    double average_time_dec = 0;
    double average_time_PDec_TDec = 0;
    double average_time_HE_Addition = 0;
    double average_time_HE_ScalarMul = 0;
    double average_time_HE_Subtraction = 0;

    int k = 112; //BIT SECURITY, 64、80、104、112、124、128
    KeyGen keyGen(k, sigma);
    Paillier pai(keyGen.pk, keyGen.sk, sigma);

    PaillierThirdParty *psk = keyGen.pai_third_parties;

    // 初始化CP
    int sock_to_cp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8081);
    mpz_t L_k;
    mpz_init(L_k);
    mpz_set_si(L_k, keyGen.pk.L_k);

    // 发给CP <sk1, pk , ex, ey>
    gmp_printf("Send (PK,SK) to CP, sk1 = %Zd\n", psk[0].partial_key);
    send_mpz(sock_to_cp, psk[0].partial_key);
    send_mpz(sock_to_cp, keyGen.pk.N); // N
    send_mpz(sock_to_cp, keyGen.pk.h); // h
    send_mpz(sock_to_cp, keyGen.pk.h_N); // h_N
    send_mpz(sock_to_cp, keyGen.pk.L); // L
    send_mpz(sock_to_cp, L_k); //L_k
    send_mpz(sock_to_cp, keyGen.sk.alpha); //用于验证

    // 初始化CSP
    int sock_to_csp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8082);

    // 发给CSP <sk2, pk>
    gmp_printf("Send (PK,SK) to CSP, sk2 = %Zd\n", psk[1].partial_key);
    send_mpz(sock_to_csp, psk[1].partial_key);
    send_mpz(sock_to_csp, keyGen.pk.N); // N
    send_mpz(sock_to_csp, keyGen.pk.h); // h
    send_mpz(sock_to_csp, keyGen.pk.h_N); // h_N
    send_mpz(sock_to_csp, keyGen.pk.L); // L
    send_mpz(sock_to_csp, L_k); //L_k

    gmp_printf("------ N = %Zd\n", keyGen.pk.N);
    gmp_printf("------ h = %Zd\n", keyGen.pk.h);
    gmp_printf("------ h_N = %Zd\n", keyGen.pk.h_N);
    gmp_printf("------ L = %Zd\n", keyGen.pk.L);
    gmp_printf("------ L_k = %Zd\n", L_k);

    // 收到结束标志
    close(sock_to_cp); // printf("关闭CP的连接\n");
    close(sock_to_csp); // printf("关闭CSP的连接\n");

    #ifdef _WIN32 
    WSACleanup();
    #endif

    return 0;
}