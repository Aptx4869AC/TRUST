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

int main()
{

#ifdef _WIN32
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2, 2), &wsaData) != 0)
    {
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

    int k = 80; // BIT SECURITY, 64、80、104、112、124、128
    KeyGen keyGen(k, sigma);
    Paillier pai(keyGen.pk, keyGen.sk, sigma);

    PaillierThirdParty *psk = keyGen.pai_third_parties;

    // 初始化CP
    int sock_to_cp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8081);
    mpz_t L_k;
    mpz_init(L_k);
    mpz_set_si(L_k, keyGen.pk.L_k);

    // 发给CP <sk1, pk , ex, ey>
    gmp_printf("Send (PK, SK_1) to CP, sk1 = %Zd\n", psk[0].partial_key);
    send_mpz(sock_to_cp, psk[0].partial_key);
    send_mpz(sock_to_cp, keyGen.pk.N);     // N
    send_mpz(sock_to_cp, keyGen.pk.h);     // h
    send_mpz(sock_to_cp, keyGen.pk.h_N);   // h_N
    send_mpz(sock_to_cp, keyGen.pk.L);     // L
    send_mpz(sock_to_cp, L_k);             // L_k
    send_mpz(sock_to_cp, keyGen.sk.alpha); // 用于验证

    // 初始化CSP
    int sock_to_csp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8082);

    // 发给CSP <sk2, pk>
    gmp_printf("Send (PK, SK_2) to CSP, sk2 = %Zd\n", psk[1].partial_key);
    send_mpz(sock_to_csp, psk[1].partial_key);
    send_mpz(sock_to_csp, keyGen.pk.N);   // N
    send_mpz(sock_to_csp, keyGen.pk.h);   // h
    send_mpz(sock_to_csp, keyGen.pk.h_N); // h_N
    send_mpz(sock_to_csp, keyGen.pk.L);   // L
    send_mpz(sock_to_csp, L_k);           // L_k

    gmp_printf("------ N = %Zd\n", keyGen.pk.N);
    gmp_printf("------ h = %Zd\n", keyGen.pk.h);
    gmp_printf("------ h_N = %Zd\n", keyGen.pk.h_N);
    gmp_printf("------ L = %Zd\n", keyGen.pk.L);
    gmp_printf("------ L_k = %Zd\n", L_k);

    close(sock_to_cp);  // printf("关闭CP的初始化连接\n");
    close(sock_to_csp); // printf("关闭CSP的初始化连接\n");

    printf("----------------------------------------------------------\n");
    ifstream file("../Data/training_dataADD.csv"); // 打开CSV文件
    string line, value;

    if (!file.is_open())
    {
        cerr << "Failed to open the file training_dataADD.csv" << '\n';
        return 1;
    }
    // 跳过表头
    getline(file, line);

    // 外包计算任务，重新与CP连接
    sock_to_cp = client_socket_initial(const_cast<char *>("127.0.0.1"), 8084);

    int epoch = 1;

    // Trust 协议
    while (getline(file, line)) // 逐行读取CSV文件的内容
    {
        stringstream ss(line);
        getline(ss, value, ',');
        int value1 = stoi(value); // 转换为整数
        getline(ss, value, ',');
        int value2 = stoi(value); // 转换为整数

        mpz_t m_1, m_2, m_3, em_1, em_2, em_3, eres, res;
        mpz_inits(m_1, m_2, m_3, em_1, em_2, em_3, eres, res, NULL);
        mpz_set_si(m_1, value1);
        mpz_set_si(m_2, value2);
        mpz_rrandomb(m_3, randstate, 8);
        gmp_printf("------ m_1 = %Zd\n", m_1);
        gmp_printf("------ m_2 = %Zd\n", m_2);
        gmp_printf("------ m_3 = %Zd\n", m_3);
        if (mpz_cmp_si(m_1, 0) < 0)
        {
            mpz_add(m_1, m_1, pai.pk.N);
        }
        if (mpz_cmp_si(m_2, 0) < 0)
        {
            mpz_add(m_2, m_2, pai.pk.N);
        }
        pai.encrypt(em_1, m_1);
        pai.encrypt(em_2, m_2);
        pai.encrypt(em_3, m_3);
        gmp_printf("------ em_1 = %Zd\n", em_1);
        gmp_printf("------ em_2 = %Zd\n", em_2);
        gmp_printf("------ em_3 = %Zd\n", em_3);

        printf("Select an operation to execute:\n");
        printf("--- 1. FMUL \n");
        printf("--- 2. FCMP \n");
        printf("--- 3. FEQL \n");
        printf("--- 4. FABS \n");
        printf("--- 5. FTRN \n");

        int operation_type;
        cin >> operation_type;
        mpz_t mpz_type;
        mpz_init(mpz_type);
        mpz_set_si(mpz_type, operation_type);
        switch (operation_type)
        {
        case 1: // FMUL
            printf("Send (em_1, em_2) to CP, operation_type = FMUL\n");
            send_mpz(sock_to_cp, mpz_type);
            send_mpz(sock_to_cp, em_1);
            send_mpz(sock_to_cp, em_2);

            recv_mpz(sock_to_cp, eres);             // 接数据
            check_smul(epoch, eres, m_1, m_2, pai); // 验证正确性

            break;
        case 2: // FCMP
            printf("Send (em_1, em_2) to CP, operation_type = FCMP\n");
            send_mpz(sock_to_cp, mpz_type);
            send_mpz(sock_to_cp, em_1);
            send_mpz(sock_to_cp, em_2);

            recv_mpz(sock_to_cp, eres);             // 接数据
            check_scmp(epoch, eres, m_1, m_2, pai); // 验证正确性

            break;
        case 3: // FEQL
            printf("Send (em_1, em_2) to CP, operation_type = FEQL\n");
            send_mpz(sock_to_cp, mpz_type);
            send_mpz(sock_to_cp, em_1);
            send_mpz(sock_to_cp, em_2);

            recv_mpz(sock_to_cp, eres);             // 接数据
            check_feql(epoch, eres, m_1, m_2, pai); // 验证正确性

            break;
        case 4: // FABS
            printf("Send (em_1, em_2) to CP, operation_type = FABS\n");
            send_mpz(sock_to_cp, mpz_type);
            send_mpz(sock_to_cp, em_1);
            send_mpz(sock_to_cp, em_2);

            recv_mpz(sock_to_cp, eres);        // 接数据
            check_fabs(epoch, eres, m_1, pai); // 验证正确性

            break;
        case 5: // FTRN
            printf("Send (em_1, em_2, em_3) to CP, operation_type = FTRN\n");
            mpz_set_si(m_1, 1);
            pai.encrypt(em_1, m_1);
            send_mpz(sock_to_cp, mpz_type);
            send_mpz(sock_to_cp, em_1);
            send_mpz(sock_to_cp, em_2);
            send_mpz(sock_to_cp, em_3);

            recv_mpz(sock_to_cp, eres);                  // 接数据
            check_ftrn(epoch, eres, m_1, m_2, m_3, pai); // 验证正确性

            break;
        default:
            printf("[WARNING] Invalid operation selected\n");
            break;
        }
        epoch++;

        int input;
        cout << "Enter 9 to continue, 0 to exit: ";
        cin >> input;
        if (input == 0)
        {
            break;
        }
    }

#ifdef _WIN32
    WSACleanup();
#endif

    return 0;
}