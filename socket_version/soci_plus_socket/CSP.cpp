/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-01
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include "include/fastPai.h"
#include "include/csp_protocol.h"
#include "include/network.h"

using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;

#define PORT 8083
#define MAX_PENDING_CONNECTIONS 20


struct ThreadArgs
{
    int cp_socket;
    mpz_t flag;
    PaillierThdDec csp;
};


/**
 * 线程执行函数
 * @param arg
 * @return
 */
void *handle_client(void *arg)
{
    // Cast argument back to the correct type
    ThreadArgs *args = static_cast<ThreadArgs *>(arg);
    int cp_socket = args->cp_socket;
    PaillierThdDec csp = args->csp;

    if (mpz_cmp_si(args->flag, 1) == 0)
    {
        CSP_SMUL(csp, cp_socket);
//        printf("Processing Protocol SMUL\n");
    } else if (mpz_cmp_si(args->flag, 2) == 0)
    {
        CSP_SCMP(csp, cp_socket);
//        printf("Processing Protocol SCMP\n");
    } else if (mpz_cmp_si(args->flag, 3) == 0)
    {
        CSP_PDec(csp, cp_socket);
//        printf("Processing Protocol PDec\n");
    } else if (mpz_cmp_si(args->flag, 4) == 0)
    {
        CSP_TEE_SMUL(csp, cp_socket);

    } else if (mpz_cmp_si(args->flag, 5) == 0)
    {
        CSP_TEE_SCMP_version1(csp, cp_socket);

    } else if (mpz_cmp_si(args->flag, 6) == 0)
    {
        CSP_TEE_SCMP_version2(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 7) == 0)
    {
        CSP_TEE_SABS(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 8) == 0)
    {
        CSP_Trust_FEQL(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 9) == 0)
    {
        CSP_Trust_FABS(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 10) == 0)
    {
        CSP_Trust_FTRN_version1(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 11) == 0)
    {
        CSP_Trust_FTRN_version2(csp, cp_socket);
    }

//    printf("----------------------------------------------------------\n");

    close(cp_socket);
    pthread_exit(NULL);
}

int main()
{

    int server_sock, data_owner_sock;
    server_socket_initial(server_sock, data_owner_sock, 8082);

    // 接收数据
    mpz_t sk2, pk_N, pk_h, pk_h_N, pk_L, pk_L_k;
    mpz_inits(sk2, pk_N, pk_h, pk_h_N, pk_L, pk_L_k, NULL);
    recv_mpz(data_owner_sock, sk2);
    recv_mpz(data_owner_sock, pk_N); // N
    recv_mpz(data_owner_sock, pk_h); // h
    recv_mpz(data_owner_sock, pk_h_N); // h_N
    recv_mpz(data_owner_sock, pk_L); // L
    recv_mpz(data_owner_sock, pk_L_k); //L_k

    close(data_owner_sock);
    close(server_sock);

    gmp_printf("------ sk2 = %Zd\n", sk2);
    gmp_printf("------ N = %Zd\n", pk_N);
    gmp_printf("------ h = %Zd\n", pk_h);
    gmp_printf("------ h_N = %Zd\n", pk_h_N);
    gmp_printf("------ L = %Zd\n", pk_L);
    gmp_printf("------ L_k = %Zd\n", pk_L_k);

    // 配置自身
    setrandom(&randstate);
    PublicKey pk(pk_N, pk_h, pk_h_N, pk_L, mpz_get_si(pk_L_k));
    PaillierThdDec csp = PaillierThdDec(pk_N, sk2, pk, sigma);
    printf("----------------------------------------------------------\n");

    int server_socket, cp_socket;
    struct sockaddr_in server_addr, client_addr;
    socklen_t addr_len = sizeof(client_addr);
    pthread_t thread;

    // 创建服务器 socket
    if ((server_socket = socket(AF_INET, SOCK_STREAM, 0)) == 0)
    {
        perror("socket creation failed");
        exit(EXIT_FAILURE);
    }

    // 设置服务器地址结构
    server_addr.sin_family = AF_INET;
    server_addr.sin_addr.s_addr = INADDR_ANY;
    server_addr.sin_port = htons(PORT);

    // 绑定服务器 socket 到指定地址和端口
    if (bind(server_socket, (struct sockaddr *) &server_addr, sizeof(server_addr)) < 0)
    {
        perror("bind failed");
        exit(EXIT_FAILURE);
    }

    // 开始监听连接
    if (listen(server_socket, MAX_PENDING_CONNECTIONS) < 0)
    {
        perror("listen failed");
        exit(EXIT_FAILURE);
    }

    printf("Server listening on port %d...\n", PORT);

    while (1)
    {
        // 等待客户端连接
        if ((cp_socket = accept(server_socket, (struct sockaddr *) &client_addr, &addr_len)) < 0)
        {
            perror("accept failed");
            continue;
        }

//        printf("Connection accepted from %s:%d\n", inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));

        // 协议标志
        mpz_t flag;
        mpz_init(flag);
        recv_mpz(cp_socket, flag);

        ThreadArgs *args = new ThreadArgs;
        args->cp_socket = cp_socket;
        args->csp = csp;
        mpz_init(args->flag);
        mpz_set(args->flag, flag);

        // 创建一个新线程来处理客户端请求
        if (pthread_create(&thread, NULL, handle_client, static_cast<void *>(args)) != 0)
        {
            perror("pthread_create failed");
            close(cp_socket);
            continue;
        }

        // 分离线程，不等待线程结束，避免僵尸线程
        pthread_detach(thread);
    }

    close(server_socket);

    return 0;
}