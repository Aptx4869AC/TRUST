/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include "include/pcpd.h"
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
    Paillier_Third_Party csp;
};

/**
 * 线程处理函数
 * @param arg
 * @return
 */
void *handle_client(void *arg)
{
    // Cast argument back to the correct type
    ThreadArgs *args = static_cast<ThreadArgs *>(arg);
    int cp_socket = args->cp_socket;
    Paillier_Third_Party csp = args->csp;

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
        CSP_SSBA(csp, cp_socket);

    } else if (mpz_cmp_si(args->flag, 4) == 0)
    {

    } else if (mpz_cmp_si(args->flag, 5) == 0)
    {
        CSP_PDec(csp, cp_socket);

    } else if (mpz_cmp_si(args->flag, 6) == 0)
    {
        CSP_Encrypted_LSB(csp, cp_socket);

    } else if (mpz_cmp_si(args->flag, 7) == 0)
    {
        CSP_SVR(csp, cp_socket);

    }
//    printf("----------------------------------------------------------\n");

    close(cp_socket);
    pthread_exit(NULL);
}

int main()
{
    double start_time, end_time;
    double average_time_SMUL = 0;
    double average_time_SCMP = 0;
    double average_time_SSBA = 0;
    double average_time_SDIV = 0;

    int server_sock, data_owner_sock;
    server_socket_initial(server_sock, data_owner_sock, 8082);

    mpz_inits(zero, one, two, neg_one, N, N_square, N_half, NULL);
    mpz_set_si(zero, 0);
    mpz_set_si(one, 1);
    mpz_set_si(two, 2);
    mpz_set_si(neg_one, -1);

    // 接收数据
    mpz_t lambda, lambda_reverse;
    mpz_inits(lambda, lambda_reverse, N, NULL);
    recv_mpz(data_owner_sock, lambda);
    recv_mpz(data_owner_sock, lambda_reverse);
    recv_mpz(data_owner_sock, N);
    close(data_owner_sock);
    close(server_sock);

    gmp_printf("------ sk2 = %Zd\n", lambda);
    gmp_printf("------ n = %Zd\n", N);

    PublicKey public_key(N, sigma);
    PrivateKey private_key_2(lambda, lambda_reverse, N);
    Paillier_Third_Party csp(public_key, private_key_2);

    mpz_set(N, csp.public_key.N);
    mpz_set(N_square, csp.public_key.N_square);
    mpz_set(N_half, csp.public_key.N_half);
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

//        int new_send_buffer_size = 1024; // 1 GB
//        int new_recv_buffer_size = 1024; // 1 GB
//        if (!setSocketBufferSize(server_socket, new_send_buffer_size, new_recv_buffer_size)) {
//            std::cerr << "Failed to set socket buffer sizes\n";
//            return 1;
//        }
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