/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/


#include "include/paillier.h"
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
    mpz_t sk2;
    mpz_t n;
    mpz_t nsquare;
};

// 线程执行函数
void *handle_client(void *arg)
{
    // Cast argument back to the correct type
    double start_time, end_time;
    start_time = omp_get_wtime();
    ThreadArgs *args = static_cast<ThreadArgs *>(arg);
    int cp_socket = args->cp_socket;
    mpz_t sk2, n, nsquare;
    mpz_inits(sk2, n, nsquare, NULL);
    mpz_set(sk2, args->sk2);
    mpz_set(n, args->n);
    mpz_set(nsquare, args->nsquare);

    setrandom(&randstate);
    PaillierKey pubkey(n);
    PaillierThdPrivateKey share_part2(sk2, n, nsquare);
    PaillierThdDec csp(share_part2, pubkey);


    if (mpz_cmp_si(args->flag, 1) == 0)
    {
        CSP_SMUL(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 2) == 0)
    {
        CSP_SCMP(csp, cp_socket);
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
    end_time = omp_get_wtime();

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

    // 接收数据
    mpz_t sk2, n, nsquare;
    mpz_inits(sk2, n, nsquare, NULL);
    recv_mpz(data_owner_sock, sk2);
    recv_mpz(data_owner_sock, n);
    server_socket_close(data_owner_sock, server_sock);


    mpz_mul(nsquare, n, n);
    gmp_printf("------ sk2 = %Zd\n", sk2);
    gmp_printf("------ n = %Zd\n", n);
    gmp_printf("------ nsquare = %Zd\n", nsquare);

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
    mpz_set_si(neg_one, -1);
    mpz_set_si(one, 1);
    mpz_set_si(zero, 0);

    ThreadArgs *args = new ThreadArgs;
    mpz_inits(args->flag, args->sk2, args->n, args->nsquare, NULL);
    mpz_set(args->sk2, sk2);
    mpz_set(args->n, n);
    mpz_set(args->nsquare, nsquare);

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
        args->cp_socket = cp_socket;
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