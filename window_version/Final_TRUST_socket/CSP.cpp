/******************************************************************************
 * Author: Aptx4869AC 
 * Created: 2024-12-01 
 * GitHub: https://github.com/Aptx4869AC  
 * Last Modified: 2025-05-13 (Fixed pthread compatibility for MinGW)
 *****************************************************************************/
 
#include "include/fastPai.h"
#include "include/csp_protocol.h"
#include "include/network.h"
 
// 跨平台线程头文件 
#ifdef _WIN32 
#include <windows.h>
#include <process.h>  // Windows线程API 
#define pthread_t HANDLE 
#define pthread_exit(status) { _endthreadex(0); return NULL; }
#else 
#include <pthread.h>  // POSIX线程API 
#endif 
 
using namespace std;
using namespace PHESPACE;
using namespace NETWORKSPACE;
 
#define PORT 8083 
#define MAX_PENDING_CONNECTIONS 20 
 
struct ThreadArgs {
    int cp_socket;
    mpz_t flag;
    PaillierThdDec csp;
};
 
/**
 * 线程执行函数（跨平台兼容版本）
 */
void* handle_client(void* arg) {
    ThreadArgs* args = static_cast<ThreadArgs*>(arg);
    int cp_socket = args->cp_socket;
    PaillierThdDec csp = args->csp;
 
    // 协议选择逻辑（保持不变）
    if (mpz_cmp_si(args->flag, 1) == 0) {
        CSP_SMUL(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 2) == 0) {
        CSP_SCMP(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 3) == 0) {
        CSP_PDec(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 4) == 0) {
        CSP_TEE_SMUL(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 6) == 0) {
        CSP_TEE_SCMP_version2(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 7) == 0) {
        CSP_TEE_SABS(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 8) == 0) {
        CSP_Trust_FEQL(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 9) == 0) {
        CSP_Trust_FABS(csp, cp_socket);
    } else if (mpz_cmp_si(args->flag, 10) == 0) {
        CSP_Trust_FTRN_version1(csp, cp_socket);
    }
 
    // 资源清理 
    close(cp_socket);
    delete args;  // 防止内存泄漏 
 
    // 跨平台线程退出 
#ifdef _WIN32 
    return 0;
#else 
    pthread_exit(NULL);
#endif 
}
 

int main() {
    #ifdef _WIN32 
    WSADATA wsaData;
    if (WSAStartup(MAKEWORD(2,2), &wsaData) != 0) {
        std::cerr << "WSAStartup failed: " << WSAGetLastError() << std::endl;
        return EXIT_FAILURE;
    }
    #endif


    gmp_printf("CSP run success\n");

    // 初始化服务器socket（端口8082）
    int server_sock, data_owner_sock;
    server_socket_initial(server_sock, data_owner_sock, 8082);
 
    // 接收加密参数 
    mpz_t sk2, pk_N, pk_h, pk_h_N, pk_L, pk_L_k;
    mpz_inits(sk2, pk_N, pk_h, pk_h_N, pk_L, pk_L_k, NULL);
    recv_mpz(data_owner_sock, sk2);
    recv_mpz(data_owner_sock, pk_N);
    recv_mpz(data_owner_sock, pk_h);
    recv_mpz(data_owner_sock, pk_h_N);
    recv_mpz(data_owner_sock, pk_L);
    recv_mpz(data_owner_sock, pk_L_k);
 
    // 打印接收到的参数 
    gmp_printf("------ sk2 = %Zd\n", sk2);
    gmp_printf("------ N = %Zd\n", pk_N);
    gmp_printf("------ h = %Zd\n", pk_h);
    gmp_printf("------ h_N = %Zd\n", pk_h_N);
    gmp_printf("------ L = %Zd\n", pk_L);
    gmp_printf("------ L_k = %Zd\n", pk_L_k);
 
    // 初始化CSP 
    setrandom(&randstate);
    PublicKey pk(pk_N, pk_h, pk_h_N, pk_L, mpz_get_si(pk_L_k));
    PaillierThdDec csp = PaillierThdDec(pk_N, sk2, pk, sigma);
    printf("----------------------------------------------------------\n");
 
    // 主服务socket（端口8083）
    int server_socket, cp_socket;
    struct sockaddr_in server_addr, client_addr;
    socklen_t addr_len = sizeof(client_addr);
 
    // 创建主服务socket 
    if ((server_socket = socket(AF_INET, SOCK_STREAM, 0)) == 0) {
        perror("socket creation failed");
        exit(EXIT_FAILURE);
    }
 
    // 配置服务器地址 
    server_addr.sin_family  = AF_INET;
    server_addr.sin_addr.s_addr  = INADDR_ANY;
    server_addr.sin_port  = htons(PORT);
 
    // 绑定并监听 
    if (bind(server_socket, (struct sockaddr*)&server_addr, sizeof(server_addr)) < 0) {
        perror("bind failed");
        exit(EXIT_FAILURE);
    }
    if (listen(server_socket, MAX_PENDING_CONNECTIONS) < 0) {
        perror("listen failed");
        exit(EXIT_FAILURE);
    }
    printf("Server listening on port %d...\n", PORT);
 
    // 主服务循环 
    while (1) {
        // 接受客户端连接 
        if ((cp_socket = accept(server_socket, (struct sockaddr*)&client_addr, &addr_len)) < 0) {
            perror("accept failed");
            continue;
        }
 
        // 接收协议标志 
        mpz_t flag;
        mpz_init(flag);
        recv_mpz(cp_socket, flag);
 
        // 创建线程参数 
        ThreadArgs* args = new ThreadArgs;
        args->cp_socket = cp_socket;
        args->csp = csp;
        mpz_init(args->flag);
        mpz_set(args->flag, flag);
 
        // 跨平台线程创建 
#ifdef _WIN32 
        HANDLE thread = (HANDLE)_beginthreadex(NULL, 0, 
            (unsigned(__stdcall*)(void*))handle_client, args, 0, NULL);
        if (thread == NULL) {
#else 
        pthread_t thread;
        if (pthread_create(&thread, NULL, handle_client, args) != 0) {
#endif 
            perror("Thread creation failed");
            close(cp_socket);
            delete args;
            continue;
        }
 
        // 线程分离（Windows自动分离）
#ifndef _WIN32 
        pthread_detach(thread);
#endif 
    }
 
    close(server_socket);

    #ifdef _WIN32 
    WSACleanup();
    #endif

    return 0;
}