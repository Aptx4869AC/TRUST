/******************************************************************************
 * Author: Aptx4869AC 
 * Created: 2024-12-01 
 * GitHub: https://github.com/Aptx4869AC  
 *****************************************************************************/

#ifndef NETWORK_H
#define NETWORK_H
#pragma once

#include "fastPai.h"

using namespace std;

namespace NETWORKSPACE {
    void send_mpz(int sock, const mpz_t num) {
        size_t count = (mpz_sizeinbase(num, 2) + 7) / 8;
        unsigned char *buffer = new unsigned char[count];
        mpz_export(buffer, &count, 1, 1, 0, 0, num);
        // 发送数据长度
        // send(sock, &count, sizeof(size_t), 0);
        send(sock, reinterpret_cast<const char*>(&count), sizeof(size_t), 0);  // 强制类型转换 
        // 发送数据内容
        // send(sock, buffer, count, 0);
        send(sock, reinterpret_cast<const char*>(buffer), count, 0);  // 强制类型转换 
        delete[] buffer;
    }

    void recv_mpz(int sock, mpz_t num) {
        mpz_init(num);
        size_t count;
        recv(sock, reinterpret_cast<char*>(&count), sizeof(size_t), 0);
        // recv(sock, &count, sizeof(size_t), 0);
        unsigned char *buffer = new unsigned char[count];
        recv(sock, reinterpret_cast<char*>(buffer), count, 0);
        //recv(sock, buffer, count, 0);
        mpz_import(num, count, 1, 1, 0, 0, buffer);
        delete[] buffer;
    }


    int setServerSock(int port) {
        // 创建服务器 socket
        int server_sock = socket(AF_INET, SOCK_STREAM, 0);
        if (server_sock < 0) {
            cerr << "Socket creation error" << '\n';
            return -1;
        }

        // 配置服务器地址结构
        sockaddr_in server_addr;
        server_addr.sin_family = AF_INET;
        server_addr.sin_port = htons(port); // 指定监听的端口号
        server_addr.sin_addr.s_addr = inet_addr("127.0.0.1"); // 监听所有本地 IP 地址

        // 绑定 socket 到指定的端口和地址
        if (bind(server_sock, (struct sockaddr *) &server_addr, sizeof(server_addr)) < 0) {
            cerr << "Bind failed" << '\n';
            close(server_sock);
            return -1;
        }

        // 开始监听端口，等待客户端（DO）连接
        if (listen(server_sock, 3) < 0) {
            cerr << "Listen failed" << '\n';
            close(server_sock);
            return -1;
        }
        printf("Server is listening on port %d...\n", port);
        return server_sock;

    }


    int connectTo(int port) {
        // 连接CSP
        int sock = socket(AF_INET, SOCK_STREAM, 0);
        if (sock < 0) {
            cerr << "Socket creation error" << '\n';
            return -1;
        }

        // 配置CSP地址架构
        sockaddr_in addr;
        addr.sin_family = AF_INET;
        addr.sin_port = htons(port); // 待连接端口号
        addr.sin_addr.s_addr = inet_addr("127.0.0.1"); // 目的IP地址

        // 连接到CSP
        if (connect(sock, (struct sockaddr *) &addr, sizeof(addr)) < 0) {
            cerr << "Connection failed" << '\n';
            close(sock);
            return -1;
        }

        printf("Connection ok\n");
        return sock; // 返回连接成功的socket描述符
    }

    int acceptFrom(int server_sock) {
        // 接受客户端连接
        sockaddr_in addr;
        socklen_t addr_len = sizeof(addr);
        int sock = accept(server_sock, (struct sockaddr *) &addr, &addr_len);
        if (sock < 0) {
            cerr << "Accept failed" << '\n';
            close(server_sock);
            return -1;
        }
        cout << "Connection accepted" << '\n';
        return sock;
    }

    /**
      * initialize server socket
      *
      * @param server_socket server socket
      * @param client_socket client socket
      * @param host_port host port
      */
    void server_socket_initial(int &server_socket, int &client_socket, int host_port = 7766) {
        server_socket = socket(AF_INET, SOCK_STREAM, 0);
        if (server_socket == -1) {
            cout << "Error: Failed to create socket." << endl;
            exit(EXIT_FAILURE);
        }

        struct sockaddr_in server_address{};
        server_address.sin_family = AF_INET;
        server_address.sin_port = htons(host_port);
        server_address.sin_addr.s_addr = INADDR_ANY;

        int opt = 1;
        // Set socket options to allow address reuse
        setsockopt(server_socket, SOL_SOCKET, SO_REUSEADDR, reinterpret_cast<const char*>(&opt), sizeof(opt));
        // setsockopt(server_socket, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
        if (bind(server_socket, (struct sockaddr *) &server_address, sizeof(server_address)) < 0) {
            cout << "Error: Failed to bind the socket." << endl;
            exit(EXIT_FAILURE);
        }


        if (listen(server_socket, 1) < 0) {
            cout << "Error: Failed to listen on the socket." << endl;
            exit(EXIT_FAILURE);
        }

        struct sockaddr_in client_address{};
        socklen_t client_address_size = sizeof(client_address);
        client_socket = accept(server_socket, (struct sockaddr *) &client_address, &client_address_size);

        if (client_socket < 0) {
            cout << "Error: Failed to accept client connection." << endl;
            exit(EXIT_FAILURE);
        }
    }


    /**
     * initialize client socket
     *
     * @param ip server IP
     * @param port server port
     * @return client socket
     */
    int client_socket_initial(char *ip, int port = 7766) {
        // create socket
        int sock;
        do {
            sock = socket(AF_INET, SOCK_STREAM, 0);
        } while (sock == -1);

        // set up server address and port
        struct sockaddr_in server_addr{};
        server_addr.sin_family = AF_INET;
        server_addr.sin_addr.s_addr = inet_addr(ip); // server IP
        server_addr.sin_port = htons(port); // server port

        // connect to server
        int connect_stats;
        do {
            connect_stats = connect(sock, (struct sockaddr *) &server_addr, sizeof(server_addr));
        } while (connect_stats == -1);

        // Get client socket information
        struct sockaddr_in client_addr{};
        socklen_t addr_len = sizeof(client_addr);
        getsockname(sock, (struct sockaddr *) &client_addr, &addr_len);

        // printf("Local socket address: %s:%d\n", inet_ntoa(client_addr.sin_addr), ntohs(client_addr.sin_port));

        // Set TCP_QUICKACK option to send ACK responses quickly
        int sockopt = 1;
        // int ret = setsockopt(sock, IPPROTO_TCP, TCP_QUICKACK, &sockopt, sizeof(sockopt));
        int ret = setsockopt(sock, IPPROTO_TCP, TCP_NODELAY, reinterpret_cast<const char*>(&sockopt), sizeof(sockopt));
        if (ret == -1) {
            perror("Couldn't setsockopt(TCP_QUICKACK)");
            exit(EXIT_FAILURE);
        }

        return sock;
    }

    void server_socket_close(int client_socket, int server_socket) {
        close(client_socket);
        close(server_socket);
    }

    void client_socket_close(int sock) {
        close(sock);
    }


}
#endif

