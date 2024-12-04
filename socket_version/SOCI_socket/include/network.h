/******************************************************************************
 * Author: Aptx4869AC
 * Created: 2024-12-04
 * GitHub: https://github.com/Aptx4869AC
 *****************************************************************************/

#ifndef NETWORK_H
#define NETWORK_H
#pragma once

#include "paillier.h"

using namespace std;

namespace NETWORKSPACE {

    /**
     * 获取当前发送缓冲区大小
     * @param sockfd
     */
    void show_send_buffer_size(int sockfd) {
        int current_send_buffer_size;
        socklen_t optlen = sizeof(current_send_buffer_size);
        if (getsockopt(sockfd, SOL_SOCKET, SO_SNDBUF, &current_send_buffer_size, &optlen) == -1) {
            perror("Failed to get send buffer size");
            close(sockfd);
            return;
        }
        cout << "Current send buffer size: " << current_send_buffer_size << " bytes\n";

    }

    /**
     * 获取当前接收缓冲区大小
     * @param sockfd
     */
    void show_recv_buffer_size(int sockfd) {
        int current_recv_buffer_size;
        socklen_t optlen = sizeof(current_recv_buffer_size);
        if (getsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, &current_recv_buffer_size, &optlen) == -1) {
            perror("Failed to get receive buffer size");
            close(sockfd);
            return;
        }
        cout << "Current receive buffer size: " << current_recv_buffer_size << " bytes\n";
    }

    /**
     * 设置新的发送缓冲区大小和设置新的接收缓冲区大小
     * @param sockfd
     * @param new_send_buffer_size
     * @param new_recv_buffer_size
     * @return
     */
    bool setSocketBufferSize(int sockfd, int new_send_buffer_size, int new_recv_buffer_size) {

        if (setsockopt(sockfd, SOL_SOCKET, SO_SNDBUF, &new_send_buffer_size, sizeof(new_send_buffer_size)) == -1) {
            perror("Failed to set send buffer size");
            close(sockfd);
            return false;
        }

        if (setsockopt(sockfd, SOL_SOCKET, SO_RCVBUF, &new_recv_buffer_size, sizeof(new_recv_buffer_size)) == -1) {
            perror("Failed to set receive buffer size");
            close(sockfd);
            return false;
        }

        show_send_buffer_size(sockfd);
        show_recv_buffer_size(sockfd);
        return true;
    }

    /**
     * 传输到文件以记录开销
     * @param filename
     * @param src
     * @return
     */
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

    /**
     * 序列化mpz
     * @param sock
     * @param num
     */
    void send_mpz(int sock, const mpz_t num) {
        size_t count = (mpz_sizeinbase(num, 2) + 7) / 8;
        unsigned char *buffer = new unsigned char[count];
        mpz_export(buffer, &count, 1, 1, 0, 0, num);
        // 发送数据长度
        send(sock, &count, sizeof(size_t), 0);
        // 发送数据内容
        send(sock, buffer, count, 0);
        delete[] buffer;
    }

    /**
     * 反序列化mpz
     * @param sock
     * @param num
     */
    void recv_mpz(int sock, mpz_t num) {
        mpz_init(num);
        size_t count;
        recv(sock, &count, sizeof(size_t), 0);
        unsigned char *buffer = new unsigned char[count];
        recv(sock, buffer, count, 0);
        mpz_import(num, count, 1, 1, 0, 0, buffer);
        delete[] buffer;
    }

    /**
     * 以mpz数组形式序列化
     * @param sock
     * @param src
     * @param num_count
     */
    void send_mpz_array(int sock, const vector <mpz_ptr> &src, size_t num_count) {
        double start_time, end_time;
        start_time = omp_get_wtime();

        size_t *sizes = new size_t[num_count];

        // 计算每个元素的大小并总计所有数据的总大小
        size_t total_data_size = 0;
        for (size_t i = 0; i < num_count; ++i) {
            sizes[i] = (mpz_sizeinbase(src[i], 2) + 7) / 8;
            total_data_size += sizes[i];
        }

        // 第一个缓冲区包含num_count和sizes
        size_t buffer1_size = num_count * sizeof(size_t);
        unsigned char *buffer1 = new unsigned char[buffer1_size];
        unsigned char *ptr1 = buffer1;

        // 将每个元素的大小拷贝到缓冲区
        for (size_t i = 0; i < num_count; ++i) {
            memcpy(ptr1, &sizes[i], sizeof(size_t));
            ptr1 += sizeof(size_t);
        }

        // 第二个缓冲区包含所有的mpz_t数据
        unsigned char *buffer2 = new unsigned char[total_data_size];
        unsigned char *ptr2 = buffer2;

        // 将每个元素的数据拷贝到缓冲区
        for (size_t i = 0; i < num_count; ++i) {
            mpz_export(ptr2, &sizes[i], 1, 1, 0, 0, src[i]);
            ptr2 += sizes[i];
        }

        // 发送第一个缓冲区
        send(sock, buffer1, buffer1_size, 0);

        // 发送第二个缓冲区
        send(sock, buffer2, total_data_size, 0);

        delete[] sizes;
        delete[] buffer1;
        delete[] buffer2;

        end_time = omp_get_wtime();
    }

    /**
     * 以mpz数组形式反序列化
     * @param sock
     * @param dst
     * @param num_count
     */
    void receive_mpz_array(int sock, vector <mpz_ptr> &dst, size_t num_count) {
        double start_time, end_time;
        start_time = omp_get_wtime();

        size_t *sizes = new size_t[num_count];
        recv(sock, sizes, num_count * sizeof(size_t), 0);

        // 初始化 dst 并为每个 mpz_t 分配空间
        dst.resize(num_count);
        for (size_t i = 0; i < num_count; ++i) {
            mpz_init(dst[i]);
        }

        // 计算第二个缓冲区的大小并接收它
        size_t total_data_size = 0;
        for (size_t i = 0; i < num_count; ++i) {
            total_data_size += sizes[i];
        }

        unsigned char *buffer2 = new unsigned char[total_data_size];
        recv(sock, buffer2, total_data_size, 0);

        // 将每个元素的数据从缓冲区导入到 mpz_t
        unsigned char *ptr2 = buffer2;
        for (size_t i = 0; i < num_count; ++i) {
            mpz_init(dst[i]);
            mpz_import(dst[i], sizes[i], 1, 1, 0, 0, ptr2);
            ptr2 += sizes[i];
        }

        delete[] sizes;
        delete[] buffer2;

        end_time = omp_get_wtime();
//        printf("接收开销 ------  %f ms\n", (end_time - start_time) * 1000);
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
        setsockopt(server_socket, SOL_SOCKET, SO_REUSEADDR, &opt, sizeof(opt));
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
        int ret = setsockopt(sock, IPPROTO_TCP, TCP_QUICKACK, &sockopt, sizeof(sockopt));
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

