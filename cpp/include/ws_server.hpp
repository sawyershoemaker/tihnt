#pragma once

#include <cstdint>
#include <atomic>
#include <functional>
#include <mutex>
#include <string>
#include <thread>
#include <vector>

#ifdef _WIN32
#  include <winsock2.h>
#  include <ws2tcpip.h>
#endif

namespace net {

struct WsMessage {
    std::string text;
};

class WebSocketServer {
public:
    using MessageCallback = std::function<void(const WsMessage&)>;
    using ConnectionCallback = std::function<void(bool)>;

    WebSocketServer();
    ~WebSocketServer();

    bool start(uint16_t port);
    void stop();

    void set_on_message(MessageCallback cb);
    void set_on_connection(ConnectionCallback cb);
    bool send_text(const std::string& data);

    bool is_connected() const {
#ifdef _WIN32
        std::lock_guard<std::mutex> lock(client_mutex_);
        return client_socket_ != INVALID_SOCKET;
#else
        return false;
#endif
    }

private:
    std::atomic<bool> running_{false};
    MessageCallback on_message_;
    ConnectionCallback on_connection_;
    std::mutex callback_mutex_;
    std::mutex send_mutex_;
    void notify_connection(bool connected);

#ifdef _WIN32
    SOCKET listen_socket_ = INVALID_SOCKET;
    SOCKET client_socket_ = INVALID_SOCKET;
    SOCKET handshake_socket_ = INVALID_SOCKET;
    bool winsock_ready_ = false;
    mutable std::mutex client_mutex_;
    bool send_frame(SOCKET s, uint8_t opcode, const std::string& data);
    bool read_message(SOCKET s, std::string& buffered, std::string& text);
#endif

    std::thread accept_thread_;
    std::thread client_thread_;

    void accept_loop();
#ifdef _WIN32
    bool perform_handshake(SOCKET s, const std::string& http_request);
    void client_loop(SOCKET s, std::string buffered);
#endif
};

}
