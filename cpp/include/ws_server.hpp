#pragma once

#include <cstdint>
#include <functional>
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

    WebSocketServer();
    ~WebSocketServer();

    bool start(uint16_t port);
    void stop();

    void set_on_message(MessageCallback cb);
    bool send_text(const std::string& data);

    bool is_connected() const {
#ifdef _WIN32
        return client_socket_ != INVALID_SOCKET;
#else
        return false;
#endif
    }

private:
    bool running_;
    MessageCallback on_message_;

#ifdef _WIN32
    SOCKET listen_socket_ = INVALID_SOCKET;
    SOCKET client_socket_ = INVALID_SOCKET;
#endif

    std::thread accept_thread_;

    bool perform_handshake(SOCKET s, const std::string& http_request);
    void accept_loop();
    void client_loop();
};

}
