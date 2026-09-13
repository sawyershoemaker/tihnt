#include "ws_server.hpp"
#include <atomic>
#include <iostream>
#include <stdexcept>
#include <string>

int main(int argc, char** argv) {
    if(argc != 2) return 2;
    net::WebSocketServer server;
    std::atomic<bool> throwConnections{false};
    server.set_on_connection([&](bool) {
        if(throwConnections) throw 42;
    });
    const auto port = static_cast<uint16_t>(std::stoi(argv[1]));
    if(!server.start(port)) return 3;
    {
        net::WebSocketServer other;
        if(other.start(port)) return 4;
        if(!other.start(0)) return 5; // A failed bind must be recoverable.
        other.stop(); other.stop();
    }
    server.set_on_message([&](const net::WsMessage& message) {
        if(message.text == "throw") throw std::runtime_error("test callback error");
        if(message.text == "throw nonstandard") throw 42;
        if(message.text == "throw connections") throwConnections=true;
        if(message.text == "flood") {
            const std::string payload(1024*1024,'x');
            std::cout<<"SENDING\n"<<std::flush;
            for(int i=0;i<64;++i) if(!server.send_text(payload)) break;
            std::cout<<"FLOOD_DONE\n"<<std::flush;
            return;
        }
        server.send_text(message.text);
    });
    std::cout << "READY\n" << std::flush;
    std::string command;
    std::getline(std::cin, command);
    server.stop(); server.stop();
    return 0;
}
