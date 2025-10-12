#include "ws_server.hpp"

#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <mutex>

namespace {
struct Sha1Ctx { uint32_t h[5]; uint64_t len; uint8_t buf[64]; size_t idx; };

static uint32_t rol(uint32_t v, int s){ return (v<<s) | (v>>(32-s)); }

static void sha1_init(Sha1Ctx& c){ c.h[0]=0x67452301; c.h[1]=0xEFCDAB89; c.h[2]=0x98BADCFE; c.h[3]=0x10325476; c.h[4]=0xC3D2E1F0; c.len=0; c.idx=0; }

static void sha1_block(Sha1Ctx& c){
    uint32_t w[80];
    for(int i=0;i<16;++i){ w[i] = (c.buf[i*4]<<24)|(c.buf[i*4+1]<<16)|(c.buf[i*4+2]<<8)|(c.buf[i*4+3]); }
    for(int i=16;i<80;++i){ w[i] = rol(w[i-3]^w[i-8]^w[i-14]^w[i-16],1); }
    uint32_t a=c.h[0],b=c.h[1],c2=c.h[2],d=c.h[3],e=c.h[4];
    for(int i=0;i<80;++i){
        uint32_t f,k;
        if(i<20){ f=(b&c2)|((~b)&d); k=0x5A827999; }
        else if(i<40){ f=b^c2^d; k=0x6ED9EBA1; }
        else if(i<60){ f=(b&c2)|(b&d)|(c2&d); k=0x8F1BBCDC; }
        else { f=b^c2^d; k=0xCA62C1D6; }
        uint32_t t = rol(a,5) + f + e + k + w[i];
        e=d; d=c2; c2=rol(b,30); b=a; a=t;
    }
    c.h[0]+=a; c.h[1]+=b; c.h[2]+=c2; c.h[3]+=d; c.h[4]+=e;
}

static void sha1_update(Sha1Ctx& c, const uint8_t* data, size_t len){
    c.len += len*8;
    for(size_t i=0;i<len;++i){
        c.buf[c.idx++] = data[i];
        if(c.idx==64){ sha1_block(c); c.idx=0; }
    }
}

static void sha1_final(Sha1Ctx& c, uint8_t out[20]){
    c.buf[c.idx++] = 0x80;
    if(c.idx>56){ while(c.idx<64) c.buf[c.idx++]=0; sha1_block(c); c.idx=0; }
    while(c.idx<56) c.buf[c.idx++]=0;
    for(int i=7;i>=0;--i){ c.buf[c.idx++] = (uint8_t)((c.len>>(i*8))&0xFF); }
    sha1_block(c);
    for(int i=0;i<5;++i){ out[i*4]=(c.h[i]>>24)&0xFF; out[i*4+1]=(c.h[i]>>16)&0xFF; out[i*4+2]=(c.h[i]>>8)&0xFF; out[i*4+3]=c.h[i]&0xFF; }
}

static std::string sha1_base64(const std::string& s){
    Sha1Ctx c; sha1_init(c); sha1_update(c, reinterpret_cast<const uint8_t*>(s.data()), s.size()); uint8_t d[20]; sha1_final(c,d);
    static const char* B64 = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    std::string out; out.reserve(28);
    int i=0; for(; i+2<20; i+=3){ uint32_t v=(d[i]<<16)|(d[i+1]<<8)|d[i+2]; out.push_back(B64[(v>>18)&63]); out.push_back(B64[(v>>12)&63]); out.push_back(B64[(v>>6)&63]); out.push_back(B64[v&63]); }
    int rem = 20 - i;
    if(rem==1){ uint32_t v = (d[i] << 16); out.push_back(B64[(v>>18)&63]); out.push_back(B64[(v>>12)&63]); out.push_back('='); out.push_back('='); }
    else if(rem==2){ uint32_t v = (d[i] << 16) | (d[i+1] << 8); out.push_back(B64[(v>>18)&63]); out.push_back(B64[(v>>12)&63]); out.push_back(B64[(v>>6)&63]); out.push_back('='); }
    return out;
}
}

namespace net {

WebSocketServer::WebSocketServer() : running_(false) {
#ifdef _WIN32
    WSADATA wsaData; WSAStartup(MAKEWORD(2,2), &wsaData);
#endif
}

WebSocketServer::~WebSocketServer() {
    stop();
#ifdef _WIN32
    WSACleanup();
#endif
}

void WebSocketServer::set_on_message(MessageCallback cb){ on_message_ = std::move(cb); }

bool WebSocketServer::start(uint16_t port){
    if(running_) return true;
    running_ = true;
#ifdef _WIN32
    listen_socket_ = socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if(listen_socket_ == INVALID_SOCKET) return false;

    u_long nonblock = 1; ioctlsocket(listen_socket_, FIONBIO, &nonblock);

    sockaddr_in addr{}; addr.sin_family = AF_INET; addr.sin_port = htons(port); addr.sin_addr.s_addr = htonl(INADDR_LOOPBACK);
    if(bind(listen_socket_, reinterpret_cast<sockaddr*>(&addr), sizeof(addr))==SOCKET_ERROR) return false;
    if(listen(listen_socket_, 1)==SOCKET_ERROR) return false;
#endif
    accept_thread_ = std::thread(&WebSocketServer::accept_loop, this);
    return true;
}

void WebSocketServer::stop(){
    running_ = false;
#ifdef _WIN32
    if(listen_socket_ != INVALID_SOCKET){ closesocket(listen_socket_); listen_socket_ = INVALID_SOCKET; }
    if(client_socket_ != INVALID_SOCKET){ closesocket(client_socket_); client_socket_ = INVALID_SOCKET; }
#endif
    if(accept_thread_.joinable()) accept_thread_.join();
}

bool WebSocketServer::perform_handshake(SOCKET s, const std::string& http_request){
    std::string lower = http_request; for(char& c : lower){ if(c>='A'&&c<='Z') c = char(c - 'A' + 'a'); }
    const std::string keyHdr = "sec-websocket-key:";
    auto pos = lower.find(keyHdr);
    if(pos==std::string::npos) return false;
    size_t colon = http_request.find(':', pos);
    if(colon==std::string::npos) return false;
    size_t lineEnd = http_request.find("\r\n", colon);
    if(lineEnd==std::string::npos) lineEnd = http_request.size();
    std::string key = http_request.substr(colon+1, lineEnd - (colon+1));
    size_t start = key.find_first_not_of(" \t"); size_t stop = key.find_last_not_of(" \t");
    if(start==std::string::npos) return false; key = key.substr(start, stop-start+1);
    const std::string magic = key + "258EAFA5-E914-47DA-95CA-C5AB0DC85B11";
    std::string accept_key = sha1_base64(magic);
    std::string resp =
        "HTTP/1.1 101 Switching Protocols\r\n"
        "Upgrade: websocket\r\n"
        "Connection: Upgrade\r\n"
        "Sec-WebSocket-Accept: " + accept_key + "\r\n\r\n";
    static const bool kWsDebug = false;
    int n = send(s, resp.c_str(), (int)resp.size(), 0);
    if(kWsDebug){ std::cerr << "ws: handshake resp bytes=" << n << std::endl; }
    return true;
}

void WebSocketServer::accept_loop(){
#ifdef _WIN32
    while(running_){
        fd_set rfds; FD_ZERO(&rfds); FD_SET(listen_socket_, &rfds);
        timeval tv{0, 200*1000};
        int r = select(0, &rfds, nullptr, nullptr, &tv);
        if(r>0 && FD_ISSET(listen_socket_, &rfds)){
            SOCKET s = accept(listen_socket_, nullptr, nullptr);
            if(s!=INVALID_SOCKET){
                static const bool kWsDebug = false;
                if(kWsDebug){ std::cerr << "ws: accept" << std::endl; }
                // make client socket blocking
                u_long nb0 = 0; ioctlsocket(s, FIONBIO, &nb0);
                // read HTTP request headers fully until CRLF CRLF
                std::string req; req.reserve(2048);
                char buf[1024];
                for(;;){
                    int n = recv(s, buf, sizeof(buf), 0);
                    if(n<=0) break;
                    req.append(buf, n);
                    if(req.find("\r\n\r\n") != std::string::npos) break;
                    if(req.size() > 16384) break; // sanity cap
                }
                if(kWsDebug){ std::cerr << "ws: received headers bytes=" << req.size() << std::endl; }
                if(!req.empty() && perform_handshake(s, req)){
                    if(kWsDebug){ std::cerr << "ws: handshake ok" << std::endl; }
                    client_socket_ = s;
                    std::thread(&WebSocketServer::client_loop, this).detach();
                } else {
                    if(kWsDebug){ std::cerr << "ws: handshake failed" << std::endl; }
                    closesocket(s);
                }
            }
        }
    }
#endif
}

static bool ws_read_frame(SOCKET s, std::string& out_text){
    std::string message_accum;
    uint8_t current_opcode = 0;
    bool in_fragmented = false;
    for(;;){
        uint8_t hdr[2]; int n = recv(s, reinterpret_cast<char*>(hdr), 2, 0); if(n!=2) return false;
        bool fin = (hdr[0] & 0x80)!=0; uint8_t opcode = hdr[0] & 0x0F; bool masked = (hdr[1] & 0x80)!=0; uint64_t len = hdr[1] & 0x7F;
        if(len==126){ uint8_t ext[2]; if(recv(s, reinterpret_cast<char*>(ext), 2, 0)!=2) return false; len = (ext[0]<<8)|ext[1]; }
        else if(len==127){ uint8_t ext[8]; if(recv(s, reinterpret_cast<char*>(ext), 8, 0)!=8) return false; len = 0; for(int i=0;i<8;++i){ len = (len<<8)|ext[i]; } }
        uint8_t mask[4]{}; if(masked){ if(recv(s, reinterpret_cast<char*>(mask), 4, 0)!=4) return false; }
        std::string payload; payload.resize(static_cast<size_t>(len));
        size_t got=0; while(got<len){ int m = recv(s, payload.data()+got, (int)(len-got), 0); if(m<=0) return false; got += m; }
        if(masked){ for(size_t i=0;i<len;++i){ payload[i] = payload[i] ^ mask[i%4]; } }
        // contorl opcodes
        if(opcode==0x8){ // close
            return false;
        } else if(opcode==0x9){ // ping -> reply pong
            uint8_t pong_hdr[2]; pong_hdr[0] = 0x80 | 0x0A; // FIN + pong
            if(len<126){ pong_hdr[1] = (uint8_t)len; send(s, reinterpret_cast<const char*>(pong_hdr), 2, 0); }
            else if(len<=0xFFFF){ pong_hdr[1] = 126; send(s, reinterpret_cast<const char*>(pong_hdr), 2, 0); uint8_t ext2[2]{ (uint8_t)((len>>8)&0xFF), (uint8_t)(len&0xFF) }; send(s, reinterpret_cast<const char*>(ext2), 2, 0); }
            else { pong_hdr[1] = 127; send(s, reinterpret_cast<const char*>(pong_hdr), 2, 0); uint8_t ext2[8]; for(int i=7;i>=0;--i){ ext2[7-i] = (uint8_t)((len>>(i*8))&0xFF); } send(s, reinterpret_cast<const char*>(ext2), 8, 0); }
            if(len>0) send(s, payload.data(), (int)len, 0);
            continue;
        } else if(opcode==0xA){ // pong
            continue;
        }
        // data opcodes: handle fragmentation
        if(opcode==0x1 && !in_fragmented){ // start text
            current_opcode = opcode; message_accum = std::move(payload); in_fragmented = !fin;
            if(fin){ out_text = std::move(message_accum); return true; }
        } else if(opcode==0x0 && in_fragmented){ // continuation
            message_accum += payload; if(fin){ out_text = std::move(message_accum); return true; }
        } else if(opcode==0x2 || (opcode==0x0 && !in_fragmented)){
            // ignore binary frames and unexpected continuation
            continue;
        } else if(opcode==0x1 && in_fragmented){
            // unexpected new text while fragmented; reset state and treat as new
            message_accum.clear(); in_fragmented=false; current_opcode=0; message_accum = std::move(payload); if(fin){ out_text = std::move(message_accum); return true; } else { in_fragmented=true; current_opcode=0x1; }
        }
    }
}

void WebSocketServer::client_loop(){
#ifdef _WIN32
    SOCKET s = client_socket_;
    for(;;){
        std::string text; if(!ws_read_frame(s, text)) break;
        if(on_message_) on_message_(WsMessage{std::move(text)});
    }
    closesocket(s);
    client_socket_ = INVALID_SOCKET;
#endif
}

bool WebSocketServer::send_text(const std::string& data){
#ifdef _WIN32
    if(client_socket_==INVALID_SOCKET) return false;
    std::string frame;
    frame.push_back((char)0x81); // FIN + text
    size_t len = data.size();
    if(len<126){ frame.push_back((char)len); }
    else if(len<=0xFFFF){ frame.push_back(126); frame.push_back((char)((len>>8)&0xFF)); frame.push_back((char)(len&0xFF)); }
    else { frame.push_back(127); for(int i=7;i>=0;--i) frame.push_back((char)((len>>(i*8))&0xFF)); }
    frame += data;
    int n = send(client_socket_, frame.c_str(), (int)frame.size(), 0);
    return n==(int)frame.size();
#else
    return false;
#endif
}
}