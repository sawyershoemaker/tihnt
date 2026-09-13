#include "ws_server.hpp"

#include <array>
#include <atomic>
#include <cstdint>
#include <cstring>
#include <iostream>
#include <mutex>
#include <algorithm>
#include <chrono>
#include <string_view>
#include <unordered_map>
#include "proto.hpp"
#include "utf8.hpp"

namespace {
struct Sha1Ctx { uint32_t h[5]; uint64_t len; uint8_t buf[64]; size_t idx; };

static uint32_t rol(uint32_t v, int s){ return (v<<s) | (v>>(32-s)); }

static void sha1_init(Sha1Ctx& c){ c.h[0]=0x67452301; c.h[1]=0xEFCDAB89; c.h[2]=0x98BADCFE; c.h[3]=0x10325476; c.h[4]=0xC3D2E1F0; c.len=0; c.idx=0; }

static void sha1_block(Sha1Ctx& c){
    uint32_t w[80];
    for(int i=0;i<16;++i){ w[i] = (uint32_t(c.buf[i*4])<<24)|(uint32_t(c.buf[i*4+1])<<16)|(uint32_t(c.buf[i*4+2])<<8)|uint32_t(c.buf[i*4+3]); }
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

WebSocketServer::WebSocketServer() {
#ifdef _WIN32
    WSADATA data{};
    winsock_ready_ = WSAStartup(MAKEWORD(2,2), &data) == 0;
#endif
}

WebSocketServer::~WebSocketServer() {
    stop();
#ifdef _WIN32
    if(winsock_ready_) WSACleanup();
#endif
}

void WebSocketServer::set_on_message(MessageCallback cb) {
    std::lock_guard<std::mutex> lock(callback_mutex_);
    on_message_ = std::move(cb);
}
void WebSocketServer::set_on_connection(ConnectionCallback cb) {
    std::lock_guard<std::mutex> lock(callback_mutex_);
    on_connection_ = std::move(cb);
}
void WebSocketServer::notify_connection(bool connected) {
    ConnectionCallback cb;
    { std::lock_guard<std::mutex> lock(callback_mutex_); cb=on_connection_; }
    try { if(cb) cb(connected); }
    catch(...) { std::cerr<<"WebSocket connection callback failed\n"; }
}

bool WebSocketServer::start(uint16_t port) {
#ifdef _WIN32
    if(running_) return true;
    if(!winsock_ready_) return false;
    listen_socket_=socket(AF_INET, SOCK_STREAM, IPPROTO_TCP);
    if(listen_socket_==INVALID_SOCKET) return false;
    const BOOL exclusive=TRUE;
    setsockopt(listen_socket_, SOL_SOCKET, SO_EXCLUSIVEADDRUSE, reinterpret_cast<const char*>(&exclusive), sizeof(exclusive));
    sockaddr_in addr{};
    addr.sin_family=AF_INET; addr.sin_port=htons(port); addr.sin_addr.s_addr=htonl(INADDR_LOOPBACK);
    if(bind(listen_socket_, reinterpret_cast<sockaddr*>(&addr), sizeof(addr))==SOCKET_ERROR ||
        listen(listen_socket_, SOMAXCONN)==SOCKET_ERROR) {
        closesocket(listen_socket_); listen_socket_=INVALID_SOCKET; return false;
    }
    running_=true;
    try { accept_thread_=std::thread(&WebSocketServer::accept_loop,this); }
    catch(...) { running_=false; closesocket(listen_socket_); listen_socket_=INVALID_SOCKET; throw; }
    return true;
#else
    (void)port; return false;
#endif
}

void WebSocketServer::stop() {
    running_=false;
#ifdef _WIN32
    client_stopping_=true;
    // Interrupt I/O, but let each owner close its socket after it stops
    // using it. Closing and reusing a SOCKET during recv races a new connection.
    {
        std::lock_guard<std::mutex> lock(client_mutex_);
        if(client_socket_!=INVALID_SOCKET) shutdown(client_socket_,SD_BOTH);
        if(handshake_socket_!=INVALID_SOCKET) shutdown(handshake_socket_,SD_BOTH);
    }
#endif
    if(accept_thread_.joinable()) accept_thread_.join();
    if(client_thread_.joinable()) client_thread_.join();
#ifdef _WIN32
    if(listen_socket_!=INVALID_SOCKET) { closesocket(listen_socket_); listen_socket_=INVALID_SOCKET; }
#endif
}

#ifdef _WIN32
namespace {
using Clock=std::chrono::steady_clock;

bool wait_socket(SOCKET s,bool writing,Clock::time_point deadline,const std::atomic<bool>& running,
                 const std::atomic<bool>* interrupted=nullptr) {
    while(running && (!interrupted || !*interrupted)) {
        const auto now=Clock::now();
        if(now>=deadline) return false;
        const auto remaining=deadline==Clock::time_point::max() ? std::chrono::microseconds(100000) :
            std::min(std::chrono::duration_cast<std::chrono::microseconds>(deadline-now),std::chrono::microseconds(100000));
        timeval timeout{0,static_cast<long>(remaining.count())};
        fd_set ready; FD_ZERO(&ready); FD_SET(s,&ready);
        const int result=select(0,writing ? nullptr : &ready,writing ? &ready : nullptr,nullptr,&timeout);
        if(result>0) return running && (!interrupted || !*interrupted);
        if(result==SOCKET_ERROR) return false;
    }
    return false;
}

bool send_all(SOCKET s,const std::string& data,const std::atomic<bool>& running,
              const std::atomic<bool>* interrupted=nullptr) {
    size_t sent=0;
    const auto deadline=Clock::now()+std::chrono::seconds(5);
    while(sent<data.size()) {
        if(!running || (interrupted && *interrupted) || Clock::now()>=deadline) return false;
        const int n=send(s,data.data()+sent,static_cast<int>(data.size()-sent),0);
        if(n==SOCKET_ERROR && WSAGetLastError()==WSAEWOULDBLOCK) {
            if(!wait_socket(s,true,deadline,running,interrupted)) return false;
            continue;
        }
        if(n<=0) return false;
        sent+=n;
    }
    return true;
}
std::string lower(std::string text) {
    for(char& c:text) if(c>='A' && c<='Z') c=static_cast<char>(c-'A'+'a');
    return text;
}
std::string trim(std::string text) {
    const auto first=text.find_first_not_of(" \t"), last=text.find_last_not_of(" \t");
    return first==std::string::npos ? "" : text.substr(first,last-first+1);
}
bool has_token(const std::string& text, std::string_view wanted) {
    size_t pos=0;
    do {
        auto end=text.find(',',pos);
        if(trim(lower(text.substr(pos,end-pos)))==wanted) return true;
        if(end==std::string::npos) return false;
        pos=end+1;
    } while(pos<text.size());
    return false;
}
bool extension_origin(const std::string& origin) {
    const std::string prefix="chrome-extension://";
    return origin.size()==prefix.size()+32 && origin.compare(0,prefix.size(),prefix)==0 &&
        std::all_of(origin.begin()+prefix.size(),origin.end(),[](char c){ return c>='a' && c<='p'; });
}
bool valid_key(const std::string& key) {
    if(key.size()!=24 || key.substr(22)!="==") return false;
    const std::string b64="ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789+/";
    for(size_t i=0;i<22;++i) if(b64.find(key[i])==std::string::npos) return false;
    return b64.find(key[21])%16==0; // A canonical base64 encoding of 16 bytes.
}
bool valid_close(const std::string& payload) {
    if(payload.empty()) return true;
    if(payload.size()==1) return false;
    unsigned code=(static_cast<unsigned char>(payload[0])<<8)|static_cast<unsigned char>(payload[1]);
    const bool valid=(code>=3000 && code<=4999) || (code>=1000 && code<=1014 && code!=1004 && code!=1005 && code!=1006);
    return valid && valid_utf8(std::string_view(payload).substr(2));
}
}

bool WebSocketServer::perform_handshake(SOCKET s,const std::string& request) {
    const size_t firstEnd=request.find("\r\n");
    if(firstEnd==std::string::npos || request.substr(0,firstEnd)!="GET / HTTP/1.1") return false;
    std::unordered_map<std::string,std::string> headers;
    for(size_t pos=firstEnd+2; pos<request.size();) {
        size_t end=request.find("\r\n",pos);
        if(end==std::string::npos) return false;
        if(end==pos) break;
        size_t colon=request.find(':',pos);
        if(colon==std::string::npos || colon>=end) return false;
        std::string name=lower(request.substr(pos,colon-pos));
        if(name.empty() || name.find_first_of(" \t")!=std::string::npos) return false;
        if(!headers.emplace(name,trim(request.substr(colon+1,end-colon-1))).second) return false;
        pos=end+2;
    }
    if(!extension_origin(headers["origin"]) || lower(headers["upgrade"])!="websocket" ||
       !has_token(headers["connection"],"upgrade") || headers["sec-websocket-version"]!="13" ||
       !valid_key(headers["sec-websocket-key"])) return false;
    // Validate the actual bound port, including ephemeral ports used by tests.
    sockaddr_in addr{}; int length=sizeof(addr);
    if(getsockname(s,reinterpret_cast<sockaddr*>(&addr),&length)!=0) return false;
    const std::string port=":"+std::to_string(ntohs(addr.sin_port));
    if(headers["host"]!="127.0.0.1"+port && headers["host"]!="localhost"+port) return false;
    const std::string response="HTTP/1.1 101 Switching Protocols\r\nUpgrade: websocket\r\nConnection: Upgrade\r\nSec-WebSocket-Accept: "+
        sha1_base64(headers["sec-websocket-key"]+"258EAFA5-E914-47DA-95CA-C5AB0DC85B11")+"\r\n\r\n";
    return send_all(s,response,running_);
}

void WebSocketServer::accept_loop() {
    while(running_) {
        fd_set ready; FD_ZERO(&ready); FD_SET(listen_socket_,&ready);
        timeval timeout{0,100000};
        int result=select(0,&ready,nullptr,nullptr,&timeout);
        if(result<=0 || !running_) continue;
        SOCKET s=accept(listen_socket_,nullptr,nullptr);
        if(s==INVALID_SOCKET) continue;
        {
            std::lock_guard<std::mutex> lock(client_mutex_);
            if(!running_) { closesocket(s); break; }
            handshake_socket_=s;
        }
        u_long nonblocking=1;
        if(ioctlsocket(s,FIONBIO,&nonblocking)!=0) {
            std::lock_guard<std::mutex> lock(client_mutex_);
            handshake_socket_=INVALID_SOCKET;
            closesocket(s);
            continue;
        }
        const BOOL noDelay=TRUE;
        setsockopt(s,IPPROTO_TCP,TCP_NODELAY,reinterpret_cast<const char*>(&noDelay),sizeof(noDelay));
        std::string request;
        size_t headerEnd=std::string::npos;
        const auto deadline=Clock::now()+std::chrono::seconds(3);
        while(running_ && request.size()<16384 && Clock::now()<deadline) {
            if(!wait_socket(s,false,deadline,running_)) break;
            char buffer[2048];
            int n=recv(s,buffer,static_cast<int>(std::min<size_t>(sizeof(buffer),16384-request.size())),0);
            if(n==SOCKET_ERROR && WSAGetLastError()==WSAEWOULDBLOCK) continue;
            if(n<=0) break;
            request.append(buffer,n);
            headerEnd=request.find("\r\n\r\n");
            if(headerEnd!=std::string::npos) break;
        }
        bool accepted=running_ && headerEnd!=std::string::npos && perform_handshake(s,request.substr(0,headerEnd+4));
        if(accepted) {
            client_stopping_=true;
            {
                std::lock_guard<std::mutex> lock(client_mutex_);
                if(client_socket_!=INVALID_SOCKET) shutdown(client_socket_,SD_BOTH);
            }
            if(client_thread_.joinable()) client_thread_.join();
        }
        {
            std::lock_guard<std::mutex> lock(client_mutex_);
            handshake_socket_=INVALID_SOCKET;
            if(!accepted || !running_) { closesocket(s); continue; }
            client_socket_=s;
            client_stopping_=false;
        }
        notify_connection(true);
        client_thread_=std::thread(&WebSocketServer::client_loop,this,s,request.substr(headerEnd+4));
    }
}

bool WebSocketServer::send_frame(SOCKET s,uint8_t opcode,const std::string& data) {
    std::lock_guard<std::mutex> lock(send_mutex_);
    return send_frame_locked(s,opcode,data);
}

bool WebSocketServer::send_frame_locked(SOCKET s,uint8_t opcode,const std::string& data) {
    std::string frame;
    const uint64_t length=data.size();
    frame.reserve(data.size()+10);
    frame.push_back(static_cast<char>(0x80|opcode));
    if(length<126) frame.push_back(static_cast<char>(length));
    else if(length<=65535) {
        frame.push_back(126); frame.push_back(static_cast<char>(length>>8)); frame.push_back(static_cast<char>(length));
    } else {
        frame.push_back(127);
        for(int i=7;i>=0;--i) frame.push_back(static_cast<char>(length>>(i*8)));
    }
    frame+=data;
    if(send_all(s,frame,running_,&client_stopping_)) return true;
    // A partial frame cannot be followed by another frame on this stream.
    client_stopping_=true;
    shutdown(s,SD_BOTH);
    return false;
}

bool WebSocketServer::read_message(SOCKET s,ReceiveBuffer& buffered,std::string& text) {
    text.clear();
    bool fragmented=false;
    auto deadline=Clock::time_point::max();
    auto read_exact = [&](char* dst,size_t size) {
        size_t got=0;
        while(got<size && running_ && !client_stopping_) {
            if(Clock::now()>=deadline) return false;
            size_t available=std::min(size-got,buffered.bytes.size()-buffered.consumed);
            if(available>0) {
                std::memcpy(dst+got,buffered.bytes.data()+buffered.consumed,available);
                buffered.consumed+=available; got+=available;
            } else {
                std::array<char,8192> chunk;
                const bool direct=size-got>=chunk.size();
                char* destination=direct ? dst+got : chunk.data();
                const size_t capacity=direct ? size-got : chunk.size();
                int n=recv(s,destination,static_cast<int>(capacity),0);
                if(n==SOCKET_ERROR && WSAGetLastError()==WSAEWOULDBLOCK) {
                    if(!wait_socket(s,false,deadline,running_,&client_stopping_)) return false;
                    continue;
                }
                if(n<=0) return false;
                if(direct) got+=n;
                else {
                    buffered.bytes.assign(chunk.data(),static_cast<size_t>(n));
                    buffered.consumed=0;
                }
            }
            if(deadline==Clock::time_point::max()) deadline=Clock::now()+std::chrono::seconds(5);
        }
        return got==size;
    };
    auto fail = [&](unsigned code) {
        send_frame(s,8,std::string{static_cast<char>(code>>8),static_cast<char>(code)});
        return false;
    };
    for(;;) {
        uint8_t header[2];
        if(!read_exact(reinterpret_cast<char*>(header),2)) return false;
        bool fin=(header[0]&0x80)!=0;
        const uint8_t opcode=header[0]&15;
        const bool control=(opcode&8)!=0;
        uint64_t length=header[1]&127;
        if((header[0]&0x70)!=0 || !(header[1]&0x80) ||
           (opcode!=0 && opcode!=1 && opcode!=8 && opcode!=9 && opcode!=10) ||
           (control && (!fin || length>125))) return fail(1002);
        const unsigned marker=static_cast<unsigned>(length);
        if(marker==126 || marker==127) {
            uint8_t ext[8]{};
            const size_t bytes=marker==126 ? 2 : 8;
            if(!read_exact(reinterpret_cast<char*>(ext),bytes)) return false;
            if(bytes==8 && (ext[0]&0x80)) return fail(1002);
            length=0; for(size_t i=0;i<bytes;++i) length=(length<<8)|ext[i];
            if((marker==126 && length<126) || (marker==127 && length<=65535)) return fail(1002);
        }
        if(length>proto::MaxMessageBytes || (!control && length>proto::MaxMessageBytes-text.size())) return fail(1009);
        if(!control && ((opcode==0 && !fragmented) || (opcode==1 && fragmented))) return fail(1002);
        uint8_t mask[4];
        if(!read_exact(reinterpret_cast<char*>(mask),4)) return false;
        std::string payload(static_cast<size_t>(length),'\0');
        if(!read_exact(payload.data(),payload.size())) return false;
        for(size_t i=0;i<payload.size();++i) payload[i]^=mask[i%4];
        if(opcode==8) {
            if(!valid_close(payload)) return fail(1002);
            send_frame(s,8,payload); return false;
        }
        if(opcode==9 || opcode==10) {
            if(opcode==9 && !send_frame(s,10,payload)) return false;
            if(!fragmented) deadline=Clock::time_point::max();
            continue;
        }
        text+=payload;
        fragmented=!fin;
        if(fin) return valid_utf8(text) ? true : fail(1007);
    }
}

void WebSocketServer::client_loop(SOCKET s,std::string buffered) {
    ReceiveBuffer input{std::move(buffered),0};
    try {
        while(running_ && !client_stopping_) {
            std::string text;
            if(!read_message(s,input,text)) break;
            MessageCallback cb;
            { std::lock_guard<std::mutex> lock(callback_mutex_); cb=on_message_; }
            if(cb && running_ && !client_stopping_) cb(WsMessage{std::move(text)});
        }
    } catch(const std::exception& e) {
        std::cerr<<"WebSocket client: "<<e.what()<<'\n';
    } catch(...) {
        std::cerr<<"WebSocket client callback failed\n";
    }
    {
        std::lock_guard<std::mutex> sendLock(send_mutex_);
        std::lock_guard<std::mutex> lock(client_mutex_);
        closesocket(s);
        if(client_socket_==s) client_socket_=INVALID_SOCKET;
    }
    notify_connection(false);
}
#else
void WebSocketServer::accept_loop() {}
#endif

bool WebSocketServer::send_text(const std::string& data) {
#ifdef _WIN32
    if(data.size()>proto::MaxMessageBytes || !valid_utf8(data)) return false;
    std::lock_guard<std::mutex> sendLock(send_mutex_);
    SOCKET s;
    {
        std::lock_guard<std::mutex> lock(client_mutex_);
        s=client_socket_;
    }
    return running_ && !client_stopping_ && s!=INVALID_SOCKET && send_frame_locked(s,1,data);
#else
    (void)data; return false;
#endif
}
} // namespace net
