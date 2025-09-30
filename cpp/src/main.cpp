#include <atomic>
#include <chrono>
#include <iostream>
#include <mutex>
#include <thread>

#include "ws_server.hpp"
#include "board.hpp"
#include "proto.hpp"
#include "solver.hpp"
#include "overlay.hpp"

#ifdef _WIN32
#include <windows.h>
#endif

using namespace std::chrono_literals;

std::atomic<int> g_mines_total{0};
std::atomic<bool> g_enable_chords{true};

int main(){
#ifdef _WIN32
    HMODULE user32 = GetModuleHandleW(L"user32.dll");
    bool dpiSet = false;
    if(user32){
        typedef BOOL (WINAPI *SetProcDpiCtxFn)(DPI_AWARENESS_CONTEXT);
        auto setCtx = reinterpret_cast<SetProcDpiCtxFn>(GetProcAddress(user32, "SetProcessDpiAwarenessContext"));
        if(setCtx){
            if(setCtx(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2)){
                dpiSet = true;
            }
        }
    }
    if(!dpiSet){ SetProcessDPIAware(); }
#endif
#ifdef _WIN32
    const int kHotkeyExit = 1;
    const int kHotkeyToggleChords = 2;
    const int kHotkeyToggleCapture = 3;
    const int kHotkeyToggleSafety = 4;
    const int kHotkeyClickNearest = 5;
    RegisterHotKey(nullptr, kHotkeyExit, MOD_CONTROL | MOD_ALT | 0x4000, 'X');
    RegisterHotKey(nullptr, kHotkeyToggleChords, MOD_CONTROL | 0x4000, 'P');
    RegisterHotKey(nullptr, kHotkeyToggleCapture, MOD_CONTROL | 0x4000, 'W');
    RegisterHotKey(nullptr, kHotkeyToggleSafety, MOD_CONTROL | MOD_ALT | 0x4000, 'S');
    RegisterHotKey(nullptr, kHotkeyClickNearest, MOD_CONTROL | 0x4000, 'Q');
#endif
    std::ios::sync_with_stdio(false);
    std::cin.tie(nullptr);
    game::Board board;
    OverlayWindow overlay;
    overlay.create();

    net::WebSocketServer server;
    std::mutex mtx;
    std::atomic<bool> dirty{false};

    OverlayGeometry geom{};

    server.set_on_message([&](const net::WsMessage& msg){
        proto::ParsedMessage parsed; if(!proto::parse_message(msg.text, parsed)) return;
        std::lock_guard<std::mutex> lock(mtx);
        if(parsed.type==proto::MsgType::Full){
            board.apply_full(parsed.full.cells, parsed.full.w, parsed.full.h);
            int mt = parsed.full.mines_total; if(mt<0) mt = 0; g_mines_total.store(mt);
            geom.board_w = parsed.full.w; geom.board_h = parsed.full.h;
            geom.rect_l = parsed.full.rect_l; geom.rect_t = parsed.full.rect_t; geom.rect_w = parsed.full.rect_w; geom.rect_h = parsed.full.rect_h;
            geom.vv_x = parsed.full.vv_x; geom.vv_y = parsed.full.vv_y; geom.vv_scale = (parsed.full.vv_scale>0?parsed.full.vv_scale:1.0);
            geom.dpr = (parsed.full.dpr>0?parsed.full.dpr:1.0);
        } else if(parsed.type==proto::MsgType::Delta){
            board.apply_updates(parsed.delta.updates);
            int mt = parsed.delta.mines_total; if(mt<0) mt = 0; if(mt!=0) g_mines_total.store(mt);
            if(parsed.delta.rect_w>0 && parsed.delta.rect_h>0){
                geom.rect_l = parsed.delta.rect_l; geom.rect_t = parsed.delta.rect_t; geom.rect_w = parsed.delta.rect_w; geom.rect_h = parsed.delta.rect_h;
                geom.vv_x = parsed.delta.vv_x; geom.vv_y = parsed.delta.vv_y; geom.vv_scale = (parsed.delta.vv_scale>0?parsed.delta.vv_scale:geom.vv_scale);
                geom.dpr = (parsed.delta.dpr>0?parsed.delta.dpr:geom.dpr);
            }
        }
        dirty = true;
    });

    if(!server.start(8765)){
        std::cerr << "Failed to start WebSocket server on 8765" << std::endl;
        return 1;
    }

    bool running = true;
    for(; running;){
        MSG msg; while(PeekMessage(&msg, nullptr, 0, 0, PM_REMOVE)){
            if(msg.message == WM_HOTKEY){
                if((int)msg.wParam == kHotkeyExit){ running = false; break; }
                if((int)msg.wParam == kHotkeyToggleChords){ g_enable_chords.store(!g_enable_chords.load()); dirty = true; continue; }
                if((int)msg.wParam == kHotkeyToggleCapture){ overlay.set_excluded_from_capture(!overlay.is_excluded_from_capture()); continue; }
                if((int)msg.wParam == kHotkeyToggleSafety){ overlay.set_safety_mode(!overlay.is_safety_mode()); continue; }
                if((int)msg.wParam == kHotkeyClickNearest){ overlay.click_nearest_safe_guess(); continue; }
            }
            TranslateMessage(&msg); DispatchMessage(&msg);
        }
        bool didWork = false;
        if(dirty.load()){
            dirty=false;
            game::Board snapshot;
            OverlayGeometry gcopy;
            {
                std::lock_guard<std::mutex> lock(mtx);
                snapshot.apply_full(board.data(), board.width(), board.height());
                gcopy = geom;
            }
            int minesTotal = g_mines_total.load();
            bool enableChords = g_enable_chords.load();
            solve::Overlay ov = solve::compute_overlay(snapshot, minesTotal, enableChords);
            overlay.update(ov.marks, gcopy, minesTotal);
            didWork = true;
        }
        if(!didWork){ std::this_thread::sleep_for(2ms); }
        if(!didWork && g_enable_chords.load() == false){
            game::Board snapshot;
            OverlayGeometry gcopy;
            {
                std::lock_guard<std::mutex> lock(mtx);
                snapshot.apply_full(board.data(), board.width(), board.height());
                gcopy = geom;
            }
            int minesTotal = g_mines_total.load();
            solve::Overlay ov = solve::compute_overlay(snapshot, minesTotal, false);
            overlay.update(ov.marks, gcopy, minesTotal);
        }
    }
#ifdef _WIN32
    UnregisterHotKey(nullptr, kHotkeyExit);
    UnregisterHotKey(nullptr, kHotkeyToggleChords);
    UnregisterHotKey(nullptr, kHotkeyToggleCapture);
    UnregisterHotKey(nullptr, kHotkeyClickNearest);
#endif
    server.stop();
    overlay.destroy();
    return 0;
}

#ifdef _WIN32
int WINAPI wWinMain(HINSTANCE, HINSTANCE, PWSTR, int){
    return main();
}
#endif



