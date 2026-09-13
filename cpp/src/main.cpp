#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <string>
#include <thread>
#include <utility>

#include "ws_server.hpp"
#include "board.hpp"
#include "proto.hpp"
#include "solver.hpp"
#include "overlay.hpp"
#include "session.hpp"

namespace {
enum HotkeyId { Exit=1, Chords=2, Capture=3, Safety=4 };
constexpr struct { int id; UINT modifiers, key; const char* name; } hotkeyCommands[] = {
    {Exit,MOD_CONTROL|MOD_ALT|MOD_NOREPEAT,'X',"Ctrl + Alt + X"},
    {Chords,MOD_CONTROL|MOD_NOREPEAT,'P',"Ctrl + P"},
    {Capture,MOD_CONTROL|MOD_ALT|MOD_SHIFT|MOD_NOREPEAT,'W',"Ctrl + Alt + Shift + W"},
    {Safety,MOD_CONTROL|MOD_ALT|MOD_NOREPEAT,'S',"Ctrl + Alt + S"}
};

struct Hotkeys {
    int registered=0;
    ~Hotkeys() { release(); }
    void release() {
        while(registered>0) UnregisterHotKey(nullptr,hotkeyCommands[--registered].id);
    }
    bool create() {
        for(const auto& command : hotkeyCommands) {
            if(!RegisterHotKey(nullptr,command.id,command.modifiers,command.key)) {
                const DWORD error=GetLastError();
                release();
                const std::string name=command.name;
                std::cerr<<"Failed to register "<<name<<", Win32 error "<<error<<'\n';
                const std::wstring message=L"Could not register "+std::wstring(name.begin(),name.end())+
                    L".\nTIHNT cannot start without its keyboard controls.";
                MessageBoxW(nullptr,message.c_str(),L"TIHNT",MB_OK|MB_ICONERROR);
                return false;
            }
            ++registered;
        }
        return true;
    }
};

OverlayGeometry geometry_for(const proto::GeometryMsg& message, int w, int h) {
    OverlayGeometry g{};
    g.board_w=w; g.board_h=h;
    g.rect_l=message.rect_l; g.rect_t=message.rect_t;
    g.rect_w=message.rect_w; g.rect_h=message.rect_h;
    g.has_clip=message.has_clip;
    g.clip_l=message.clip_l; g.clip_t=message.clip_t;
    g.clip_w=message.clip_w; g.clip_h=message.clip_h;
    g.vv_x=message.vv_x; g.vv_y=message.vv_y;
    g.vv_scale=message.vv_scale; g.dpr=message.dpr;
    return g;
}
}

int main() {
    HMODULE user32=GetModuleHandleW(L"user32.dll");
    using SetDpiContext=BOOL(WINAPI*)(DPI_AWARENESS_CONTEXT);
    auto setContext=user32 ? reinterpret_cast<SetDpiContext>(GetProcAddress(user32,"SetProcessDpiAwarenessContext")) : nullptr;
    if(!setContext || !setContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2)) SetProcessDPIAware();

    Hotkeys hotkeys;
    if(!hotkeys.create()) return 1;
    OverlayWindow overlay;
    if(!overlay.create()) return 1;
    HANDLE wake=CreateEventW(nullptr,FALSE,FALSE,nullptr);
    if(!wake) return 1;

    std::mutex stateMutex;
    app::Session session;
    std::atomic<bool> stopping{false};
    auto updatePresentation=[&] {
        if(!session.take_presentation_changed()) return;
        const auto& view=session.presentation();
        overlay.update(view.marks,geometry_for(view.geometry,view.width,view.height),view.minesTotal);
    };

    net::WebSocketServer server;
    server.set_on_connection([&](bool) {
        {
            std::lock_guard<std::mutex> lock(stateMutex);
            session.reset();
        }
        SetEvent(wake);
    });
    server.set_on_message([&](const net::WsMessage& msg) {
        proto::ParsedMessage parsed;
        if(!proto::parse_message(msg.text,parsed)) return;
        {
            std::lock_guard<std::mutex> lock(stateMutex);
            if(!session.apply(std::move(parsed))) return;
        }
        SetEvent(wake);
    });
    if(!server.start(8765)) {
        std::cerr<<"Failed to start WebSocket server on 127.0.0.1:8765\n";
        CloseHandle(wake); return 1;
    }

    std::mutex jobMutex,resultMutex;
    std::condition_variable jobCv;
    app::SolverJob pendingJob;
    app::SolverResult latestResult;
    bool jobReady=false,resultReady=false;
    std::thread solver([&] {
        for(;;) {
            app::SolverJob job;
            {
                std::unique_lock<std::mutex> lock(jobMutex);
                jobCv.wait(lock,[&] { return stopping || jobReady; });
                if(stopping) break;
                job=std::move(pendingJob); jobReady=false;
            }
            auto cancelled=[&] { return stopping.load() || session.current_revision()!=job.revision; };
            if(cancelled()) continue;
            solve::Overlay result;
            try { result=solve::compute_overlay(job.board,job.minesTotal,job.enableChords,0,cancelled); }
            catch(const std::exception& e) {
                std::cerr<<"Solver failed: "<<e.what()<<'\n';
                result.marks.assign(job.board.data().size(),solve::Mark::None);
            }
            if(cancelled()) continue;
            {
                std::lock_guard<std::mutex> lock(resultMutex);
                latestResult={std::move(result.marks),job.revision};
                resultReady=true;
            }
            SetEvent(wake);
        }
    });

    bool running=true;
    while(running) {
        MSG msg;
        for(int dispatched=0;dispatched<64 && PeekMessageW(&msg,nullptr,0,0,PM_REMOVE);++dispatched) {
            if(msg.message==WM_QUIT) { running=false; break; }
            if(msg.message==WM_HOTKEY) {
                if(msg.wParam==Exit) { running=false; break; }
                if(msg.wParam==Chords) {
                    std::lock_guard<std::mutex> lock(stateMutex);
                    session.set_chords(!session.chords_enabled());
                }
                if(msg.wParam==Capture) overlay.set_excluded_from_capture(!overlay.is_excluded_from_capture());
                if(msg.wParam==Safety) overlay.set_safety_mode(!overlay.is_safety_mode());
                continue;
            }
            TranslateMessage(&msg); DispatchMessageW(&msg);
        }
        if(!running) break;
        {
            std::lock_guard<std::mutex> lock(stateMutex);
            if(auto binding=session.take_binding()) overlay.set_target_pid(*binding);
            if(auto job=session.take_job()) {
                {
                    std::lock_guard<std::mutex> jobLock(jobMutex);
                    pendingJob=std::move(*job); jobReady=true;
                }
                jobCv.notify_one();
            }
            updatePresentation();
        }
        {
            std::lock_guard<std::mutex> lock(resultMutex);
            if(resultReady) {
                std::lock_guard<std::mutex> stateLock(stateMutex);
                session.accept(std::move(latestResult));
                updatePresentation();
                resultReady=false;
            }
        }
        overlay.tick();
        MsgWaitForMultipleObjectsEx(1,&wake,50,QS_ALLINPUT,MWMO_INPUTAVAILABLE);
    }
    {
        std::lock_guard<std::mutex> lock(jobMutex);
        stopping=true;
    }
    jobCv.notify_all();
    solver.join();
    server.stop();
    overlay.destroy();
    CloseHandle(wake);
    return 0;
}

_Use_decl_annotations_ int WINAPI wWinMain(HINSTANCE,HINSTANCE,PWSTR,int) { return main(); }
