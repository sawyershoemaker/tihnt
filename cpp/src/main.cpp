#include <atomic>
#include <condition_variable>
#include <cstdint>
#include <iostream>
#include <mutex>
#include <thread>
#include <utility>

#include "ws_server.hpp"
#include "board.hpp"
#include "proto.hpp"
#include "solver.hpp"
#include "overlay.hpp"

namespace {
struct SolverJob {
    game::Board board;
    int minesTotal=-1;
    bool enableChords=true;
    uint64_t revision=0;
};
struct SolverResult {
    std::vector<solve::Mark> marks;
    uint64_t revision=0;
};

OverlayGeometry geometry_for(const proto::GeometryMsg& message, int w, int h) {
    OverlayGeometry g{};
    g.board_w=w; g.board_h=h;
    g.rect_l=message.rect_l; g.rect_t=message.rect_t;
    g.rect_w=message.rect_w; g.rect_h=message.rect_h;
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

    OverlayWindow overlay;
    if(!overlay.create()) return 1;
    HANDLE wake=CreateEventW(nullptr,FALSE,FALSE,nullptr);
    if(!wake) return 1;

    std::mutex stateMutex;
    game::Board board;
    OverlayGeometry geometry{};
    int minesTotal=-1;
    std::atomic<uint64_t> revision{0};
    std::atomic<bool> dirty{false}, geometryDirty{false}, bindingDirty{false}, stopping{false};
    uint32_t targetPid=0;
    bool enableChords=true;

    net::WebSocketServer server;
    server.set_on_connection([&](bool) {
        {
            std::lock_guard<std::mutex> lock(stateMutex);
            board.resize(0,0); geometry={}; minesTotal=-1; targetPid=0;
            ++revision;
            dirty=true;
            bindingDirty=true;
        }
        SetEvent(wake);
    });
    server.set_on_message([&](const net::WsMessage& msg) {
        proto::ParsedMessage parsed;
        if(!proto::parse_message(msg.text,parsed)) return;
        {
            std::lock_guard<std::mutex> lock(stateMutex);
            bool changed=false;
            if(parsed.type==proto::MsgType::Full) {
                const auto& f=parsed.full;
                changed=board.width()!=f.w || board.height()!=f.h || board.data()!=f.cells || minesTotal!=f.mines_total;
                if(changed) board.apply_full(std::move(parsed.full.cells),f.w,f.h);
                minesTotal=f.mines_total;
                geometry=geometry_for(f,f.w,f.h);
                geometryDirty=true;
            } else if(parsed.type==proto::MsgType::Delta) {
                if(board.width()==0) return; // A delta never establishes a board.
                const auto& d=parsed.delta;
                for(const auto& u:d.updates) if(u.x>=board.width() || u.y>=board.height()) return;
                for(const auto& u:d.updates) {
                    if(board.at(u.x,u.y)!=u.state) { board.set(u.x,u.y,u.state); changed=true; }
                }
                if(d.has_mines_total && minesTotal!=d.mines_total) { minesTotal=d.mines_total; changed=true; }
                if(d.has_geometry) {
                    geometry=geometry_for(d,board.width(),board.height());
                    geometryDirty=true;
                }
            } else if(parsed.type==proto::MsgType::Bind) {
                targetPid=static_cast<uint32_t>(parsed.bind.pid);
                bindingDirty=true;
            }
            if(changed) { ++revision; dirty=true; }
        }
        SetEvent(wake);
    });
    if(!server.start(8765)) {
        std::cerr<<"Failed to start WebSocket server on 127.0.0.1:8765\n";
        CloseHandle(wake); return 1;
    }

    constexpr int Exit=1, Chords=2, Capture=3, Safety=4;
    RegisterHotKey(nullptr,Exit,MOD_CONTROL|MOD_ALT|MOD_NOREPEAT,'X');
    RegisterHotKey(nullptr,Chords,MOD_CONTROL|MOD_NOREPEAT,'P');
    RegisterHotKey(nullptr,Capture,MOD_CONTROL|MOD_ALT|MOD_SHIFT|MOD_NOREPEAT,'W');
    RegisterHotKey(nullptr,Safety,MOD_CONTROL|MOD_ALT|MOD_NOREPEAT,'S');

    std::mutex jobMutex,resultMutex;
    std::condition_variable jobCv;
    SolverJob pendingJob;
    SolverResult latestResult;
    bool jobReady=false,resultReady=false;
    std::thread solver([&] {
        for(;;) {
            SolverJob job;
            {
                std::unique_lock<std::mutex> lock(jobMutex);
                jobCv.wait(lock,[&] { return stopping || jobReady; });
                if(stopping) break;
                job=std::move(pendingJob); jobReady=false;
            }
            auto cancelled=[&] { return stopping.load() || revision.load()!=job.revision; };
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

    std::vector<solve::Mark> lastMarks;
    bool running=true;
    while(running) {
        MSG msg;
        while(PeekMessageW(&msg,nullptr,0,0,PM_REMOVE)) {
            if(msg.message==WM_QUIT) { running=false; break; }
            if(msg.message==WM_HOTKEY) {
                if(msg.wParam==Exit) { running=false; break; }
                if(msg.wParam==Chords) {
                    std::lock_guard<std::mutex> lock(stateMutex);
                    enableChords=!enableChords; ++revision; dirty=true;
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
            if(bindingDirty.exchange(false)) overlay.set_target_pid(targetPid);
            if(dirty.exchange(false)) {
                SolverJob job{board,minesTotal,enableChords,revision.load()};
                {
                    std::lock_guard<std::mutex> jobLock(jobMutex);
                    pendingJob=std::move(job); jobReady=true;
                }
                jobCv.notify_one();
                // Previous marks may refer to another game, even at the same size.
                lastMarks.assign(board.data().size(),solve::Mark::None);
                geometryDirty=true;
            }
            if(geometryDirty.exchange(false)) overlay.update(lastMarks,geometry,minesTotal);
        }
        {
            std::lock_guard<std::mutex> lock(resultMutex);
            if(resultReady) {
                std::lock_guard<std::mutex> stateLock(stateMutex);
                if(latestResult.revision==revision.load()) {
                    lastMarks=std::move(latestResult.marks);
                    // Use current geometry: scrolling does not invalidate a solve.
                    overlay.update(lastMarks,geometry,minesTotal);
                }
                resultReady=false;
            }
        }
        overlay.tick();
        MsgWaitForMultipleObjects(1,&wake,FALSE,50,QS_ALLINPUT);
    }
    for(int id : {Exit,Chords,Capture,Safety}) UnregisterHotKey(nullptr,id);
    stopping=true;
    jobCv.notify_all();
    solver.join();
    server.stop();
    overlay.destroy();
    CloseHandle(wake);
    return 0;
}

int WINAPI wWinMain(HINSTANCE,HINSTANCE,PWSTR,int) { return main(); }
