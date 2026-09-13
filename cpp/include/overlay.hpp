#pragma once

#include <vector>
#include <cstdint>
#include <atomic>
#include <thread>

#ifdef _WIN32
#  include <windows.h>
#endif

#include "solver.hpp"

struct OverlayGeometry {
    int board_w = 0;
    int board_h = 0;
    double rect_l = 0.0;
    double rect_t = 0.0;
    double rect_w = 0.0;
    double rect_h = 0.0;
    bool has_clip = false;
    double clip_l = 0.0;
    double clip_t = 0.0;
    double clip_w = 0.0;
    double clip_h = 0.0;
    double vv_x = 0.0;
    double vv_y = 0.0;
    double vv_scale = 1.0;
    double dpr = 1.0;
};

class OverlayWindow {
public:
    OverlayWindow();
    ~OverlayWindow();

    bool create();
    void destroy();

    void update(const std::vector<solve::Mark>& marks, const OverlayGeometry& geom, int minesTotal);

    void tick();

    void set_excluded_from_capture(bool exclude);
    bool is_excluded_from_capture() const;

    void set_visible(bool visible);
    bool is_visible() const;

    void set_safety_mode(bool enabled);
    bool is_safety_mode() const;

#ifdef _WIN32
    void set_target_pid(uint32_t pid);
#endif

private:
    friend struct OverlayTestAccess;
#ifdef _WIN32
    static LRESULT CALLBACK WndProcThunk(HWND, UINT, WPARAM, LPARAM);
    static LRESULT CALLBACK MouseHook(int, WPARAM, LPARAM);
    static OverlayWindow* hook_owner_;
    static void CALLBACK WindowEvent(HWINEVENTHOOK, DWORD, HWND, LONG, LONG, DWORD, DWORD);
    static thread_local std::vector<OverlayWindow*> event_owners_;
    LRESULT WndProc(HWND, UINT, WPARAM, LPARAM);
    void redraw(bool full = true);
    void paint_surface();
    bool ensureSurface(int w, int h);
    HWND findRenderHost();
    bool host_viewport(HWND host, POINT& origin, RECT& viewport, HRGN& region) const;
    void watch_browser(HWND browser);
    void hide();
    bool refresh_exclusion_state();
    static RECT visible_board_rect(int x, int y, int w, int h, const RECT& viewport);
    bool board_bounds(RECT& bounds) const;
    RECT board_clip(const RECT& bounds, const RECT& viewport) const;
    bool board_region(const RECT& bounds, const RECT& viewport, HRGN host, HRGN region) const;
    bool unsafe_at(POINT point, const RECT& clip) const;
    bool block_mouse(WPARAM message, POINT point);

    HWND hwnd_ = nullptr;
    HBITMAP dib_ = nullptr;
    HDC memdc_ = nullptr;
    int surf_w_ = 0;
    int surf_h_ = 0;
    unsigned char* bits_ = nullptr;
    int stride_ = 0;
    HHOOK mouse_hook_ = nullptr;
    HWINEVENTHOOK foreground_hook_ = nullptr;
    HWINEVENTHOOK location_hook_ = nullptr;
    HWINEVENTHOOK visibility_hook_ = nullptr;
    DWORD watched_pid_ = 0;
    bool refresh_pending_ = false;
    bool blocked_left_ = false;
    HWND cached_top_ = nullptr;
    HWND cached_host_ = nullptr;
    HWND (WINAPI *foreground_window_)() = &GetForegroundWindow;
    POINT last_host_origin_{};
    RECT last_host_rect_{};
    HRGN query_region_ = nullptr;
    HRGN last_host_region_ = nullptr;
    RECT visible_rect_{};
    HRGN visible_region_ = nullptr;
    bool has_clip_ = false;
#endif
    std::vector<solve::Mark> marks_;
    OverlayGeometry geom_{};
	int mines_total_ = -1;
    bool excluded_from_capture_ = false;
    bool visible_ = true;
    bool safety_mode_ = false;
    std::vector<int> cached_xEdge_;
    std::vector<int> cached_yEdge_;
    int cached_dstW_ = 0;
    int cached_dstH_ = 0;
    bool surface_painted_ = false;
    bool surface_dirty_ = true;
    std::vector<solve::Mark> painted_marks_;
    int painted_mines_total_ = -1;
    bool painted_exclusion_ = false;
    bool painted_safety_ = false;
#ifdef _WIN32
    unsigned long target_pid_ = 0;
#endif
};
