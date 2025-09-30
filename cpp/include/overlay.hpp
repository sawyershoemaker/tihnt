#pragma once

#include <vector>
#include <cstdint>

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
    bool click_nearest_safe_guess();
#endif

private:
#ifdef _WIN32
    static LRESULT CALLBACK WndProcThunk(HWND, UINT, WPARAM, LPARAM);
    LRESULT WndProc(HWND, UINT, WPARAM, LPARAM);
    void redraw();
    bool ensureSurface(int w, int h);
    HWND findRenderHost();
    void hide();
    void show_noactivate();

    HWND hwnd_ = nullptr;
    HBITMAP dib_ = nullptr;
    HDC memdc_ = nullptr;
    int surf_w_ = 0;
    int surf_h_ = 0;
    unsigned char* bits_ = nullptr;
    int stride_ = 0;
#endif
    std::vector<solve::Mark> marks_;
    OverlayGeometry geom_{};
    int mines_total_ = 0;
    bool excluded_from_capture_ = true;
    bool visible_ = true;
    bool safety_mode_ = false;
};


