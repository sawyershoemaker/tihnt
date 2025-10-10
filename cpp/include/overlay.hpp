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
    void set_target_pid(uint32_t pid);
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
    bool human_move_to(int targetScreenX, int targetScreenY);
    enum class ClickType : uint8_t { Left, Right, ChordBoth };
    bool host_click_at_screen(int screenX, int screenY, ClickType type);
    bool refresh_exclusion_state();

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
    bool excluded_from_capture_ = false;
    bool visible_ = true;
    bool safety_mode_ = false;
#ifdef _WIN32
#  if defined(_WIN32)
    unsigned long target_pid_ = 0;
#  else
    unsigned long target_pid_ = 0;
#  endif
#endif
};


