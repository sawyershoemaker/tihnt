#include "overlay.hpp"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <string>
#include <iostream>
#include <climits>

#ifdef _WIN32
#include <dwmapi.h>
#pragma comment(lib, "Dwmapi.lib")
#endif

OverlayWindow::OverlayWindow() {}
OverlayWindow::~OverlayWindow(){ destroy(); }

bool OverlayWindow::create(){
#ifdef _WIN32
    WNDCLASSW wc{}; wc.lpfnWndProc = &OverlayWindow::WndProcThunk; wc.hInstance = GetModuleHandle(nullptr); wc.lpszClassName = L"MinesOverlayWindow";
    static bool reg=false; if(!reg){ RegisterClassW(&wc); reg=true; }
    hwnd_ = CreateWindowExW(WS_EX_LAYERED|WS_EX_TOOLWINDOW|WS_EX_TOPMOST|WS_EX_TRANSPARENT,
                            wc.lpszClassName, L"Mines Overlay",
                            WS_POPUP,
                            0,0, 1,1,
                            nullptr,nullptr,wc.hInstance,this);
    if(!hwnd_) {
        std::cerr << "OverlayWindow::create: CreateWindowExW failed, GetLastError=" << GetLastError() << std::endl;
        return false;
    }
    // streamproof the overlay from most capture APIs when supported
#if defined(WDA_EXCLUDEFROMCAPTURE)
    excluded_from_capture_ = false;
    BOOL okAffinity = SetWindowDisplayAffinity(hwnd_, WDA_EXCLUDEFROMCAPTURE);
    if(!okAffinity){
        std::cerr << "OverlayWindow::create: SetWindowDisplayAffinity failed, err=" << GetLastError() << std::endl;
    }
    // reflect actual state using GetWindowDisplayAffinity when available
    refresh_exclusion_state();
#endif
#ifdef DWMWA_EXCLUDED_FROM_PEEK
    // also exclude from aero peek/thumbnail where available
    BOOL exclude = TRUE;
    DwmSetWindowAttribute(hwnd_, DWMWA_EXCLUDED_FROM_PEEK, &exclude, sizeof(exclude));
#endif
    ShowWindow(hwnd_, SW_SHOWNA);
    visible_ = true;
    // ensure always-on-top above most windows
    SetWindowPos(hwnd_, HWND_TOPMOST, 0,0,0,0, SWP_NOMOVE|SWP_NOSIZE|SWP_NOACTIVATE);
    return true;
#else
    return false;
#endif
}

void OverlayWindow::destroy(){
#ifdef _WIN32
    if(memdc_){ DeleteDC(memdc_); memdc_=nullptr; }
    if(dib_){ DeleteObject(dib_); dib_=nullptr; }
    bits_ = nullptr; stride_ = 0; surf_w_ = surf_h_ = 0;
    if(hwnd_){ DestroyWindow(hwnd_); hwnd_=nullptr; }
#endif
}

void OverlayWindow::update(const std::vector<solve::Mark>& marks, const OverlayGeometry& geom, int minesTotal){
    const double kGeomEps = 1e-3;
    auto diff = [&](double a, double b){ return std::abs(a - b) > kGeomEps; };
    bool marksChanged = marks != marks_;
    bool minesChanged = (minesTotal != mines_total_);
    bool sizeChanged = (geom.board_w != geom_.board_w) || (geom.board_h != geom_.board_h)
        || diff(geom.rect_w, geom_.rect_w) || diff(geom.rect_h, geom_.rect_h)
        || diff(geom.dpr, geom_.dpr);
    bool positionChanged = diff(geom.rect_l, geom_.rect_l) || diff(geom.rect_t, geom_.rect_t);

    if(!marksChanged && !minesChanged && !sizeChanged && !positionChanged){
        return;
    }

    marks_ = marks;
    geom_ = geom;
    mines_total_ = minesTotal;

    if(!marksChanged && !minesChanged && !sizeChanged && positionChanged){
        redraw(false);
    } else {
        redraw(true);
    }
}

void OverlayWindow::tick(){
    
}

#ifdef _WIN32
LRESULT CALLBACK OverlayWindow::WndProcThunk(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp){
    OverlayWindow* self = nullptr;
    if(msg==WM_NCCREATE){
        CREATESTRUCTW* cs = reinterpret_cast<CREATESTRUCTW*>(lp);
        self = reinterpret_cast<OverlayWindow*>(cs->lpCreateParams);
        SetWindowLongPtrW(hwnd, GWLP_USERDATA, reinterpret_cast<LONG_PTR>(self));
    } else {
        self = reinterpret_cast<OverlayWindow*>(GetWindowLongPtrW(hwnd, GWLP_USERDATA));
    }
    if(self) return self->WndProc(hwnd, msg, wp, lp);
    return DefWindowProcW(hwnd, msg, wp, lp);
}

HWND OverlayWindow::findRenderHost(){
    // anchoring (should be fixed if u ctrl+q bind)
    static HWND lastTop = nullptr;
    static HWND lastHost = nullptr;

    HWND fg = GetForegroundWindow();
    if(!fg || fg==hwnd_) return lastHost; // keep previous if available

    if(fg == lastTop && lastHost && IsWindow(lastHost)){
        return lastHost;
    }

    auto findChildByClass = [](HWND root, const wchar_t* cls)->HWND{
        for(HWND child = FindWindowExW(root, nullptr, nullptr, nullptr); child; child = FindWindowExW(root, child, nullptr, nullptr)){
            wchar_t name[256]; name[0]=0; GetClassNameW(child, name, 255);
            if(lstrcmpiW(name, cls)==0) return child;
            HWND deeper = FindWindowExW(child, nullptr, cls, nullptr);
            if(deeper) return deeper;
        }
        return nullptr;
    };

    HWND host = findChildByClass(fg, L"Chrome_RenderWidgetHostHWND");
    if(!host) host = fg;

    // if target PID set, attempt to locate a window with that PID
    if(target_pid_ != 0){
        DWORD pid=0; GetWindowThreadProcessId(host, &pid);
        if(pid != target_pid_){
            HWND best = nullptr;
            for(HWND w = GetTopWindow(nullptr); w; w = GetNextWindow(w, GW_HWNDNEXT)){
                DWORD p=0; GetWindowThreadProcessId(w, &p);
                if(p == target_pid_){ best = w; break; }
            }
            if(best){
                HWND h2 = findChildByClass(best, L"Chrome_RenderWidgetHostHWND");
                host = h2 ? h2 : best;
            }
        }
    }

    lastTop = fg;
    lastHost = host;
    return host;
}
void OverlayWindow::set_target_pid(uint32_t pid){
    target_pid_ = pid;
}

bool OverlayWindow::ensureSurface(int w, int h){
    if(w<=0 || h<=0) return false;
    if(w==surf_w_ && h==surf_h_ && dib_ && memdc_) return true;
    if(memdc_){ DeleteDC(memdc_); memdc_=nullptr; }
    if(dib_){ DeleteObject(dib_); dib_=nullptr; }
    HDC screen = GetDC(nullptr);
    // update layered window requires a top-down 32bpp DIB for per-pixel alpha (holy headache)
    BITMAPINFO bi{}; bi.bmiHeader.biSize=sizeof(BITMAPINFOHEADER); bi.bmiHeader.biWidth=w; bi.bmiHeader.biHeight=-h; bi.bmiHeader.biPlanes=1; bi.bmiHeader.biBitCount=32; bi.bmiHeader.biCompression=BI_RGB;
    void* bits=nullptr;
    dib_ = CreateDIBSection(screen, &bi, DIB_RGB_COLORS, &bits, nullptr, 0);
    memdc_ = CreateCompatibleDC(screen);
    HGDIOBJ prev = SelectObject(memdc_, dib_);
    if(prev==NULL || prev==HGDI_ERROR){
        std::cerr << "OverlayWindow::ensureSurface: SelectObject failed, GetLastError=" << GetLastError() << std::endl;
    }
    ReleaseDC(nullptr, screen);
    surf_w_ = w; surf_h_ = h;
    bits_ = reinterpret_cast<unsigned char*>(bits);
    stride_ = w * 4;
    return dib_ && memdc_ && bits_;
}

void OverlayWindow::hide(){
    if(hwnd_) ShowWindow(hwnd_, SW_HIDE);
}

void OverlayWindow::show_noactivate(){
    if(hwnd_) ShowWindow(hwnd_, SW_SHOWNA);
}

void OverlayWindow::redraw(bool full){
    if(!hwnd_) return;
    const double dpr = (geom_.dpr > 0.0 ? geom_.dpr : 1.0);
    const double scalePx = dpr;
    const double cssLeft = geom_.rect_l; // already relative to visual viewport origin
    const double cssTop  = geom_.rect_t;
    const double cssW = geom_.rect_w;
    const double cssH = geom_.rect_h;
    int dstW = (int)std::round(cssW * scalePx);
    int dstH = (int)std::round(cssH * scalePx);

    // find screen origin by converting client (0,0) of the render host
    HWND host = findRenderHost();
    if(!host){ /* we'll fallback later */ }
    // convert client origin to screen, then add device-pixel offset from CSS via DPR
    POINT origin{0,0}; int x=50, y=50;
    if(host){
        ClientToScreen(host, &origin);
        x = origin.x + (int)std::round(cssLeft * scalePx);
        y = origin.y + (int)std::round(cssTop  * scalePx);
    } else {
        // fallback position if we can't find a host window
        x = 50; y = 50;
    }

    if(!ensureSurface(dstW>0?dstW:320, dstH>0?dstH:120)){ hide(); return; }
    if(dstW<=0 || dstH<=0){ dstW = surf_w_; dstH = surf_h_; }
    // clear to fully transparent
    BLENDFUNCTION bf{}; bf.BlendOp=AC_SRC_OVER; bf.SourceConstantAlpha=255; bf.AlphaFormat=AC_SRC_ALPHA;
    // pixel writers available to subsequent draw helpers   
    auto putPixel = [&](int px, int py, uint32_t rgba){
        if(px<0||py<0||px>=dstW||py>=dstH) return;
        uint8_t a = (uint8_t)((rgba>>24)&0xFF);
        uint8_t r = (uint8_t)((rgba>>16)&0xFF);
        uint8_t g = (uint8_t)((rgba>>8)&0xFF);
        uint8_t b = (uint8_t)(rgba & 0xFF);
        uint8_t rp = (uint8_t)((r * a) / 255);
        uint8_t gp = (uint8_t)((g * a) / 255);
        uint8_t bp = (uint8_t)((b * a) / 255);
        uint32_t bgra = ((uint32_t)a<<24) | ((uint32_t)rp<<16) | ((uint32_t)gp<<8) | (uint32_t)bp;
        uint8_t* p = bits_ + py*stride_ + px*4;
        *reinterpret_cast<uint32_t*>(p) = bgra;
    };
    auto toPremulBGRA = [&](uint32_t rgba)->uint32_t{
        uint8_t a = (uint8_t)((rgba>>24)&0xFF);
        uint8_t r = (uint8_t)((rgba>>16)&0xFF);
        uint8_t g = (uint8_t)((rgba>>8)&0xFF);
        uint8_t b = (uint8_t)(rgba & 0xFF);
        uint8_t rp = (uint8_t)((r * a) / 255);
        uint8_t gp = (uint8_t)((g * a) / 255);
        uint8_t bp = (uint8_t)((b * a) / 255);
        return ((uint32_t)a<<24) | ((uint32_t)rp<<16) | ((uint32_t)gp<<8) | (uint32_t)bp;
    };
    auto putPixelV = [&](int px, int py, uint32_t bgra){
        if(px<0||py<0||px>=dstW||py>=dstH) return;
        uint8_t* p = bits_ + py*stride_ + px*4;
        *reinterpret_cast<uint32_t*>(p) = bgra;
    };

    if(full){
        std::memset(bits_, 0, (size_t)(stride_ * dstH));

        // draw overlay marks: green boxes for Safe, red X for Mine, orange box for Guess
        const int w = geom_.board_w;
        const int h = geom_.board_h;
        if(w>0 && h>0 && (int)marks_.size()==w*h){
        // precompute per-cell pixel edges to avoid repeated rounding in the hot loop
        std::vector<int> x0(w+1), y0(h+1);
        for(int xx=0; xx<=w; ++xx){ x0[xx] = (int)std::round(((double)dstW * (double)xx) / (double)w); }
        for(int yy=0; yy<=h; ++yy){ y0[yy] = (int)std::round(((double)dstH * (double)yy) / (double)h); }

        // when safety mode is enabled, block clicks for entire cells that are NOT Safe/Guess
        // by painting a minimally visible (alpha=1) mask across those cells.
        if(safety_mode_){
            const uint32_t kMaskBGRA = 0x01000000; // alpha=1, RGB=0 (premultiplied)
            for(int yy=0; yy<h; ++yy){
                const int cy0 = y0[yy];
                const int cy1 = y0[yy+1] - 1;
                for(int xx=0; xx<w; ++xx){
                    const int cx0 = x0[xx];
                    const int cx1 = x0[xx+1] - 1;
                    solve::Mark m = marks_[yy*w + xx];
                    bool allow = (m==solve::Mark::Safe || m==solve::Mark::Guess);
                    if(!allow){
                        for(int py=cy0; py<=cy1; ++py){
                            uint8_t* row = bits_ + py*stride_ + cx0*4;
                            for(int px=cx0; px<=cx1; ++px){
                                *reinterpret_cast<uint32_t*>(row) = kMaskBGRA;
                                row += 4;
                            }
                        }
                    }
                }
            }
        }

        // fast check: if no marks to draw, skip per-cell drawing
        bool anyMarks=false; for(const auto m : marks_){ if(m!=solve::Mark::None){ anyMarks=true; break; } }
        if(anyMarks){
        
        // cell size in target surface pixels
        double cellW = (double)dstW / (double)w;
        double cellH = (double)dstH / (double)h;
        const int thickness = 5; // thicker strokes for visibility
        auto drawRectStroke = [&](int x0, int y0, int x1, int y1, uint32_t color){
            uint32_t colv = toPremulBGRA(color);
            for(int t=0; t<thickness; ++t){
                int yt0 = y0 + t, yt1 = y1 - t, xt0 = x0 + t, xt1 = x1 - t;
                for(int x=xt0; x<=xt1; ++x){ putPixelV(x, yt0, colv); putPixelV(x, yt1, colv); }
                for(int y=yt0; y<=yt1; ++y){ putPixelV(xt0, y, colv); putPixelV(xt1, y, colv); }
            }
        };
        auto fillRect = [&](int x0, int y0, int x1, int y1, uint32_t color){
            uint32_t colv = toPremulBGRA(color);
            if(x1<x0) x1=x0; if(y1<y0) y1=y0;
            for(int y=y0; y<=y1; ++y){
                uint8_t* row = bits_ + y*stride_ + x0*4;
                for(int x=x0; x<=x1; ++x){ *reinterpret_cast<uint32_t*>(row) = colv; row += 4; }
            }
        };
        auto drawGlyphFScaled = [&](int cellX0, int cellY0, int cellX1, int cellY1, uint32_t color){
            uint32_t colv = toPremulBGRA(color);
            static const unsigned char glyphF[7] = {
                0x1F,
                0x10,
                0x1E,
                0x10,
                0x10,
                0x10,
                0x10
            };
            int cellWpx = std::max(1, cellX1 - cellX0 + 1);
            int cellHpx = std::max(1, cellY1 - cellY0 + 1);
            int scale = std::max(1, std::min(cellWpx / 6, cellHpx / 8));
            int gw = 5 * scale, gh = 7 * scale;
            int px0 = cellX0 + (cellWpx - gw) / 2;
            int py0 = cellY0 + (cellHpx - gh) / 2;
            for(int r=0; r<7; ++r){
                unsigned char row = glyphF[r];
                for(int c=0; c<5; ++c){
                    if(row & (1 << (4-c))){
                        int rx = px0 + c*scale;
                        int ry = py0 + r*scale;
                        for(int dy=0; dy<scale; ++dy){
                            for(int dx=0; dx<scale; ++dx){ putPixelV(rx+dx, ry+dy, colv); }
                        }
                    }
                }
            }
        };
        auto drawX = [&](int x0, int y0, int x1, int y1, uint32_t color){
            uint32_t colv = toPremulBGRA(color);
            // slightly inset the X endpoints so it isn't too long
            int inset = std::max(2, std::min((x1-x0), (y1-y0)) / 6);
            int lx0 = x0 + inset, ly0 = y0 + inset, lx1 = x1 - inset, ly1 = y1 - inset;
            int w=lx1-lx0, h=ly1-ly0; int n=std::max(w,h);
            for(int i=0;i<=n;++i){
                int px=lx0 + i*w/n; int py=ly0 + i*h/n;
                int px2=lx1 - i*w/n; int py2=ly0 + i*h/n;
                // thicker stroke
                for(int dx=-(thickness/2); dx<= (thickness/2); ++dx){ putPixelV(px+dx, py, colv); putPixelV(px2+dx, py2, colv); }
            }
        };

        for(int yy=0; yy<h; ++yy){
            const int cy0 = y0[yy];
            const int cy1 = y0[yy+1] - 1;
            for(int xx=0; xx<w; ++xx){
                solve::Mark m = marks_[yy*w + xx];
                if(m==solve::Mark::None) continue; // only overlay the specified categories
                int cx0 = x0[xx];
                int cx1 = x0[xx+1] - 1;
                if(cx1<cx0) cx1=cx0;
                if(m==solve::Mark::Safe){
                    drawRectStroke(cx0, cy0, cx1, cy1, 0xAA00FF00); // A,R,G,B packed as ARGB
                } else if(m==solve::Mark::Mine){
                    // red X for generic mine suggestions
                    drawX(cx0, cy0, cx1, cy1, 0xCCFF0000);
                } else if(m==solve::Mark::Guess){
                    drawRectStroke(cx0, cy0, cx1, cy1, 0xA0FFA500); // orange
                } else if(m==solve::Mark::Chord){
                    // blue 'F' on the numbered cell to click for chord (fully opaque)
                    drawGlyphFScaled(cx0, cy0, cx1, cy1, 0xFF1E90FF);
                } else if(m==solve::Mark::FlagForChord){
                    // blue filled square on mines that need flagging for the shown chord
                    fillRect(cx0, cy0, cx1, cy1, 0x661E90FF);
                    drawRectStroke(cx0, cy0, cx1, cy1, 0xFF1E90FF);
                }
            }
        }
        }
    }
    // draw debug/mines total text in top-left of overlay surface
    // simple 5x7 bitmap font for digits and a few chars
    auto drawChar = [&](int x0, int y0, char ch, uint32_t color){
        uint32_t colv = toPremulBGRA(color);
        static const unsigned char font5x7[10][7] = {
            {0x1E,0x11,0x13,0x15,0x19,0x11,0x1E}, // 0
            {0x04,0x0C,0x14,0x04,0x04,0x04,0x1F}, // 1
            {0x1E,0x01,0x01,0x1E,0x10,0x10,0x1F}, // 2
            {0x1E,0x01,0x01,0x0E,0x01,0x01,0x1E}, // 3
            {0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}, // 4
            {0x1F,0x10,0x10,0x1E,0x01,0x01,0x1E}, // 5
            {0x0E,0x10,0x10,0x1E,0x11,0x11,0x1E}, // 6
            {0x1F,0x01,0x02,0x04,0x08,0x08,0x08}, // 7
            {0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}, // 8
            {0x1E,0x11,0x11,0x1E,0x01,0x01,0x1E}  // 9
        };
        auto plot = [&](int px, int py){ putPixelV(px, py, colv); };
        if(ch>='0' && ch<='9'){
            int d = ch - '0';
            for(int r=0;r<7;++r){ unsigned char row = font5x7[d][r]; for(int c=0;c<5;++c){ if(row & (1<<(4-c))){ plot(x0+c, y0+r); } } }
        } else if(ch=='[' || ch==']' || ch=='='){
            // very simple brackets/equals
            if(ch=='['){ for(int r=0;r<7;++r){ plot(x0, y0+r); } for(int c=0;c<3;++c){ plot(x0+c, y0); plot(x0+c, y0+6);} }
            if(ch==']'){ for(int r=0;r<7;++r){ plot(x0+3, y0+r); } for(int c=0;c<3;++c){ plot(x0+c+1, y0); plot(x0+c+1, y0+6);} }
            if(ch=='='){ for(int c=0;c<5;++c){ plot(x0+c, y0+2); plot(x0+c, y0+4);} }
        } else if(ch=='M'){
            // crude M
            for(int r=0;r<7;++r){ plot(x0, y0+r); plot(x0+4, y0+r);} plot(x0+1,y0+1); plot(x0+2,y0+2); plot(x0+3,y0+1);
        } else if(ch=='V' || ch=='I' || ch=='S' || ch=='B' || ch=='L' || ch=='E' || ch=='A' || ch=='F'){
            // 5x7 uppercase letters needed for UI labels
            unsigned char glyph[7] = {0,0,0,0,0,0,0};
            if(ch=='V'){ unsigned char g[7] = {0x11,0x11,0x11,0x11,0x11,0x0A,0x04}; std::memcpy(glyph,g,7); }
            if(ch=='I'){ unsigned char g[7] = {0x0E,0x04,0x04,0x04,0x04,0x04,0x0E}; std::memcpy(glyph,g,7); }
            if(ch=='S'){ unsigned char g[7] = {0x0E,0x10,0x10,0x0E,0x01,0x01,0x0E}; std::memcpy(glyph,g,7); }
            if(ch=='B'){ unsigned char g[7] = {0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}; std::memcpy(glyph,g,7); }
            if(ch=='L'){ unsigned char g[7] = {0x10,0x10,0x10,0x10,0x10,0x10,0x1F}; std::memcpy(glyph,g,7); }
            if(ch=='E'){ unsigned char g[7] = {0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}; std::memcpy(glyph,g,7); }
            if(ch=='A'){ unsigned char g[7] = {0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}; std::memcpy(glyph,g,7); }
            if(ch=='F'){ unsigned char g[7] = {0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}; std::memcpy(glyph,g,7); }
            for(int r=0;r<7;++r){ unsigned char row = glyph[r]; for(int c=0;c<5;++c){ if(row & (1<<(4-c))){ plot(x0+c, y0+r); } } }
        }
    };
    auto drawText = [&](int x, int y, const std::string& s, uint32_t color){ int pen=x; for(char ch : s){ drawChar(pen,y,ch,color); pen += 6; } };
    // always draw M= text even if board dims are zero-sized surface; guard for tiny overlays
        if(dstW>=20 && dstH>=10){
            std::string label = "M=" + std::to_string(std::max(0, mines_total_));
            drawText(2, 2, label, 0xFFFFFFFF);
            if(!excluded_from_capture_){
                drawText(2, 12, std::string("VISIBLE"), 0xFFFFFFFF);
            }
            if(safety_mode_){
                drawText(2, 22, std::string("SAFE"), 0xFF00FF00);
            }
        }
    }

    // present using UpdateLayeredWindow (per-pixel alpha)
    HDC screen = GetDC(nullptr);
    SIZE siz{dstW, dstH}; POINT ptSrc{0,0}; POINT ptDst{x,y};
    BOOL updOk = UpdateLayeredWindow(hwnd_, screen, &ptDst, &siz, memdc_, &ptSrc, 0, &bf, ULW_ALPHA);
    if(!updOk){
        DWORD err = GetLastError();
        std::cerr << "OverlayWindow::redraw: UpdateLayeredWindow failed, err=" << err
                  << " w=" << dstW << " h=" << dstH << " x=" << x << " y=" << y << std::endl;
        ptDst.x = 50; ptDst.y = 50; siz.cx = 320; siz.cy = 120;
        ensureSurface(siz.cx, siz.cy);
        std::memset(bits_, 0, (size_t)(stride_ * siz.cy));
        for(int yy=0; yy<siz.cy; ++yy){ for(int xx=0; xx<siz.cx; ++xx){ putPixel(xx,yy, 0x7F000000); } }
        updOk = UpdateLayeredWindow(hwnd_, screen, &ptDst, &siz, memdc_, &ptSrc, 0, &bf, ULW_ALPHA);
        if(!updOk){
            std::cerr << "OverlayWindow::redraw: Fallback UpdateLayeredWindow still failing, err=" << GetLastError() << std::endl;
        }
    }
    ReleaseDC(nullptr, screen);
    ShowWindow(hwnd_, SW_SHOWNA);
}

LRESULT OverlayWindow::WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp){
    switch(msg){
    case WM_MOUSEWHEEL:
    case WM_MOUSEHWHEEL: {
        // while safety mode is active, forward wheel events to the underlying host (not working)
        if(safety_mode_){
            HWND host = findRenderHost();
            if(host){ PostMessage(host, msg, wp, lp); }
            return 0; // eat it so overlay doesn't interfere
        }
        break;
    }
    case WM_NCHITTEST: {
        // safety-mode selective pass-through: only allow clicks on Safe/Guess cells
        // otherwise, return HTCLIENT to block (prevent click-through)
        if(safety_mode_ && !marks_.empty() && geom_.board_w>0 && geom_.board_h>0){
            // determine mouse position in overlay client coordinates (convert from screen)
            POINTS pts = MAKEPOINTS(lp);
            POINT p{ pts.x, pts.y };
            ScreenToClient(hwnd_, &p);
            // map to cell indices using current surface size
            RECT rc; GetClientRect(hwnd_, &rc);
            int dstW = rc.right - rc.left;
            int dstH = rc.bottom - rc.top;
            if(dstW>0 && dstH>0){
                int x = p.x; int y = p.y;
                const int w = geom_.board_w;
                const int h = geom_.board_h;
                // build the same rounded edge arrays used for drawing for consistent mapping
                std::vector<int> xEdge(w+1), yEdge(h+1);
                for(int xx=0; xx<=w; ++xx){ xEdge[xx] = (int)std::round(((double)dstW * (double)xx) / (double)w); }
                for(int yy=0; yy<=h; ++yy){ yEdge[yy] = (int)std::round(((double)dstH * (double)yy) / (double)h); }
                auto locate = [&](int px, int py, int& cx, int& cy){
                    cx = -1; cy = -1;
                    // find column
                    for(int xx=0; xx<w; ++xx){ if(px >= xEdge[xx] && px <= xEdge[xx+1]-1){ cx = xx; break; } }
                    for(int yy=0; yy<h; ++yy){ if(py >= yEdge[yy] && py <= yEdge[yy+1]-1){ cy = yy; break; } }
                };
                int cx=-1, cy=-1; locate(x,y,cx,cy);
                auto isSafeGuessAt = [&](int ix, int iy)->bool{
                    if(ix<0||iy<0||ix>=w||iy>=h) return false;
                    size_t idx = (size_t)iy * (size_t)w + (size_t)ix;
                    if(idx >= marks_.size()) return false;
                    auto m = marks_[idx];
                    return (m==solve::Mark::Safe || m==solve::Mark::Guess);
                };
                // if exact cell is safe/guess, pass through
                if(cx>=0 && cy>=0 && isSafeGuessAt(cx,cy)){
                    return HTTRANSPARENT;
                }
                // otherwise, if the point lies inside any neighboring safe/guess cell rectangle, pass through
                for(int dy=-1; dy<=1; ++dy){
                    for(int dx=-1; dx<=1; ++dx){
                        int nx = (cx<0? -1 : cx+dx);
                        int ny = (cy<0? -1 : cy+dy);
                        if(!isSafeGuessAt(nx,ny)) continue;
                        int rx0 = xEdge[nx]; int rx1 = xEdge[nx+1]-1;
                        int ry0 = yEdge[ny]; int ry1 = yEdge[ny+1]-1;
                        if(x>=rx0 && x<=rx1 && y>=ry0 && y<=ry1){
                            return HTTRANSPARENT;
                        }
                    }
                }
                // all other areas: block
                return HTCLIENT;
            }
        }
        // default: transparent unless explicitly blocking
        return HTTRANSPARENT;
    }
    case WM_ERASEBKGND: return 1;
    default: break;
    }
    return DefWindowProcW(hwnd, msg, wp, lp);
}
#endif

void OverlayWindow::set_excluded_from_capture(bool exclude){
#ifdef _WIN32
#if defined(WDA_EXCLUDEFROMCAPTURE)
    if(hwnd_){
        BOOL ok = SetWindowDisplayAffinity(hwnd_, exclude ? WDA_EXCLUDEFROMCAPTURE : WDA_NONE);
        if(!ok){
            std::cerr << "OverlayWindow::set_excluded_from_capture: SetWindowDisplayAffinity failed, err=" << GetLastError() << std::endl;
        }
        // reflect actual state by querying current affinity when possible
        refresh_exclusion_state();
    }
#else
    (void)exclude;
#endif
#else
    (void)exclude;
#endif
}

bool OverlayWindow::is_excluded_from_capture() const{
    return excluded_from_capture_;
}

void OverlayWindow::set_visible(bool visible){
#ifdef _WIN32
    if(!hwnd_) { visible_ = visible; return; }
    if(visible){ show_noactivate(); }
    else { hide(); }
    visible_ = visible;
#else
    (void)visible;
#endif
}

bool OverlayWindow::is_visible() const{
    return visible_;
}

void OverlayWindow::set_safety_mode(bool enabled){
    safety_mode_ = enabled;
#ifdef _WIN32
    if(hwnd_){
        LONG_PTR ex = GetWindowLongPtr(hwnd_, GWL_EXSTYLE);
        if(enabled){ ex &= ~WS_EX_TRANSPARENT; } else { ex |= WS_EX_TRANSPARENT; }
        SetWindowLongPtr(hwnd_, GWL_EXSTYLE, ex);
        SetWindowPos(hwnd_, HWND_TOPMOST, 0,0,0,0, SWP_NOMOVE|SWP_NOSIZE|SWP_NOACTIVATE|SWP_FRAMECHANGED);
        InvalidateRect(hwnd_, nullptr, FALSE);
    }
#endif
}

bool OverlayWindow::is_safety_mode() const{
    return safety_mode_;
}


#ifdef _WIN32
bool OverlayWindow::refresh_exclusion_state(){
#if defined(WDA_EXCLUDEFROMCAPTURE)
    if(!hwnd_) { excluded_from_capture_ = false; return false; }
    DWORD affinity = 0;
    BOOL ok = GetWindowDisplayAffinity(hwnd_, &affinity);
    if(!ok){
        std::cerr << "OverlayWindow::refresh_exclusion_state: GetWindowDisplayAffinity failed, err=" << GetLastError() << std::endl;
        return false;
    }
    excluded_from_capture_ = (affinity == WDA_EXCLUDEFROMCAPTURE);
    return true;
#else
    excluded_from_capture_ = false;
    return false;
#endif
}
#endif