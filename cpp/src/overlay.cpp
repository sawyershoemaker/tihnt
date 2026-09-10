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

#ifdef _WIN32
OverlayWindow* OverlayWindow::hook_owner_ = nullptr;
#endif

OverlayWindow::OverlayWindow() {}
OverlayWindow::~OverlayWindow(){ destroy(); }

bool OverlayWindow::create(){
#ifdef _WIN32
    WNDCLASSW wc{}; wc.lpfnWndProc = &OverlayWindow::WndProcThunk; wc.hInstance = GetModuleHandle(nullptr); wc.lpszClassName = L"MinesOverlayWindow";
    static bool reg=false; if(!reg){ RegisterClassW(&wc); reg=true; }
    hwnd_ = CreateWindowExW(WS_EX_LAYERED|WS_EX_TOOLWINDOW|WS_EX_TOPMOST|WS_EX_TRANSPARENT|WS_EX_NOACTIVATE,
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
    visible_ = true; // Remain hidden until a valid board and browser host arrive.
    // ensure always-on-top above most windows
    SetWindowPos(hwnd_, HWND_TOPMOST, 0,0,0,0, SWP_NOMOVE|SWP_NOSIZE|SWP_NOACTIVATE);
    return true;
#else
    return false;
#endif
}

void OverlayWindow::destroy(){
#ifdef _WIN32
    if(mouse_hook_){ UnhookWindowsHookEx(mouse_hook_); mouse_hook_=nullptr; }
    if(hook_owner_==this) hook_owner_=nullptr;
    blocked_left_=false;
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
        || diff(geom.dpr, geom_.dpr) || diff(geom.vv_scale, geom_.vv_scale);
    bool positionChanged = diff(geom.rect_l, geom_.rect_l) || diff(geom.rect_t, geom_.rect_t)
        || diff(geom.vv_x, geom_.vv_x) || diff(geom.vv_y, geom_.vv_y);

    if(!marksChanged && !minesChanged && !sizeChanged && !positionChanged){
        return;
    }

    if(marksChanged) marks_ = marks;
    geom_ = geom;
    mines_total_ = minesTotal;

    if(!marksChanged && !minesChanged && !sizeChanged && positionChanged){
        redraw(false);
    } else {
        redraw(true);
    }
}

void OverlayWindow::tick(){
#ifdef _WIN32
    if(!hwnd_) return;
    HWND host=findRenderHost();
    if(!visible_ || !host || geom_.board_w<=0 || geom_.board_h<=0 || geom_.rect_w<=0 || geom_.rect_h<=0) { hide(); return; }
    POINT origin{0,0};
    if(!ClientToScreen(host,&origin)) { hide(); return; }
    if(!IsWindowVisible(hwnd_)) redraw(true);
    else if(origin.x!=last_host_origin_.x || origin.y!=last_host_origin_.y) redraw(false);
#endif
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
    HWND foreground=GetForegroundWindow();
    if(!foreground || foreground==hwnd_ || IsIconic(foreground)) return nullptr;
    wchar_t cls[128]{};
    GetClassNameW(foreground,cls,128);
    if(lstrcmpW(cls,L"Chrome_WidgetWin_1")!=0) return nullptr;
    if(foreground!=cached_top_ || !cached_host_ || !IsWindow(cached_host_) || !IsWindowVisible(cached_host_)) {
        cached_top_=foreground; cached_host_=nullptr;
        EnumChildWindows(foreground,[](HWND child,LPARAM param)->BOOL {
            wchar_t name[128]{};
            GetClassNameW(child,name,128);
            if(lstrcmpW(name,L"Chrome_RenderWidgetHostHWND")==0 && IsWindowVisible(child)) {
                *reinterpret_cast<HWND*>(param)=child; return FALSE;
            }
            return TRUE;
        },reinterpret_cast<LPARAM>(&cached_host_));
    }
    if(target_pid_ && cached_host_) {
        DWORD topPid=0,hostPid=0;
        GetWindowThreadProcessId(foreground,&topPid);
        GetWindowThreadProcessId(cached_host_,&hostPid);
        if(topPid!=target_pid_ && hostPid!=target_pid_) return nullptr;
    }
    return cached_host_;
}
void OverlayWindow::set_target_pid(uint32_t pid){
    if(target_pid_==pid) return;
    target_pid_=pid; cached_host_=cached_top_=nullptr;
    redraw(true);
}

bool OverlayWindow::ensureSurface(int w, int h){
    if(w<=0 || h<=0 || w>8192 || h>8192 || static_cast<uint64_t>(w)*h>16*1024*1024) return false;
    if(w==surf_w_ && h==surf_h_ && dib_ && memdc_ && bits_) return true;
    if(memdc_){ DeleteDC(memdc_); memdc_=nullptr; }
    if(dib_){ DeleteObject(dib_); dib_=nullptr; }
    bits_=nullptr; surf_w_=surf_h_=stride_=0;
    HDC screen=GetDC(nullptr);
    if(!screen) return false;
    BITMAPINFO bi{};
    bi.bmiHeader.biSize=sizeof(BITMAPINFOHEADER); bi.bmiHeader.biWidth=w; bi.bmiHeader.biHeight=-h;
    bi.bmiHeader.biPlanes=1; bi.bmiHeader.biBitCount=32; bi.bmiHeader.biCompression=BI_RGB;
    void* bits=nullptr;
    dib_=CreateDIBSection(screen,&bi,DIB_RGB_COLORS,&bits,nullptr,0);
    memdc_=CreateCompatibleDC(screen);
    ReleaseDC(nullptr,screen);
    HGDIOBJ selected=(dib_ && memdc_) ? SelectObject(memdc_,dib_) : nullptr;
    if(!selected || selected==HGDI_ERROR || !bits) {
        if(memdc_) { DeleteDC(memdc_); memdc_=nullptr; }
        if(dib_) { DeleteObject(dib_); dib_=nullptr; }
        return false;
    }
    surf_w_=w; surf_h_=h; stride_=w*4; bits_=static_cast<unsigned char*>(bits);
    return true;
}

void OverlayWindow::hide(){
    if(hwnd_) ShowWindow(hwnd_, SW_HIDE);
}

void OverlayWindow::show_noactivate(){
    if(hwnd_) ShowWindow(hwnd_, SW_SHOWNA);
}

void OverlayWindow::redraw(bool full){
    if(!hwnd_) return;
    HWND host=findRenderHost();
    if(!visible_ || !host || !game::valid_dimensions(geom_.board_w,geom_.board_h) ||
        geom_.board_w==0 || geom_.board_h==0) { hide(); return; }
    const double scalePx=geom_.dpr*geom_.vv_scale;
    const double width=geom_.rect_w*scalePx, height=geom_.rect_h*scalePx;
    const double left=(geom_.rect_l-geom_.vv_x)*scalePx, top=(geom_.rect_t-geom_.vv_y)*scalePx;
    if(!std::isfinite(width) || !std::isfinite(height) || !std::isfinite(left) || !std::isfinite(top) ||
        width<0.5 || height<0.5 || width>8192 || height>8192 || std::abs(left)>1000000 || std::abs(top)>1000000) { hide(); return; }
    const int dstW=static_cast<int>(std::round(width)), dstH=static_cast<int>(std::round(height));
    POINT origin{0,0};
    if(!ClientToScreen(host,&origin)) { hide(); return; }
    last_host_origin_=origin;
    const int x=origin.x+static_cast<int>(std::round(left)), y=origin.y+static_cast<int>(std::round(top));
    const bool resized=dstW!=surf_w_ || dstH!=surf_h_;
    if(!ensureSurface(dstW,dstH)) { hide(); return; }
    full |= resized;
    if(full) paint_surface();
    BLENDFUNCTION bf{}; bf.BlendOp=AC_SRC_OVER; bf.SourceConstantAlpha=255; bf.AlphaFormat=AC_SRC_ALPHA;
    // present using UpdateLayeredWindow (per-pixel alpha)
    HDC screen = GetDC(nullptr);
    SIZE siz{dstW, dstH}; POINT ptSrc{0,0}; POINT ptDst{x,y};
    BOOL updOk = UpdateLayeredWindow(hwnd_, screen, &ptDst, &siz, memdc_, &ptSrc, 0, &bf, ULW_ALPHA);
    if(!updOk){
        DWORD err = GetLastError();
        std::cerr << "OverlayWindow::redraw: UpdateLayeredWindow failed, err=" << err
                  << " w=" << dstW << " h=" << dstH << " x=" << x << " y=" << y << std::endl;
        ReleaseDC(nullptr,screen);
        hide();
        return;
    }
    ReleaseDC(nullptr, screen);
    ShowWindow(hwnd_, SW_SHOWNA);
}

void OverlayWindow::paint_surface(){
    const int dstW=surf_w_, dstH=surf_h_;
    // clear to fully transparent
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


        std::memset(bits_, 0, (size_t)(stride_ * dstH));

        // draw overlay marks: green boxes for Safe, red X for Mine, orange box for Guess
        const int w = geom_.board_w;
        const int h = geom_.board_h;
        if(w>0 && h>0 && (int)marks_.size()==w*h){
        std::vector<int> x0(w+1), y0(h+1);
        for(int xx=0; xx<=w; ++xx){ x0[xx] = (int)std::round(((double)dstW * (double)xx) / (double)w); }
        for(int yy=0; yy<=h; ++yy){ y0[yy] = (int)std::round(((double)dstH * (double)yy) / (double)h); }
        cached_xEdge_ = x0; cached_yEdge_ = y0;
        cached_dstW_ = dstW; cached_dstH_ = dstH;

        // fast check: if no marks to draw, skip per-cell drawing
        bool anyMarks=false; for(const auto m : marks_){ if(m!=solve::Mark::None){ anyMarks=true; break; } }
        if(anyMarks){
        
        // cell size in target surface pixels
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
            x0=std::max(0,x0); y0=std::max(0,y0);
            x1=std::min(dstW-1,x1); y1=std::min(dstH-1,y1);
            if(x1<x0 || y1<y0) return;
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
            if(n<=0) { putPixelV((x0+x1)/2,(y0+y1)/2,colv); return; }
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
                if(cx1<cx0 || cy1<cy0) continue;
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
                } else if(m==solve::Mark::ChordReady){
                    // yellow 'F' when all supporting flags are already placed
                    drawGlyphFScaled(cx0, cy0, cx1, cy1, 0xFFFFD700);
                } else if(m==solve::Mark::FlagForChord){
                    // blue filled square on mines that need flagging for the shown chord
                    fillRect(cx0, cy0, cx1, cy1, 0x661E90FF);
                    drawRectStroke(cx0, cy0, cx1, cy1, 0xFF1E90FF);
                } else if(m==solve::Mark::FlagForChordReady){
                    // yellow overlay for chord-supporting tiles that are already flagged
                    fillRect(cx0, cy0, cx1, cy1, 0x66FFD700);
                    drawRectStroke(cx0, cy0, cx1, cy1, 0xFFFFD700);
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
        } else if(ch=='?'){
            const unsigned char glyph[7]={0x0E,0x11,0x01,0x02,0x04,0x00,0x04};
            for(int r=0;r<7;++r) for(int c=0;c<5;++c) if(glyph[r]&(1<<(4-c))) plot(x0+c,y0+r);
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
            std::string label = "M=" + (mines_total_ < 0 ? std::string("?") : std::to_string(mines_total_));
            drawText(2, 2, label, 0xFFFFFFFF);
            if(!excluded_from_capture_){
                drawText(2, 12, std::string("VISIBLE"), 0xFFFFFFFF);
            }
            if(safety_mode_){
                drawText(2, 22, std::string("SAFE"), 0xFF00FF00);
            }
        }



}

LRESULT CALLBACK OverlayWindow::MouseHook(int code, WPARAM wp, LPARAM lp) {
    auto* self=hook_owner_;
    if(code>=0 && self && self->safety_mode_) {
        if(wp==WM_LBUTTONUP && self->blocked_left_) {
            self->blocked_left_=false; return 1;
        }
        if(wp==WM_LBUTTONDOWN) {
            self->blocked_left_=false;
            if(IsWindowVisible(self->hwnd_) && self->findRenderHost()) {
                const auto* event=reinterpret_cast<const MSLLHOOKSTRUCT*>(lp);
                POINT p=event->pt; ScreenToClient(self->hwnd_,&p);
                const auto& xs=self->cached_xEdge_;
                const auto& ys=self->cached_yEdge_;
                if(xs.size()>1 && ys.size()>1 && p.x>=0 && p.y>=0 && p.x<xs.back() && p.y<ys.back()) {
                    const int x=static_cast<int>(std::upper_bound(xs.begin(),xs.end(),p.x)-xs.begin())-1;
                    const int y=static_cast<int>(std::upper_bound(ys.begin(),ys.end(),p.y)-ys.begin())-1;
                    const size_t i=static_cast<size_t>(y)*self->geom_.board_w+x;
                    const auto mark=i<self->marks_.size() ? self->marks_[i] : solve::Mark::None;
                    if(mark!=solve::Mark::Safe && mark!=solve::Mark::Guess && mark!=solve::Mark::ChordReady) {
                        self->blocked_left_=true; return 1;
                    }
                }
            }
        }
    }
    return CallNextHookEx(nullptr,code,wp,lp);
}

LRESULT OverlayWindow::WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp){
    switch(msg) {
    case WM_NCHITTEST: return HTTRANSPARENT;
    case WM_MOUSEACTIVATE: return MA_NOACTIVATE;
    case WM_ERASEBKGND: return 1;
    default: return DefWindowProcW(hwnd,msg,wp,lp);
    }
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
        redraw(true);
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
    visible_ = visible;
    if(visible){ redraw(true); }
    else { hide(); }
#else
    (void)visible;
#endif
}

bool OverlayWindow::is_visible() const{
    return visible_;
}

void OverlayWindow::set_safety_mode(bool enabled){
#ifdef _WIN32
    if(enabled==safety_mode_) return;
    if(enabled) {
        if(hook_owner_ && hook_owner_!=this) return;
        mouse_hook_=SetWindowsHookExW(WH_MOUSE_LL,&OverlayWindow::MouseHook,GetModuleHandleW(nullptr),0);
        if(!mouse_hook_) { std::cerr<<"Failed to enable mouse safety: "<<GetLastError()<<'\n'; return; }
        hook_owner_=this;
    } else {
        if(mouse_hook_) UnhookWindowsHookEx(mouse_hook_);
        mouse_hook_=nullptr; hook_owner_=nullptr; blocked_left_=false;
    }
    // Keep WS_EX_TRANSPARENT so marked pixels never capture scroll/right-click
    // events. A hook blocks only unsafe left-button gestures across processes.
    safety_mode_=enabled;
    redraw(true);
#else
    safety_mode_=enabled;
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
