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
thread_local std::vector<OverlayWindow*> OverlayWindow::event_owners_;
namespace {
constexpr UINT RefreshWindow=WM_APP+1;
struct Region {
    HRGN handle;
    explicit Region(HRGN value=nullptr):handle(value) {}
    Region(const Region&)=delete;
    Region& operator=(const Region&)=delete;
    ~Region() { if(handle) DeleteObject(handle); }
    HRGN release() { HRGN value=handle; handle=nullptr; return value; }
};

bool client_rect(HWND window,RECT& rect) {
    if(!GetClientRect(window,&rect)) return false;
    SetLastError(0);
    return MapWindowPoints(window,nullptr,reinterpret_cast<POINT*>(&rect),2)!=0 || GetLastError()==0;
}

HRGN reflect_region(HRGN region,LONG width) {
    const DWORD bytes=GetRegionData(region,0,nullptr);
    if(bytes<sizeof(RGNDATAHEADER)) return nullptr;
    std::vector<unsigned char> buffer(bytes);
    auto* data=reinterpret_cast<RGNDATA*>(buffer.data());
    if(GetRegionData(region,bytes,data)!=bytes) return nullptr;
    auto* rectangles=reinterpret_cast<RECT*>(data->Buffer);
    for(DWORD i=0;i<data->rdh.nCount;++i) {
        const long long left=static_cast<long long>(width)-rectangles[i].right;
        const long long right=static_cast<long long>(width)-rectangles[i].left;
        if(left<LONG_MIN || right>LONG_MAX) return nullptr;
        rectangles[i].left=static_cast<LONG>(left); rectangles[i].right=static_cast<LONG>(right);
    }
    const LONG left=width-data->rdh.rcBound.right;
    data->rdh.rcBound.right=width-data->rdh.rcBound.left; data->rdh.rcBound.left=left;
    std::sort(rectangles,rectangles+data->rdh.nCount,[](const RECT& a,const RECT& b) {
        return a.top!=b.top ? a.top<b.top : a.left<b.left;
    });
    return ExtCreateRegion(nullptr,bytes,data);
}
}
#endif

OverlayWindow::OverlayWindow() {}
OverlayWindow::~OverlayWindow(){ destroy(); }

bool OverlayWindow::create(){
#ifdef _WIN32
    if(hwnd_) return true;
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
    const LONG_PTR style=GetWindowLongPtrW(hwnd_,GWL_EXSTYLE);
    if(style&WS_EX_LAYOUTRTL) {
        SetLastError(0);
        if(!SetWindowLongPtrW(hwnd_,GWL_EXSTYLE,style&~WS_EX_LAYOUTRTL) && GetLastError()!=0) {
            DestroyWindow(hwnd_); hwnd_=nullptr; return false;
        }
    }
    query_region_=CreateRectRgn(0,0,0,0);
    if(!query_region_) { DestroyWindow(hwnd_); hwnd_=nullptr; return false; }
    event_owners_.push_back(this);
    foreground_hook_=SetWinEventHook(EVENT_SYSTEM_FOREGROUND,EVENT_SYSTEM_FOREGROUND,nullptr,
        &OverlayWindow::WindowEvent,0,0,WINEVENT_OUTOFCONTEXT);
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
    event_owners_.erase(std::remove(event_owners_.begin(),event_owners_.end(),this),event_owners_.end());
    for(auto hook : {foreground_hook_,location_hook_,visibility_hook_}) if(hook) UnhookWinEvent(hook);
    foreground_hook_=location_hook_=visibility_hook_=nullptr;
    watched_pid_=0; refresh_pending_=false;
    if(mouse_hook_){ UnhookWindowsHookEx(mouse_hook_); mouse_hook_=nullptr; }
    if(hook_owner_==this) hook_owner_=nullptr;
    blocked_left_=false;
    if(memdc_){ DeleteDC(memdc_); memdc_=nullptr; }
    if(dib_){ DeleteObject(dib_); dib_=nullptr; }
    bits_ = nullptr; stride_ = 0; surf_w_ = surf_h_ = 0;
    surface_painted_=false; surface_dirty_=true; painted_marks_.clear();
    cached_xEdge_.clear(); cached_yEdge_.clear();
    cached_host_=cached_top_=nullptr;
    if(query_region_) DeleteObject(query_region_);
    if(last_host_region_) DeleteObject(last_host_region_);
    if(visible_region_) DeleteObject(visible_region_);
    query_region_=last_host_region_=visible_region_=nullptr;
    has_clip_=false;
    safety_mode_=false;
    if(hwnd_){ DestroyWindow(hwnd_); hwnd_=nullptr; }
#endif
}

void OverlayWindow::update(const std::vector<solve::Mark>& marks, const OverlayGeometry& geom, int minesTotal){
    bool marksChanged = marks != marks_;
    bool minesChanged = (minesTotal != mines_total_);
    bool sizeChanged = (geom.board_w != geom_.board_w) || (geom.board_h != geom_.board_h)
        || geom.rect_w!=geom_.rect_w || geom.rect_h!=geom_.rect_h
        || geom.dpr!=geom_.dpr || geom.vv_scale!=geom_.vv_scale;
    bool positionChanged = geom.rect_l!=geom_.rect_l || geom.rect_t!=geom_.rect_t
        || geom.vv_x!=geom_.vv_x || geom.vv_y!=geom_.vv_y;
    const bool clipChanged = geom.has_clip!=geom_.has_clip || (geom.has_clip &&
        (geom.clip_l!=geom_.clip_l || geom.clip_t!=geom_.clip_t ||
            geom.clip_w!=geom_.clip_w || geom.clip_h!=geom_.clip_h));

    if(!marksChanged && !minesChanged && !sizeChanged && !positionChanged && !clipChanged){
        return;
    }

    if(marksChanged) marks_ = marks;
    geom_ = geom;
    mines_total_ = minesTotal;

    if(!marksChanged && !minesChanged && !sizeChanged){
        redraw(false);
    } else {
        redraw(true);
    }
}

void OverlayWindow::tick(){
#ifdef _WIN32
    if(!hwnd_) return;
    HWND host=findRenderHost();
    if(host) watch_browser(cached_top_);
    if(!visible_ || !host || geom_.board_w<=0 || geom_.board_h<=0 || geom_.rect_w<=0 || geom_.rect_h<=0) { hide(); return; }
    POINT origin{0,0};
    RECT viewport{};
    Region region;
    if(!host_viewport(host,origin,viewport,region.handle)) { hide(); return; }
    if(!IsWindowVisible(hwnd_)) redraw(false);
    else if(origin.x!=last_host_origin_.x || origin.y!=last_host_origin_.y ||
        !EqualRect(&viewport,&last_host_rect_) || bool(region.handle)!=bool(last_host_region_) ||
        (region.handle && !EqualRgn(region.handle,last_host_region_))) redraw(false);
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
    HWND foreground=foreground_window_();
    if(!foreground || foreground==hwnd_ || IsIconic(foreground)) return nullptr;
    wchar_t cls[128]{};
    GetClassNameW(foreground,cls,128);
    if(lstrcmpW(cls,L"Chrome_WidgetWin_1")!=0) return nullptr;
    if(foreground!=cached_top_ || !cached_host_ || !IsWindow(cached_host_) ||
        !IsWindowVisible(cached_host_) || !IsChild(foreground,cached_host_)) {
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

bool OverlayWindow::host_viewport(HWND host, POINT& origin, RECT& viewport, HRGN& region) const {
    if(!host || !cached_top_ || !query_region_ || !client_rect(host,viewport)) return false;
    origin={viewport.left,viewport.top};
    HWND current=host;
    for(int depth=0;depth<128;++depth) {
        RECT bounds=current==host ? viewport : RECT{},intersection{};
        if(!current || !IsWindowVisible(current) || (current!=host && !client_rect(current,bounds)) ||
            !IntersectRect(&intersection,&viewport,&bounds)) return false;
        viewport=intersection;
        const int kind=GetWindowRgn(current,query_region_);
        if(kind==NULLREGION) return false;
        // ERROR also means that an ordinary window has no custom region.
        if(kind!=ERROR) {
            if(!GetWindowRect(current,&bounds)) return false;
            Region mirrored;
            HRGN shape=query_region_;
            if(GetWindowLongPtrW(current,GWL_EXSTYLE)&WS_EX_LAYOUTRTL) {
                const long long width=static_cast<long long>(bounds.right)-bounds.left;
                if(width<=0 || width>LONG_MAX) return false;
                mirrored.handle=reflect_region(shape,static_cast<LONG>(width));
                if(!mirrored.handle) return false;
                shape=mirrored.handle;
            }
            if(!region) region=CreateRectRgn(viewport.left,viewport.top,viewport.right,viewport.bottom);
            if(!region || OffsetRgn(shape,bounds.left,bounds.top)<=NULLREGION ||
                CombineRgn(region,region,shape,RGN_AND)<=NULLREGION) return false;
        }
        if(current==cached_top_) {
            if(!region) { OffsetRect(&viewport,-origin.x,-origin.y); return true; }
            return SetRectRgn(query_region_,viewport.left,viewport.top,viewport.right,viewport.bottom) &&
                CombineRgn(region,region,query_region_,RGN_AND)>NULLREGION &&
                OffsetRgn(region,-origin.x,-origin.y)>NULLREGION && GetRgnBox(region,&viewport)>NULLREGION;
        }
        current=GetAncestor(current,GA_PARENT);
    }
    return false;
}

void OverlayWindow::watch_browser(HWND browser) {
    DWORD pid=0;
    GetWindowThreadProcessId(browser,&pid);
    if(pid==watched_pid_) return;
    if(location_hook_) UnhookWinEvent(location_hook_);
    if(visibility_hook_) UnhookWinEvent(visibility_hook_);
    watched_pid_=pid;
    location_hook_=visibility_hook_=nullptr;
    if(!pid) return;
    location_hook_=SetWinEventHook(EVENT_OBJECT_LOCATIONCHANGE,EVENT_OBJECT_LOCATIONCHANGE,nullptr,
        &OverlayWindow::WindowEvent,pid,0,WINEVENT_OUTOFCONTEXT);
    visibility_hook_=SetWinEventHook(EVENT_OBJECT_DESTROY,EVENT_OBJECT_HIDE,nullptr,
        &OverlayWindow::WindowEvent,pid,0,WINEVENT_OUTOFCONTEXT);
}

void CALLBACK OverlayWindow::WindowEvent(HWINEVENTHOOK hook,DWORD event,HWND window,LONG object,LONG child,DWORD,DWORD) {
    // Out-of-context callbacks run on the registering thread. Defer and coalesce
    // presentation so reentrant notifications never paint inside this callback.
    for(auto* self : event_owners_) {
        if(hook!=self->foreground_hook_ && hook!=self->location_hook_ && hook!=self->visibility_hook_) continue;
        if(window==self->hwnd_) return;
        if(event!=EVENT_SYSTEM_FOREGROUND) {
            if(object!=OBJID_WINDOW || child!=CHILDID_SELF) return;
            if(window!=self->cached_top_ && window!=self->cached_host_ && !IsChild(window,self->cached_host_) &&
                !(event==EVENT_OBJECT_SHOW && IsChild(self->cached_top_,window))) return;
        }
        if(!self->refresh_pending_) self->refresh_pending_=PostMessageW(self->hwnd_,RefreshWindow,0,0)!=FALSE;
        return;
    }
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
    surface_painted_=false;
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

RECT OverlayWindow::visible_board_rect(int x, int y, int w, int h, const RECT& viewport){
    return {std::clamp(static_cast<int>(viewport.left)-x,0,w),std::clamp(static_cast<int>(viewport.top)-y,0,h),
        std::clamp(static_cast<int>(viewport.right)-x,0,w),std::clamp(static_cast<int>(viewport.bottom)-y,0,h)};
}

bool OverlayWindow::board_bounds(RECT& bounds) const {
    if(!game::valid_dimensions(geom_.board_w,geom_.board_h) || geom_.board_w==0 || geom_.board_h==0) return false;
    const double scalePx=geom_.dpr*geom_.vv_scale;
    const double width=geom_.rect_w*scalePx, height=geom_.rect_h*scalePx;
    const double left=(geom_.rect_l-geom_.vv_x)*scalePx, top=(geom_.rect_t-geom_.vv_y)*scalePx;
    if(geom_.dpr<=0 || geom_.vv_scale<=0 || !std::isfinite(width) || !std::isfinite(height) || !std::isfinite(left) || !std::isfinite(top) ||
        width<0.5 || height<0.5 || width>8192 || height>8192 || std::abs(left)>1000000 || std::abs(top)>1000000) return false;
    const int dstW=static_cast<int>(std::round(width)), dstH=static_cast<int>(std::round(height));
    if(static_cast<uint64_t>(dstW)*dstH>16*1024*1024) return false;
    const int x=static_cast<int>(std::round(left)), y=static_cast<int>(std::round(top));
    bounds={x,y,x+dstW,y+dstH};
    return true;
}

RECT OverlayWindow::board_clip(const RECT& bounds, const RECT& viewport) const {
    const int width=bounds.right-bounds.left, height=bounds.bottom-bounds.top;
    const RECT native=visible_board_rect(bounds.left,bounds.top,width,height,viewport);
    if(!geom_.has_clip) return native;
    const double scale=geom_.dpr*geom_.vv_scale;
    const double left=(geom_.clip_l-geom_.vv_x)*scale-bounds.left;
    const double top=(geom_.clip_t-geom_.vv_y)*scale-bounds.top;
    const double right=(geom_.clip_l+geom_.clip_w-geom_.vv_x)*scale-bounds.left;
    const double bottom=(geom_.clip_t+geom_.clip_h-geom_.vv_y)*scale-bounds.top;
    if(geom_.clip_w<0 || geom_.clip_h<0 || scale<=0 || !std::isfinite(left) ||
        !std::isfinite(top) || !std::isfinite(right) || !std::isfinite(bottom)) return {};
    // Keep only complete device pixels inside the browser's CSS intersection.
    const RECT dom{static_cast<LONG>(std::ceil(std::clamp(left,0.0,double(width)))),
        static_cast<LONG>(std::ceil(std::clamp(top,0.0,double(height)))),
        static_cast<LONG>(std::floor(std::clamp(right,0.0,double(width)))),
        static_cast<LONG>(std::floor(std::clamp(bottom,0.0,double(height))))};
    RECT clip{};
    IntersectRect(&clip,&native,&dom);
    return clip;
}

bool OverlayWindow::board_region(const RECT& bounds, const RECT& viewport, HRGN host, HRGN region) const {
    const RECT clip=board_clip(bounds,viewport);
    return region && !IsRectEmpty(&clip) &&
        SetRectRgn(region,clip.left+bounds.left,clip.top+bounds.top,clip.right+bounds.left,clip.bottom+bounds.top) &&
        (!host || CombineRgn(region,region,host,RGN_AND)>NULLREGION) &&
        OffsetRgn(region,-bounds.left,-bounds.top)>NULLREGION;
}

void OverlayWindow::redraw(bool full){
    surface_dirty_|=full;
    if(!hwnd_) return;
    HWND host=findRenderHost();
    if(host) watch_browser(cached_top_);
    RECT bounds{};
    if(!visible_ || !host || !board_bounds(bounds)) { hide(); return; }
    const int dstW=bounds.right-bounds.left, dstH=bounds.bottom-bounds.top;
    POINT origin{0,0};
    RECT viewport{};
    Region hostRegion,clip(CreateRectRgn(0,0,0,0));
    if(!host_viewport(host,origin,viewport,hostRegion.handle)) { hide(); return; }
    last_host_origin_=origin;
    last_host_rect_=viewport;
    if(last_host_region_) DeleteObject(last_host_region_);
    last_host_region_=hostRegion.release();
    const int boardX=bounds.left, boardY=bounds.top;
    if(!board_region(bounds,viewport,last_host_region_,clip.handle)) { hide(); return; }
    if(!visible_region_ || !EqualRgn(clip.handle,visible_region_)) {
        Region applied(CreateRectRgn(0,0,0,0));
        if(!applied.handle || CombineRgn(applied.handle,clip.handle,nullptr,RGN_COPY)==ERROR ||
            !SetWindowRgn(hwnd_,applied.handle,TRUE)) { hide(); return; }
        // Windows owns the region after a successful SetWindowRgn call.
        applied.release();
        if(visible_region_) DeleteObject(visible_region_);
        visible_region_=clip.release();
        GetRgnBox(visible_region_,&visible_rect_); has_clip_=true;
    }
    const int x=origin.x+boardX, y=origin.y+boardY;
    const bool resized=dstW!=surf_w_ || dstH!=surf_h_;
    if(!ensureSurface(dstW,dstH)) { hide(); return; }
    full = surface_dirty_ || resized || !surface_painted_;
    if(!full && SetWindowPos(hwnd_,nullptr,x,y,0,0,SWP_NOSIZE|SWP_NOZORDER|SWP_NOACTIVATE)) {
        ShowWindow(hwnd_,SW_SHOWNA);
        return;
    }
    if(full) { paint_surface(); surface_dirty_=false; }
    BLENDFUNCTION bf{}; bf.BlendOp=AC_SRC_OVER; bf.SourceConstantAlpha=255; bf.AlphaFormat=AC_SRC_ALPHA;
    // present using UpdateLayeredWindow (per-pixel alpha)
    HDC screen = GetDC(nullptr);
    SIZE siz{dstW, dstH}; POINT ptSrc{0,0}; POINT ptDst{x,y};
    BOOL updOk = UpdateLayeredWindow(hwnd_, screen, &ptDst, &siz, memdc_, &ptSrc, 0, &bf, ULW_ALPHA);
    if(!updOk){
        surface_dirty_=true;
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
    const int w=geom_.board_w, h=geom_.board_h;
    if(!bits_ || dstW<=0 || dstH<=0) return;
    const bool validBoard=game::valid_dimensions(w,h) && w>0 && h>0;
    const bool geometryChanged=cached_dstW_!=dstW || cached_dstH_!=dstH ||
        (validBoard && (cached_xEdge_.size()!=static_cast<size_t>(w+1) || cached_yEdge_.size()!=static_cast<size_t>(h+1)));
    if(validBoard && geometryChanged) {
        cached_xEdge_.resize(w+1); cached_yEdge_.resize(h+1);
        for(int x=0;x<=w;++x) cached_xEdge_[x]=static_cast<int>(std::round(double(dstW)*x/w));
        for(int y=0;y<=h;++y) cached_yEdge_[y]=static_cast<int>(std::round(double(dstH)*y/h));
    } else if(!validBoard) {
        cached_xEdge_.clear(); cached_yEdge_.clear();
    }
    cached_dstW_=dstW; cached_dstH_=dstH;
    const bool validMarks=validBoard && marks_.size()==static_cast<size_t>(w)*h;
    bool full=!surface_painted_ || geometryChanged || !validMarks || painted_marks_.size()!=marks_.size() ||
        painted_mines_total_!=mines_total_ || painted_exclusion_!=excluded_from_capture_ || painted_safety_!=safety_mode_;
    if(!full) {
        size_t changed=0;
        for(size_t i=0;i<marks_.size();++i) changed+=marks_[i]!=painted_marks_[i];
        if(changed==0) return;
        full=changed>marks_.size()/3;
    }
    if(full) std::memset(bits_,0,static_cast<size_t>(stride_)*dstH);
    int clipLeft=0,clipTop=0,clipRight=dstW,clipBottom=dstH;
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
        if(px<clipLeft || py<clipTop || px>=clipRight || py>=clipBottom) return;
        uint8_t* p = bits_ + py*stride_ + px*4;
        *reinterpret_cast<uint32_t*>(p) = bgra;
    };
        // draw overlay marks: green boxes for Safe, red X for Mine, orange box for Guess
        if(validMarks){
        const auto& x0=cached_xEdge_;
        const auto& y0=cached_yEdge_;
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
                const size_t index=static_cast<size_t>(yy)*w+xx;
                solve::Mark m = marks_[index];
                if(!full && m==painted_marks_[index]) continue;
                int cx0 = x0[xx];
                int cx1 = x0[xx+1] - 1;
                if(cx1<cx0 || cy1<cy0) continue;
                if(!full) for(int y=cy0;y<=cy1;++y)
                    std::memset(bits_+y*stride_+cx0*4,0,static_cast<size_t>(cx1-cx0+1)*4);
                if(m==solve::Mark::None) continue;
                // Marks must stay inside their cell, even below the stroke/glyph size.
                clipLeft=cx0; clipRight=cx1+1; clipTop=cy0; clipBottom=cy1+1;
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
    clipLeft=clipTop=0; clipRight=dstW; clipBottom=dstH;
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
    painted_marks_=marks_;
    painted_mines_total_=mines_total_;
    painted_exclusion_=excluded_from_capture_;
    painted_safety_=safety_mode_;
    surface_painted_=true;
}

bool OverlayWindow::unsafe_at(POINT point, const RECT& clip) const {
    if(!PtInRect(&clip,point)) return false;
    const auto& xs=cached_xEdge_;
    const auto& ys=cached_yEdge_;
    if(xs.size()<2 || ys.size()<2 || point.x<0 || point.y<0 || point.x>=xs.back() || point.y>=ys.back()) return false;
    const int x=static_cast<int>(std::upper_bound(xs.begin(),xs.end(),point.x)-xs.begin())-1;
    const int y=static_cast<int>(std::upper_bound(ys.begin(),ys.end(),point.y)-ys.begin())-1;
    const size_t i=static_cast<size_t>(y)*geom_.board_w+x;
    const auto mark=i<marks_.size() ? marks_[i] : solve::Mark::None;
    return mark!=solve::Mark::Safe && mark!=solve::Mark::Guess && mark!=solve::Mark::ChordReady;
}

LRESULT CALLBACK OverlayWindow::MouseHook(int code, WPARAM wp, LPARAM lp) {
    auto* self=hook_owner_;
    if(code>=0 && self) {
        const auto* event=reinterpret_cast<const MSLLHOOKSTRUCT*>(lp);
        if(self->block_mouse(wp,event->pt)) return 1;
    }
    return CallNextHookEx(nullptr,code,wp,lp);
}

bool OverlayWindow::block_mouse(WPARAM message, POINT point) {
    if(!safety_mode_) return false;
    if(message==WM_LBUTTONUP && blocked_left_) {
        blocked_left_=false; return true;
    }
    if(message==WM_LBUTTONDOWN) {
        blocked_left_=false;
        HWND host=IsWindowVisible(hwnd_) ? findRenderHost() : nullptr;
        const HWND target=host ? WindowFromPoint(point) : nullptr;
        if(!target || (target!=host && !IsChild(host,target))) return false;
        POINT origin{};
        RECT viewport{},bounds{};
        Region hostRegion;
        if(host_viewport(host,origin,viewport,hostRegion.handle) && board_bounds(bounds) &&
            (!hostRegion.handle || PtInRegion(hostRegion.handle,point.x-origin.x,point.y-origin.y))) {
            // The browser can move or resize between presentation ticks.
            // Use its current coordinates and clip, not the overlay's last frame.
            point.x-=origin.x+bounds.left; point.y-=origin.y+bounds.top;
            const RECT clip=board_clip(bounds,viewport);
            if(unsafe_at(point,clip)) {
                blocked_left_=true; return true;
            }
        }
    }
    return false;
}

LRESULT OverlayWindow::WndProc(HWND hwnd, UINT msg, WPARAM wp, LPARAM lp){
    switch(msg) {
    case RefreshWindow: refresh_pending_=false; tick(); return 0;
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
