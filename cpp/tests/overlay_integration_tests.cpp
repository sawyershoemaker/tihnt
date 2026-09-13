#include "overlay.hpp"
#include <imm.h>

#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <stdexcept>
#include <string>
#include <thread>
#include <vector>

namespace {
_Post_satisfies_(value)
void require(bool value,const char* message) {
    if(!value) throw std::runtime_error(std::string(message)+" (Win32 error "+std::to_string(GetLastError())+")");
}

void check_isolated(HDESK desktop) {
    BOOL receivesInput=TRUE;
    require(GetUserObjectInformationW(desktop,UOI_IO,&receivesInput,sizeof(receivesInput),nullptr)!=FALSE,
        "Cannot verify desktop isolation");
    require(!receivesInput,"Test desktop unexpectedly receives user input");
}

struct Windows {
    HINSTANCE instance=GetModuleHandleW(nullptr);
    std::vector<HWND> handles;
    std::vector<const wchar_t*> classes;
    ~Windows() {
        for(auto it=handles.rbegin();it!=handles.rend();++it) if(IsWindow(*it)) DestroyWindow(*it);
        for(const wchar_t* name : classes) UnregisterClassW(name,instance);
    }
    void register_class(const wchar_t* name) {
        WNDCLASSW cls{}; cls.lpfnWndProc=DefWindowProcW; cls.hInstance=instance; cls.lpszClassName=name;
        require(RegisterClassW(&cls)!=0,"Synthetic window class registration failed");
        classes.push_back(name);
    }
    HWND create(const wchar_t* name,HWND parent,int x,int y,int width,int height,bool visible=true) {
        const DWORD style=(parent ? WS_CHILD : WS_POPUP)|(visible ? WS_VISIBLE : 0);
        HWND hwnd=CreateWindowExW(parent ? 0 : WS_EX_TOOLWINDOW,name,L"TIHNT isolated integration test",style,
            x,y,width,height,parent,nullptr,instance,nullptr);
        require(hwnd!=nullptr,"Synthetic browser window creation failed");
        handles.push_back(hwnd);
        return hwnd;
    }
};

struct TestRegion {
    HRGN handle;
    explicit TestRegion(RECT rect={}):handle(CreateRectRgn(rect.left,rect.top,rect.right,rect.bottom)) {
        require(handle!=nullptr,"Cannot allocate a test region");
    }
    TestRegion(const TestRegion&)=delete;
    TestRegion& operator=(const TestRegion&)=delete;
    ~TestRegion() { if(handle) DeleteObject(handle); }
    void apply(HWND window) {
        require(SetWindowRgn(window,handle,TRUE)!=0,"Cannot apply a native test region");
        handle=nullptr;
    }
};

void activate(HWND window,HDESK desktop) {
    check_isolated(desktop);
    // An inactive desktop has no foreground window. Its thread-local active
    // window supplies the test seam without touching the user's foreground.
    SetActiveWindow(window);
    require(GetActiveWindow()==window,"Cannot activate the isolated test window");
}

void move(HWND window,int x,int y,int width,int height) {
    require(SetWindowPos(window,nullptr,x,y,width,height,SWP_NOZORDER|SWP_NOACTIVATE)!=FALSE,
        "Synthetic host move or resize failed");
}

void drain_messages() {
    MSG message{};
    while(PeekMessageW(&message,nullptr,0,0,PM_REMOVE)) {
        TranslateMessage(&message);
        DispatchMessageW(&message);
    }
}

template<class Predicate>
void wait_for_window_event(Predicate ready,const char* message) {
    const ULONGLONG deadline=GetTickCount64()+2000;
    while(!ready() && GetTickCount64()<deadline) {
        drain_messages();
        if(!ready()) MsgWaitForMultipleObjectsEx(0,nullptr,10,QS_ALLINPUT,MWMO_INPUTAVAILABLE);
    }
    require(ready(),message);
}
}

struct OverlayTestAccess {
    static std::vector<unsigned char> pixels(const OverlayWindow& overlay) {
        require(overlay.bits_!=nullptr,"Presented overlay has no backing pixels");
        return {overlay.bits_,overlay.bits_+static_cast<size_t>(overlay.stride_)*overlay.surf_h_};
    }

    static void native_regions(HDESK desktop) {
        Windows windows;
        windows.register_class(L"Chrome_WidgetWin_1");
        windows.register_class(L"Chrome_RenderWidgetHostHWND");
        windows.register_class(L"TIHNT_ViewportContainer");
        const HWND browser=windows.create(L"Chrome_WidgetWin_1",nullptr,100,100,400,300);
        const HWND container=windows.create(L"TIHNT_ViewportContainer",browser,20,30,300,220);
        const HWND host=windows.create(L"Chrome_RenderWidgetHostHWND",container,10,20,240,180);
        activate(browser,desktop);
        OverlayWindow overlay; overlay.foreground_window_=&GetActiveWindow;
        DWORD layout=0;
        require(GetProcessDefaultLayout(&layout)!=FALSE && SetProcessDefaultLayout(LAYOUT_RTL)!=FALSE,
            "Cannot configure the isolated process layout");
        const bool created=overlay.create();
        require(SetProcessDefaultLayout(layout)!=FALSE,"Cannot restore the isolated process layout");
        require(created && !(GetWindowLongPtrW(overlay.hwnd_,GWL_EXSTYLE)&WS_EX_LAYOUTRTL),
            "Overlay inherited a mirrored region coordinate system");
        OverlayGeometry geometry{}; geometry.board_w=4; geometry.board_h=3;
        geometry.rect_w=240; geometry.rect_h=180;
        std::vector<solve::Mark> marks(12,solve::Mark::Mine);
        overlay.update(marks,geometry,12);
        const HBITMAP bitmap=overlay.dib_; const HDC dc=overlay.memdc_;
        auto originalPixels=pixels(overlay);
        auto same_surface=[&] {
            require(overlay.dib_==bitmap && overlay.memdc_==dc && pixels(overlay)==originalPixels && !overlay.surface_dirty_,
                "Native region change replaced or repainted the surface");
        };
        auto host_client=[&] {
            RECT client{};
            require(GetClientRect(host,&client)!=FALSE,"Cannot inspect the shaped host client");
            SetLastError(0);
            require(MapWindowPoints(host,nullptr,reinterpret_cast<POINT*>(&client),2)!=0 || GetLastError()==0,
                "Cannot map the shaped host client");
            return client;
        };
        auto cutout=[&](HWND target,RECT hole,bool disjoint=false) {
            const RECT client=host_client(); RECT window{};
            require(GetWindowRect(target,&window)!=FALSE,"Cannot inspect a shaped ancestor window");
            const LONG width=window.right-window.left,height=window.bottom-window.top;
            OffsetRect(&hole,client.left-window.left,client.top-window.top);
            if(disjoint) { hole.top=0; hole.bottom=height; }
            if(GetWindowLongPtrW(target,GWL_EXSTYLE)&WS_EX_LAYOUTRTL) {
                const LONG left=width-hole.right; hole.right=width-hole.left; hole.left=left;
            }
            TestRegion region({0,0,width,height}),removed(hole);
            require(CombineRgn(region.handle,region.handle,removed.handle,RGN_DIFF)==COMPLEXREGION,
                "Synthetic region did not retain a hole or disjoint spans");
            region.apply(target);
        };
        auto check=[&](HRGN expected) {
            TestRegion actual;
            require(IsWindowVisible(overlay.hwnd_) && GetWindowRgn(overlay.hwnd_,actual.handle)!=ERROR &&
                EqualRgn(actual.handle,expected),"Overlay lost the exact native region intersection");
            RECT bounds{}; const RECT client=host_client();
            require(GetWindowRect(overlay.hwnd_,&bounds)!=FALSE && bounds.left==client.left && bounds.top==client.top &&
                bounds.right==client.left+240 && bounds.bottom==client.top+180,
                "Shaped overlay did not use the physical host client origin");
            require(overlay.visible_region_ && EqualRgn(overlay.visible_region_,actual.handle),
                "Owned presentation region differs from the window-owned region");
            same_surface();
        };
        auto mouse=[&](POINT point,bool expected) {
            const RECT client=host_client(); point.x+=client.left; point.y+=client.top;
            overlay.safety_mode_=true;
            const bool down=overlay.block_mouse(WM_LBUTTONDOWN,point),up=overlay.block_mouse(WM_LBUTTONUP,point);
            overlay.safety_mode_=false;
            require(down==expected && up==expected,"Mouse safety disagrees with the fresh native region");
        };
        drain_messages();
        require(UnhookWinEvent(overlay.location_hook_)!=FALSE,"Cannot isolate region polling from event delivery");
        overlay.location_hook_=nullptr;
        for(HWND target : {host,container,browser}) for(bool border : {false,true}) for(bool rtl : {false,true}) {
            const LONG_PTR style=GetWindowLongPtrW(target,GWL_STYLE),extended=GetWindowLongPtrW(target,GWL_EXSTYLE);
            RECT window{},client{};
            require(GetWindowRect(target,&window)!=FALSE && GetClientRect(target,&client)!=FALSE,"Cannot save native styles");
            SetWindowLongPtrW(target,GWL_STYLE,border ? style|WS_BORDER : style);
            SetWindowLongPtrW(target,GWL_EXSTYLE,rtl ? extended|WS_EX_LAYOUTRTL|WS_EX_NOINHERITLAYOUT : extended);
            require(SetWindowPos(target,nullptr,0,0,0,0,SWP_NOMOVE|SWP_NOSIZE|SWP_NOZORDER|SWP_NOACTIVATE|SWP_FRAMECHANGED)!=FALSE,
                "Cannot update native region test styles");
            RECT resized{};
            require(GetClientRect(target,&resized)!=FALSE,"Cannot inspect bordered client size");
            require(SetWindowPos(target,nullptr,0,0,window.right-window.left+client.right-resized.right,
                window.bottom-window.top+client.bottom-resized.bottom,SWP_NOMOVE|SWP_NOZORDER|SWP_NOACTIVATE)!=FALSE,
                "Cannot preserve the client size while adding a border");
            for(bool disjoint : {false,true}) {
                const RECT hole{60,disjoint ? 0 : 60,120,disjoint ? 180 : 120};
                TestRegion expected({0,0,240,180}),removed(hole);
                require(CombineRgn(expected.handle,expected.handle,removed.handle,RGN_DIFF)==COMPLEXREGION,
                    "Cannot construct the expected native intersection");
                cutout(target,hole,disjoint); overlay.tick(); check(expected.handle);
                mouse({30,90},true); mouse({90,90},false); mouse({150,90},true);
                const RECT movedHole{120,hole.top,180,hole.bottom};
                cutout(target,movedHole,disjoint);
                mouse({90,90},true); mouse({150,90},false);
                drain_messages();
                require(overlay.visible_region_ && EqualRgn(overlay.visible_region_,expected.handle),
                    "Region polling control received an unexpected refresh");
                TestRegion shifted({0,0,240,180}),shiftedRemoved(movedHole);
                require(CombineRgn(shifted.handle,shifted.handle,shiftedRemoved.handle,RGN_DIFF)==COMPLEXREGION,
                    "Cannot construct the shifted native intersection");
                overlay.tick(); check(shifted.handle);
                require(SetWindowRgn(target,nullptr,TRUE)!=0,"Cannot clear native test region");
                overlay.tick(); TestRegion whole({0,0,240,180}); check(whole.handle);
            }
            SetWindowLongPtrW(target,GWL_STYLE,style); SetWindowLongPtrW(target,GWL_EXSTYLE,extended);
            require(SetWindowPos(target,nullptr,0,0,window.right-window.left,window.bottom-window.top,
                SWP_NOMOVE|SWP_NOZORDER|SWP_NOACTIVATE|SWP_FRAMECHANGED)!=FALSE,"Cannot restore native window styles");
            overlay.tick();
        }
        const RECT holes[]={{60,60,120,120},{140,20,180,160},{80,0,100,180}};
        cutout(host,holes[0]); cutout(container,holes[1]); cutout(browser,holes[2],true);
        geometry.has_clip=true; geometry.clip_l=40.25; geometry.clip_t=30.5;
        geometry.clip_w=120.5; geometry.clip_h=120.25;
        overlay.update(marks,geometry,12);
        TestRegion combined({41,31,160,150});
        for(const RECT& hole : holes) {
            TestRegion removed(hole);
            require(CombineRgn(combined.handle,combined.handle,removed.handle,RGN_DIFF)>NULLREGION,
                "Cannot combine native and DOM test regions");
        }
        check(combined.handle); mouse({50,40},true); mouse({100,90},false); mouse({150,100},false); mouse({200,50},false);
        for(HWND target : {host,container,browser}) require(SetWindowRgn(target,nullptr,TRUE)!=0,"Cannot clear combined native regions");
        geometry.has_clip=false; overlay.update(marks,geometry,12);
        {
            TestRegion left({0,0,60,180}); left.apply(host);
            RECT parent{}; const RECT client=host_client();
            require(GetWindowRect(container,&parent)!=FALSE,"Cannot map disjoint ancestor regions");
            const LONG x=client.left-parent.left,y=client.top-parent.top;
            TestRegion right({x+180,y,x+240,y+180}); right.apply(container);
            overlay.tick();
            require(!IsWindowVisible(overlay.hwnd_),"Nonempty disjoint native regions retained visible hints");
            same_surface(); mouse({30,90},false);
            require(SetWindowRgn(host,nullptr,TRUE)!=0 && SetWindowRgn(container,nullptr,TRUE)!=0,
                "Cannot restore disjoint native regions");
            overlay.tick(); TestRegion whole({0,0,240,180}); check(whole.handle);
        }
        {
            TestRegion empty; empty.apply(host); overlay.tick();
            require(!IsWindowVisible(overlay.hwnd_) && IsWindowVisible(host),"Explicit empty native region retained hints");
            mouse({30,90},false);
            marks[0]=solve::Mark::Safe; overlay.update(marks,geometry,11);
            require(overlay.surface_dirty_ && pixels(overlay)==originalPixels,"Hidden native-region update lost pending marks");
            require(SetWindowRgn(host,nullptr,TRUE)!=0,"Cannot restore an empty native region");
            overlay.tick();
            require(overlay.painted_marks_==marks && pixels(overlay)!=originalPixels && !overlay.surface_dirty_,
                "Native-region restoration reused stale marks");
            originalPixels=pixels(overlay);
            TestRegion whole({0,0,240,180}); check(whole.handle);
        }
        cutout(host,{60,60,120,120}); overlay.tick();
        const DWORD gdiBefore=GetGuiResources(GetCurrentProcess(),GR_GDIOBJECTS);
        for(int i=0;i<160;++i) {
            cutout(host,{60+i%30,60,120+i%30,120}); overlay.tick(); same_surface();
        }
        require(GetGuiResources(GetCurrentProcess(),GR_GDIOBJECTS)<=gdiBefore+1,"Region replacement leaked GDI resources");
        const HRGN query=overlay.query_region_,hostRegion=overlay.last_host_region_,visible=overlay.visible_region_;
        require(query && hostRegion && visible,"Native region ownership was not retained");
        overlay.destroy();
        require(!overlay.query_region_ && !overlay.last_host_region_ && !overlay.visible_region_ &&
            !GetObjectType(query) && !GetObjectType(hostRegion) && !GetObjectType(visible),"Overlay destruction retained owned regions");
        check_isolated(desktop);
    }

    static void native_region_coordinates(HDESK desktop) {
        Windows windows;
        windows.register_class(L"Chrome_WidgetWin_1");
        windows.register_class(L"Chrome_RenderWidgetHostHWND");
        windows.register_class(L"TIHNT_ViewportContainer");
        const HWND browser=windows.create(L"Chrome_WidgetWin_1",nullptr,100,100,400,180);
        const HWND container=windows.create(L"TIHNT_ViewportContainer",browser,0,0,400,180);
        const HWND host=windows.create(L"Chrome_RenderWidgetHostHWND",container,0,0,400,180);
        activate(browser,desktop);
        for(HWND target : {host,container,browser}) for(bool rtl : {false,true}) {
            const LONG_PTR style=GetWindowLongPtrW(target,GWL_EXSTYLE);
            SetWindowLongPtrW(target,GWL_EXSTYLE,rtl ? style|WS_EX_LAYOUTRTL|WS_EX_NOINHERITLAYOUT : style&~WS_EX_LAYOUTRTL);
            RECT bounds{};
            require(GetWindowRect(target,&bounds)!=FALSE,"Cannot inspect region coordinate target");
            TestRegion region({0,0,120,bounds.bottom-bounds.top}); region.apply(target);
            TestRegion query; RECT raw{};
            require(GetWindowRgn(target,query.handle)!=ERROR && GetRgnBox(query.handle,&raw)!=ERROR,"Cannot inspect region coordinates");
            const HWND left=WindowFromPoint({bounds.left+60,bounds.top+90});
            const HWND right=WindowFromPoint({bounds.right-60,bounds.top+90});
            const RECT expected{0,0,120,bounds.bottom-bounds.top};
            require(EqualRect(&raw,&expected) && bool(left==target || IsChild(target,left))!=rtl &&
                bool(right==target || IsChild(target,right))==rtl,"Native RTL region coordinates do not match physical hit testing");
            require(SetWindowRgn(target,nullptr,TRUE)!=0,"Cannot clear asymmetric region probe");
            SetWindowLongPtrW(target,GWL_EXSTYLE,style);
        }
        check_isolated(desktop);
    }

    static void check_presentation(const OverlayWindow& overlay,HWND host,const OverlayGeometry& geometry) {
        require(IsWindowVisible(overlay.hwnd_)!=FALSE,"Valid browser board did not show its overlay");
        POINT origin{}; RECT viewport{},actual{};
        require(ClientToScreen(host,&origin)!=FALSE && GetClientRect(host,&viewport)!=FALSE &&
            GetWindowRect(overlay.hwnd_,&actual)!=FALSE,"Cannot read presented window geometry");
        const double scale=geometry.dpr*geometry.vv_scale;
        const int x=static_cast<int>(std::round((geometry.rect_l-geometry.vv_x)*scale));
        const int y=static_cast<int>(std::round((geometry.rect_t-geometry.vv_y)*scale));
        const int width=static_cast<int>(std::round(geometry.rect_w*scale));
        const int height=static_cast<int>(std::round(geometry.rect_h*scale));
        const RECT expected{origin.x+x,origin.y+y,origin.x+x+width,origin.y+y+height};
        require(EqualRect(&actual,&expected)!=FALSE,"Overlay HWND bounds do not track the render host");
        require(overlay.surf_w_==width && overlay.surf_h_==height,"Presented surface size disagrees with window bounds");

        HRGN region=CreateRectRgn(0,0,0,0);
        require(region!=nullptr,"Cannot allocate region query");
        const int regionType=GetWindowRgn(overlay.hwnd_,region);
        RECT actualClip{};
        const int boundsType=GetRgnBox(region,&actualClip);
        RECT expectedClip=expected;
        if(geometry.has_clip) {
            const RECT dom{origin.x+static_cast<LONG>(std::ceil((geometry.clip_l-geometry.vv_x)*scale)),
                origin.y+static_cast<LONG>(std::ceil((geometry.clip_t-geometry.vv_y)*scale)),
                origin.x+static_cast<LONG>(std::floor((geometry.clip_l+geometry.clip_w-geometry.vv_x)*scale)),
                origin.y+static_cast<LONG>(std::floor((geometry.clip_t+geometry.clip_h-geometry.vv_y)*scale))};
            require(IntersectRect(&expectedClip,&expected,&dom)!=FALSE,"Expected board is clipped by DOM geometry");
        }
        const HWND root=GetAncestor(host,GA_ROOT);
        for(HWND ancestor=host;ancestor;ancestor=GetAncestor(ancestor,GA_PARENT)) {
            RECT client{}; POINT offset{};
            require(GetClientRect(ancestor,&client)!=FALSE && ClientToScreen(ancestor,&offset)!=FALSE,
                "Cannot inspect ancestor client geometry");
            OffsetRect(&client,offset.x,offset.y);
            RECT intersection{};
            require(IntersectRect(&intersection,&expectedClip,&client)!=FALSE,"Expected board is clipped by an ancestor");
            expectedClip=intersection;
            if(ancestor==root) break;
        }
        OffsetRect(&expectedClip,-expected.left,-expected.top);
        const bool equal=EqualRect(&actualClip,&expectedClip)!=FALSE;
        const bool simple=regionType==SIMPLEREGION && boundsType==SIMPLEREGION;
        DeleteObject(region);
        require(simple && equal,"Actual window region does not clip to the browser viewport");
    }

    static void dom_clipping(HDESK desktop) {
        Windows windows;
        windows.register_class(L"Chrome_WidgetWin_1");
        windows.register_class(L"Chrome_RenderWidgetHostHWND");
        windows.register_class(L"TIHNT_ViewportContainer");
        const HWND browser=windows.create(L"Chrome_WidgetWin_1",nullptr,100,100,400,300);
        const HWND container=windows.create(L"TIHNT_ViewportContainer",browser,50,40,100,100);
        const HWND host=windows.create(L"Chrome_RenderWidgetHostHWND",container,-20,-30,240,180);
        activate(browser,desktop);
        OverlayWindow overlay; overlay.foreground_window_=&GetActiveWindow;
        require(overlay.create(),"DOM clipping overlay creation failed");
        OverlayGeometry geometry{}; geometry.board_w=4; geometry.board_h=3;
        geometry.rect_w=240; geometry.rect_h=180;
        std::vector<solve::Mark> marks(12,solve::Mark::Mine);
        overlay.update(marks,geometry,12);
        const HBITMAP bitmap=overlay.dib_; const HDC dc=overlay.memdc_;
        const auto originalPixels=pixels(overlay);
        auto same_surface=[&] {
            require(overlay.dib_==bitmap && overlay.memdc_==dc && pixels(overlay)==originalPixels && !overlay.surface_dirty_,
                "Clip-only change replaced, altered, or dirtied the surface");
        };
        auto clipped=[&](const RECT& expected) {
            check_presentation(overlay,host,geometry);
            require(EqualRect(&overlay.visible_rect_,&expected)!=FALSE,"DOM and HWND intersections disagree");
            same_surface();
        };
        auto mouse=[&](POINT point,bool expected) {
            require(ClientToScreen(host,&point)!=FALSE,"Cannot map DOM-clipped mouse point");
            overlay.safety_mode_=true;
            const bool down=overlay.block_mouse(WM_LBUTTONDOWN,point);
            const bool up=overlay.block_mouse(WM_LBUTTONUP,point);
            overlay.safety_mode_=false;
            require(down==expected && up==expected,"Safety disagrees with the DOM window region");
        };
        clipped({20,30,120,130});
        geometry.has_clip=true; geometry.clip_l=40.25; geometry.clip_t=20.5;
        geometry.clip_w=70.25; geometry.clip_h=100.25;
        overlay.update(marks,geometry,12); clipped({41,30,110,120});
        mouse({40,50},false); mouse({41,50},true); mouse({109,119},true); mouse({110,119},false);
        geometry.clip_l=60.75; geometry.clip_t=30.25; geometry.clip_w=70.5; geometry.clip_h=60.5;
        overlay.update(marks,geometry,12); clipped({61,31,120,90});
        mouse({60,50},false); mouse({61,50},true);
        move(container,50,40,70,70);
        mouse({100,50},false);
        overlay.tick(); clipped({61,31,90,90});
        move(container,50,40,100,100); overlay.tick(); clipped({61,31,120,90});

        geometry.clip_w=0;
        overlay.update(marks,geometry,12);
        require(!IsWindowVisible(overlay.hwnd_),"Empty DOM clip retained visible hints");
        mouse({80,50},false);
        for(int i=0;i<3;++i) overlay.tick();
        same_surface();
        geometry.clip_w=70.5;
        overlay.update(marks,geometry,12); clipped({61,31,120,90});
        geometry.has_clip=false;
        overlay.update(marks,geometry,12); clipped({20,30,120,130});

        geometry.has_clip=true; geometry.clip_l=200; geometry.clip_t=10;
        overlay.update(marks,geometry,12); overlay.tick();
        require(!IsWindowVisible(overlay.hwnd_),"Disjoint DOM and native clips retained hints");
        same_surface();
        geometry.clip_l=40.25; geometry.clip_t=20.5; geometry.clip_w=0;
        overlay.update(marks,geometry,12);
        marks[5]=solve::Mark::Safe;
        overlay.update(marks,geometry,11);
        require(overlay.surface_dirty_ && pixels(overlay)==originalPixels,"Hidden content update painted or lost pending changes");
        geometry.clip_w=70.25;
        overlay.update(marks,geometry,11);
        check_presentation(overlay,host,geometry);
        require(overlay.painted_marks_==marks && !overlay.surface_dirty_ && pixels(overlay)!=originalPixels,
            "Clip-only restoration reused stale marks from before the board was hidden");

        geometry.dpr=1.5; geometry.vv_scale=1.25;
        geometry.rect_l=-7.2; geometry.rect_t=3; geometry.rect_w=100; geometry.rect_h=80;
        geometry.vv_x=2.4; geometry.vv_y=7;
        geometry.clip_l=20.25; geometry.clip_t=25.5; geometry.clip_w=50.1; geometry.clip_h=30.25;
        overlay.update(marks,geometry,11);
        check_presentation(overlay,host,geometry);
        const RECT fractionalClip{52,43,138,99};
        require(EqualRect(&overlay.visible_rect_,&fractionalClip)!=FALSE,
            "Fractional DOM clip did not combine viewport offset, DPR, and zoom");
        geometry.dpr=geometry.vv_scale=1; geometry.vv_x=geometry.vv_y=0;
        geometry.rect_l=0.4999; geometry.rect_t=0; geometry.rect_w=239.4999; geometry.rect_h=180;
        overlay.update(marks,geometry,11); check_presentation(overlay,host,geometry);
        geometry.rect_l=0.5001;
        overlay.update(marks,geometry,11); check_presentation(overlay,host,geometry);
        geometry.rect_w=239.5001;
        overlay.update(marks,geometry,11); check_presentation(overlay,host,geometry);
        overlay.destroy(); check_isolated(desktop);
    }

    static void ancestor_clipping(HDESK desktop) {
        Windows windows;
        windows.register_class(L"Chrome_WidgetWin_1");
        windows.register_class(L"Chrome_RenderWidgetHostHWND");
        windows.register_class(L"TIHNT_ViewportContainer");
        const HWND browser=windows.create(L"Chrome_WidgetWin_1",nullptr,100,100,400,300);
        const HWND container=windows.create(L"TIHNT_ViewportContainer",browser,50,40,100,100);
        const HWND host=windows.create(L"Chrome_RenderWidgetHostHWND",container,-20,-30,240,180);
        activate(browser,desktop);
        OverlayWindow overlay; overlay.foreground_window_=&GetActiveWindow;
        require(overlay.create(),"Ancestor clipping overlay creation failed");
        OverlayGeometry geometry{}; geometry.board_w=4; geometry.board_h=3;
        geometry.rect_w=240; geometry.rect_h=180;
        overlay.update(std::vector<solve::Mark>(12,solve::Mark::Mine),geometry,12);
        check_presentation(overlay,host,geometry);
        const HBITMAP bitmap=overlay.dib_; const HDC dc=overlay.memdc_;
        const auto originalPixels=pixels(overlay);
        auto clipped=[&](const RECT& expected) {
            require(EqualRect(&overlay.visible_rect_,&expected)!=FALSE,"Ancestor clipping has incorrect edges");
            check_presentation(overlay,host,geometry);
            require(overlay.dib_==bitmap && overlay.memdc_==dc && pixels(overlay)==originalPixels,
                "Ancestor clipping replaced or repainted an unchanged surface");
        };
        clipped({20,30,120,130});
        auto mouse=[&](POINT point,bool expected) {
            require(ClientToScreen(host,&point)!=FALSE,"Cannot map ancestor-clipped mouse point");
            overlay.safety_mode_=true;
            const bool blocked=overlay.block_mouse(WM_LBUTTONDOWN,point);
            const bool released=overlay.block_mouse(WM_LBUTTONUP,point);
            overlay.safety_mode_=false;
            require(blocked==expected && released==expected,"Mouse safety disagrees with ancestor clipping");
        };
        mouse({2,2},false); mouse({50,50},true);
        POINT before{},after{}; ClientToScreen(host,&before);
        move(container,50,40,60,60);
        ClientToScreen(host,&after);
        require(before.x==after.x && before.y==after.y,"Parent resize unexpectedly moved the host");
        mouse({100,100},false);
        wait_for_window_event([&] { return overlay.visible_rect_.right==80 && overlay.visible_rect_.bottom==90; },
            "Ancestor-only resize waited for a periodic overlay tick");
        clipped({20,30,80,90}); mouse({50,50},true);
        move(container,50,40,0,60);
        wait_for_window_event([&] { return !IsWindowVisible(overlay.hwnd_); },
            "Empty ancestor did not hide the overlay");
        move(container,50,40,100,100);
        wait_for_window_event([&] { return IsWindowVisible(overlay.hwnd_)!=FALSE; },
            "Expanded ancestor did not restore the overlay");
        clipped({20,30,120,130});

        const HWND nested=windows.create(L"TIHNT_ViewportContainer",container,-10,20,120,90);
        require(SetParent(host,nested)==container,"Cannot nest the render host");
        move(host,-10,-50,240,180); overlay.tick();
        clipped({20,50,120,130});
        move(nested,20,-10,120,90);
        wait_for_window_event([&] { return overlay.visible_rect_.right==90 && overlay.visible_rect_.top==60; },
            "Nested ancestor movement did not update clipping");
        clipped({10,60,90,140});
        require(SetParent(host,container)==nested,"Cannot restore the render host parent");
        ShowWindow(nested,SW_HIDE);
        move(container,350,250,200,200); move(host,0,0,240,180); overlay.tick();
        clipped({0,0,50,50});
        mouse({25,25},true); mouse({60,60},false);
        move(container,-40,-30,200,200); overlay.tick();
        clipped({40,30,200,180});
        move(container,500,400,200,200);
        wait_for_window_event([&] { return !IsWindowVisible(overlay.hwnd_); },
            "Host outside the top client retained an overlay");
        require(IsWindowVisible(host)!=FALSE,"Fully clipped child lost its visible style");
        move(container,50,40,100,100); move(host,-20,-30,240,180);
        wait_for_window_event([&] { return IsWindowVisible(overlay.hwnd_)!=FALSE; },
            "Host returning inside the browser stayed hidden");
        clipped({20,30,120,130});

        // Border pixels are not part of the ancestor's client viewport.
        SetWindowLongPtrW(container,GWL_STYLE,GetWindowLongPtrW(container,GWL_STYLE)|WS_BORDER);
        require(SetWindowPos(container,nullptr,0,0,0,0,
            SWP_NOMOVE|SWP_NOSIZE|SWP_NOZORDER|SWP_NOACTIVATE|SWP_FRAMECHANGED)!=FALSE,
            "Cannot add an ancestor border");
        overlay.tick(); check_presentation(overlay,host,geometry);
        const HWND otherBrowser=windows.create(L"Chrome_WidgetWin_1",nullptr,600,100,200,200);
        activate(browser,desktop);
        require(SetParent(host,otherBrowser)==container,"Cannot reparent the cached host");
        overlay.tick();
        require(!overlay.cached_host_ && !IsWindowVisible(overlay.hwnd_),
            "Reparented host remained attached to the former browser");
        require(SetParent(host,container)==otherBrowser,"Cannot return the reparented host");
        move(host,-20,-30,240,180); overlay.tick();
        check_presentation(overlay,host,geometry);

        drain_messages();
        require(UnhookWinEvent(overlay.location_hook_)!=FALSE,"Cannot disable ancestor movement notifications");
        overlay.location_hook_=nullptr;
        const RECT oldClip=overlay.visible_rect_;
        move(container,50,40,70,70); drain_messages();
        require(EqualRect(&oldClip,&overlay.visible_rect_)!=FALSE,"Ancestor fallback still received notifications");
        overlay.tick(); check_presentation(overlay,host,geometry);
        require(!EqualRect(&oldClip,&overlay.visible_rect_),"Polling missed ancestor-only clipping changes");
        overlay.destroy(); check_isolated(desktop);
    }

    static void presentation(HDESK desktop) {
        Windows windows;
        windows.register_class(L"Chrome_WidgetWin_1");
        windows.register_class(L"Chrome_RenderWidgetHostHWND");
        windows.register_class(L"TIHNT_ViewportContainer");
        windows.register_class(L"TIHNT_UnrelatedWindow");
        const HWND browser=windows.create(L"Chrome_WidgetWin_1",nullptr,120,100,600,500);
        windows.create(L"Chrome_RenderWidgetHostHWND",browser,0,0,20,20,false);
        const HWND container=windows.create(L"TIHNT_ViewportContainer",browser,0,0,600,500);
        HWND host=windows.create(L"Chrome_RenderWidgetHostHWND",container,30,40,320,200);
        const HWND unrelated=windows.create(L"TIHNT_UnrelatedWindow",nullptr,150,150,200,100);
        activate(browser,desktop);

        OverlayWindow overlay;
        overlay.foreground_window_=&GetActiveWindow;
        require(overlay.create(),"Isolated overlay creation failed");
        require(!IsWindowVisible(overlay.hwnd_),"New overlay appeared before receiving board geometry");
        OverlayGeometry geometry;
        geometry.board_w=4; geometry.board_h=3;
        geometry.rect_l=-20; geometry.rect_t=-15;
        geometry.rect_w=240; geometry.rect_h=144;
        std::vector<solve::Mark> marks(12,solve::Mark::Safe);
        marks[5]=solve::Mark::Mine;
        overlay.update(marks,geometry,1);
        check_presentation(overlay,host,geometry);
        require(overlay.foreground_hook_ && overlay.location_hook_ && overlay.visibility_hook_,
            "Window event hooks were not installed");
        require(GetActiveWindow()==browser,"Overlay presentation stole activation");
        require(overlay.cached_host_==host,"Hidden render host was selected");

        const HBITMAP originalBitmap=overlay.dib_;
        const HDC originalDc=overlay.memdc_;
        const auto originalPixels=pixels(overlay);
        auto same_surface=[&] {
            require(overlay.dib_==originalBitmap && overlay.memdc_==originalDc && pixels(overlay)==originalPixels,
                "Position-only presentation replaced or altered the bitmap");
        };
        auto board_point=[&](int x,int y) {
            POINT point{x,y};
            require(ClientToScreen(host,&point)!=FALSE,"Cannot map mouse test point");
            point.x+=static_cast<LONG>(std::round((geometry.rect_l-geometry.vv_x)*geometry.dpr*geometry.vv_scale));
            point.y+=static_cast<LONG>(std::round((geometry.rect_t-geometry.vv_y)*geometry.dpr*geometry.vv_scale));
            return point;
        };
        auto mouse_down=[&](POINT point,bool expected,const char* message) {
            // Exercise the hook's policy directly without installing a hook or injecting input.
            overlay.safety_mode_=true;
            const bool blocked=overlay.block_mouse(WM_LBUTTONDOWN,point);
            const bool released=overlay.block_mouse(WM_LBUTTONUP,point);
            overlay.safety_mode_=false;
            require(blocked==expected && released==expected,message);
        };
        geometry.rect_l=10; geometry.rect_t=8;
        overlay.update(marks,geometry,1);
        check_presentation(overlay,host,geometry); same_surface();
        const POINT oldMine=board_point(90,72);
        require(WindowFromPoint(oldMine)==host,"Transparent overlay changed the mouse input target");
        mouse_down(oldMine,true,"Known mine was not blocked before browser movement");
        const HWND popup=windows.create(L"TIHNT_UnrelatedWindow",nullptr,oldMine.x-10,oldMine.y-10,20,20,false);
        SetWindowLongPtrW(popup,GWL_EXSTYLE,GetWindowLongPtrW(popup,GWL_EXSTYLE)|WS_EX_NOACTIVATE|WS_EX_TOPMOST);
        require(SetWindowPos(popup,HWND_TOPMOST,0,0,0,0,SWP_NOMOVE|SWP_NOSIZE|SWP_NOACTIVATE|SWP_SHOWWINDOW)!=FALSE,
            "Cannot present isolated popup");
        require(GetActiveWindow()==browser && WindowFromPoint(oldMine)==popup,"Popup did not cover the browser without activation");
        mouse_down(oldMine,false,"Safety blocked a native popup above a mine");
        ShowWindow(popup,SW_HIDE);
        mouse_down(oldMine,true,"Dismissed popup suppressed board safety");
        move(browser,210,170,600,500);
        mouse_down(board_point(90,72),true,"Browser movement before tick permitted a mine click");
        mouse_down(oldMine,false,"Browser movement before tick blocked a formerly unsafe screen point");
        overlay.tick();
        check_presentation(overlay,host,geometry); same_surface();

        move(browser,270,190,600,500);
        wait_for_window_event([&] {
            POINT origin{}; RECT actual{};
            ClientToScreen(host,&origin); GetWindowRect(overlay.hwnd_,&actual);
            return actual.left==origin.x+10 && actual.top==origin.y+8;
        },"Browser movement waited for a periodic overlay tick");
        check_presentation(overlay,host,geometry); same_surface();
        for(int i=0;i<60;++i) move(browser,270+i,190+i%7,600,500);
        wait_for_window_event([&] {
            POINT origin{}; RECT actual{};
            ClientToScreen(host,&origin); GetWindowRect(overlay.hwnd_,&actual);
            return actual.left==origin.x+10 && actual.top==origin.y+8;
        },"Burst window movement did not settle at its latest position");
        check_presentation(overlay,host,geometry); same_surface();
        move(container,8,6,600,500);
        wait_for_window_event([&] {
            POINT origin{}; RECT actual{};
            ClientToScreen(host,&origin); GetWindowRect(overlay.hwnd_,&actual);
            return actual.left==origin.x+10 && actual.top==origin.y+8;
        },"Intermediate host ancestor movement waited for a periodic overlay tick");
        check_presentation(overlay,host,geometry); same_surface();
        move(host,30,40,80,60);
        mouse_down(board_point(90,72),false,"Viewport resize before tick blocked a click outside the browser");
        move(host,50,60,80,60);
        overlay.tick();
        check_presentation(overlay,host,geometry); same_surface();
        move(host,30,40,320,200);
        mouse_down(board_point(90,72),true,"Viewport expansion before tick permitted a newly visible mine");
        overlay.tick();
        check_presentation(overlay,host,geometry); same_surface();

        // Repeated region replacement must preserve the same actual
        // layered-window surface allocation and pixels.
        for(int i=0;i<140;++i) {
            geometry.rect_l=(i%7)*13-40;
            geometry.rect_t=(i%5)*9-24;
            overlay.update(marks,geometry,1);
            check_presentation(overlay,host,geometry); same_surface();
            drain_messages();
        }

        activate(unrelated,desktop);
        mouse_down(board_point(90,72),false,"Inactive browser blocked mouse input before tick");
        // The isolated desktop substitutes GetActiveWindow for foreground lookup.
        // Send its corresponding accessibility notification without switching desktops.
        NotifyWinEvent(EVENT_SYSTEM_FOREGROUND,unrelated,OBJID_WINDOW,CHILDID_SELF);
        wait_for_window_event([&] { return !IsWindowVisible(overlay.hwnd_); },
            "Foreground notification did not hide the old overlay");
        require(!IsWindowVisible(overlay.hwnd_),"Overlay stayed visible over an unrelated application");
        activate(browser,desktop);
        NotifyWinEvent(EVENT_SYSTEM_FOREGROUND,browser,OBJID_WINDOW,CHILDID_SELF);
        wait_for_window_event([&] { return IsWindowVisible(overlay.hwnd_)!=FALSE; },
            "Foreground notification did not restore the overlay");
        check_presentation(overlay,host,geometry);
        overlay.safety_mode_=true;
        const POINT mine=board_point(90,72);
        require(!overlay.block_mouse(WM_RBUTTONDOWN,mine) && !overlay.block_mouse(WM_RBUTTONUP,mine) &&
            !overlay.block_mouse(WM_MOUSEWHEEL,mine) && !overlay.block_mouse(WM_MOUSEMOVE,mine),
            "Left-button safety swallowed a right click, scroll, or movement");
        require(overlay.block_mouse(WM_LBUTTONDOWN,mine),"Known mine did not begin a blocked gesture");
        activate(unrelated,desktop);
        require(overlay.block_mouse(WM_LBUTTONUP,{0,0}),"Focus change leaked the release of a blocked gesture");
        require(!overlay.block_mouse(WM_LBUTTONUP,mine),"Blocked release was consumed more than once");
        overlay.safety_mode_=false;
        activate(browser,desktop);
        overlay.set_visible(false); overlay.tick();
        require(!IsWindowVisible(overlay.hwnd_),"Visibility preference was ignored by tick");
        overlay.set_visible(true); check_presentation(overlay,host,geometry);
        overlay.set_target_pid(GetCurrentProcessId()+1);
        require(!IsWindowVisible(overlay.hwnd_),"Overlay ignored a target-process mismatch");
        overlay.set_target_pid(GetCurrentProcessId()); check_presentation(overlay,host,geometry);

        ShowWindow(host,SW_HIDE);
        wait_for_window_event([&] { return !IsWindowVisible(overlay.hwnd_); },
            "Hidden host waited for a periodic overlay tick");
        require(!IsWindowVisible(overlay.hwnd_),"Hidden render host retained an overlay");
        ShowWindow(host,SW_SHOWNOACTIVATE);
        wait_for_window_event([&] { return IsWindowVisible(overlay.hwnd_)!=FALSE; },
            "Shown host waited for a periodic overlay tick");
        check_presentation(overlay,host,geometry);
        ShowWindow(container,SW_HIDE);
        wait_for_window_event([&] { return !IsWindowVisible(overlay.hwnd_); },
            "Hidden host ancestor waited for a periodic overlay tick");
        ShowWindow(container,SW_SHOWNOACTIVATE);
        wait_for_window_event([&] { return IsWindowVisible(overlay.hwnd_)!=FALSE; },
            "Shown host ancestor waited for a periodic overlay tick");
        check_presentation(overlay,host,geometry);
        ShowWindow(browser,SW_SHOWMINNOACTIVE); overlay.tick();
        require(!IsWindowVisible(overlay.hwnd_),"Minimized browser retained an overlay");
        ShowWindow(browser,SW_SHOWNOACTIVATE); activate(browser,desktop);
        overlay.tick(); check_presentation(overlay,host,geometry);

        geometry.dpr=1.5; geometry.vv_scale=1.25;
        geometry.rect_l=-7.2; geometry.rect_t=3;
        geometry.vv_x=2.4; geometry.vv_y=7;
        overlay.update(marks,geometry,1); check_presentation(overlay,host,geometry);
        geometry.rect_l=10000;
        overlay.update(marks,geometry,1); overlay.tick();
        require(!IsWindowVisible(overlay.hwnd_),"Fully clipped board retained an overlay");
        geometry.rect_l=0;
        overlay.update(marks,geometry,1); check_presentation(overlay,host,geometry);
        require(DestroyWindow(host)!=FALSE,"Synthetic render host destruction failed");
        wait_for_window_event([&] { return !IsWindowVisible(overlay.hwnd_); },
            "Destroyed host waited for a periodic overlay tick");
        require(!IsWindowVisible(overlay.hwnd_),"Destroyed render host retained an overlay");
        host=windows.create(L"Chrome_RenderWidgetHostHWND",browser,12,18,280,180);
        wait_for_window_event([&] { return IsWindowVisible(overlay.hwnd_)!=FALSE; },
            "Replacement host waited for a periodic overlay tick");
        check_presentation(overlay,host,geometry);

        drain_messages();
        require(UnhookWinEvent(overlay.location_hook_)!=FALSE,"Cannot simulate unavailable movement notifications");
        overlay.location_hook_=nullptr;
        RECT beforeFallback{},afterFallback{};
        GetWindowRect(overlay.hwnd_,&beforeFallback);
        move(browser,420,220,600,500);
        drain_messages(); GetWindowRect(overlay.hwnd_,&afterFallback);
        require(EqualRect(&beforeFallback,&afterFallback)!=FALSE,"Movement fallback test still received notifications");
        overlay.tick(); check_presentation(overlay,host,geometry);
        overlay.watched_pid_=0; overlay.watch_browser(browser);

        const HWND overlayHandle=overlay.hwnd_;
        const HBITMAP finalBitmap=overlay.dib_;
        const HDC finalDc=overlay.memdc_;
        const auto foregroundHook=overlay.foreground_hook_,locationHook=overlay.location_hook_,visibilityHook=overlay.visibility_hook_;
        overlay.destroy();
        require(!IsWindow(overlayHandle) && !overlay.hwnd_ && !overlay.dib_ && !overlay.memdc_ && !overlay.bits_,
            "Overlay destruction retained a window or drawing resources");
        require(!overlay.foreground_hook_ && !overlay.location_hook_ && !overlay.visibility_hook_ &&
            !UnhookWinEvent(foregroundHook) && !UnhookWinEvent(locationHook) && !UnhookWinEvent(visibilityHook),
            "Overlay destruction retained its window event hooks");
        drain_messages();
        BITMAP bitmap{};
        const int bitmapSize=GetObjectW(finalBitmap,sizeof(bitmap),&bitmap);
        const HGDIOBJ selected=GetCurrentObject(finalDc,OBJ_BITMAP);
        if(bitmapSize || selected) std::cerr<<"Destroyed bitmap size="<<bitmapSize<<" selected="<<selected<<'\n';
        require(!bitmapSize && !selected,
            "Overlay destruction lost handles without releasing GDI resources");
        check_isolated(desktop);
    }
};

int main() {
    HDESK desktop=nullptr;
    std::exception_ptr failure;
    std::thread worker([&] {
        try {
            const std::wstring name=L"TIHNT_overlay_test_"+std::to_wstring(GetCurrentProcessId())+L"_"+std::to_wstring(GetTickCount64());
            // Request no DESKTOP_SWITCHDESKTOP right. All HWNDs are created on
            // this fresh thread after moving it off the user's input desktop.
            desktop=CreateDesktopW(name.c_str(),nullptr,nullptr,0,
                DESKTOP_CREATEWINDOW|DESKTOP_READOBJECTS|DESKTOP_WRITEOBJECTS,nullptr);
            require(desktop!=nullptr,"Cannot create an isolated test desktop");
            require(SetThreadDesktop(desktop)!=FALSE,"Cannot bind the test thread to its isolated desktop");
            check_isolated(desktop);
            // Reparenting can start IME helper threads that retain this desktop.
            // No input method is needed on the isolated test thread.
            require(ImmDisableIME(GetCurrentThreadId())!=FALSE,"Cannot disable the isolated thread's input method");
            SetThreadDpiAwarenessContext(DPI_AWARENESS_CONTEXT_PER_MONITOR_AWARE_V2);
            OverlayTestAccess::native_region_coordinates(desktop);
            OverlayTestAccess::native_regions(desktop);
            OverlayTestAccess::presentation(desktop);
            OverlayTestAccess::ancestor_clipping(desktop);
            OverlayTestAccess::dom_clipping(desktop);
            check_isolated(desktop);
        } catch(...) { failure=std::current_exception(); }
    });
    worker.join();
    // Closing from the parent after the UI thread exits releases the thread's
    // desktop association as well as all window and DWM references.
    if(desktop && !CloseDesktop(desktop) && !failure)
        failure=std::make_exception_ptr(std::runtime_error("Isolated desktop cleanup failed (Win32 error "+
            std::to_string(GetLastError())+")"));
    try {
        if(failure) std::rethrow_exception(failure);
        std::cout<<"Isolated Win32 overlay presentation checks passed\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
