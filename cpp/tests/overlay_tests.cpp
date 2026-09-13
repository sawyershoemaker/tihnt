#include "overlay.hpp"
#include <algorithm>
#include <chrono>
#include <cmath>
#include <iostream>
#include <random>
#include <stdexcept>

struct OverlayTestAccess {
    static void dom_clipping() {
        OverlayWindow overlay;
        auto& g=overlay.geom_;
        g.board_w=3; g.board_h=2;
        g.rect_l=-2.25; g.rect_t=3.5; g.rect_w=17.5; g.rect_h=11.25;
        const double clips[][4]={{-2.25,3.5,17.5,11.25},{0.2,4.75,8.5,5.6},{-100,-100,200,200},
            {2.25,4.25,0,4},{2.25,4.25,4,0},{2.25,4.25,0.1,0.1},{-100,-100,1,1},{100,100,1,1}};
        for(double dpr : {0.5,1.0,1.25,2.0}) for(double zoom : {0.5,1.0,1.5})
        for(double offset : {0.0,0.2,3.5}) for(const auto& cssClip : clips) {
            g.dpr=dpr; g.vv_scale=zoom; g.vv_x=offset; g.vv_y=offset*2;
            g.clip_l=cssClip[0]; g.clip_t=cssClip[1]; g.clip_w=cssClip[2]; g.clip_h=cssClip[3];
            RECT bounds{};
            if(!overlay.board_bounds(bounds)) throw std::runtime_error("Valid fractional board was rejected");
            const int width=bounds.right-bounds.left, height=bounds.bottom-bounds.top;
            for(const RECT viewport : {RECT{-100,-100,100,100},RECT{2,3,13,15},RECT{4,3,4,3}}) {
                g.has_clip=false;
                const RECT legacy=overlay.board_clip(bounds,viewport);
                const RECT native=OverlayWindow::visible_board_rect(bounds.left,bounds.top,width,height,viewport);
                if(!EqualRect(&legacy,&native)) throw std::runtime_error("Legacy geometry was clipped by absent DOM metadata");
                g.has_clip=true;
                const RECT clip=overlay.board_clip(bounds,viewport);
                const double scale=dpr*zoom;
                for(int y=-1;y<=height;++y) for(int x=-1;x<=width;++x) {
                    const POINT host{x+bounds.left,y+bounds.top};
                    const bool expected=x>=0 && y>=0 && x<width && y<height && PtInRect(&viewport,host) &&
                        host.x>=(g.clip_l-g.vv_x)*scale && host.y>=(g.clip_t-g.vv_y)*scale &&
                        host.x+1<=(g.clip_l+g.clip_w-g.vv_x)*scale && host.y+1<=(g.clip_t+g.clip_h-g.vv_y)*scale;
                    if(bool(PtInRect(&clip,{x,y}))!=expected)
                        throw std::runtime_error("DOM clip includes a fractional outside pixel or loses a complete visible pixel");
                }
            }
        }
    }
    static void geometry_and_hit_testing() {
        for(const RECT viewport : {RECT{0,0,13,7},RECT{3,2,13,7},RECT{4,3,4,3}})
        for(int x : {-1000000,-13,-3,0,3,13,1000000}) for(int y : {-1000000,-13,-3,0,3,13,1000000}) {
            const RECT clip=OverlayWindow::visible_board_rect(x,y,17,11,viewport);
            for(int py=-1;py<=11;++py) for(int px=-1;px<=17;++px) {
                const bool expected=px>=0 && py>=0 && px<17 && py<11 &&
                    px+x>=viewport.left && py+y>=viewport.top && px+x<viewport.right && py+y<viewport.bottom;
                if(bool(PtInRect(&clip,{px,py}))!=expected)
                    throw std::runtime_error("Viewport clipping includes a browser control or hides a visible pixel");
            }
        }
        OverlayWindow overlay;
        std::vector<unsigned char> pixels(31*17*4);
        overlay.bits_=pixels.data(); overlay.stride_=31*4;
        overlay.surf_w_=31; overlay.surf_h_=17;
        overlay.geom_.board_w=30; overlay.geom_.board_h=16;
        overlay.marks_.assign(480,solve::Mark::Mine);
        overlay.visible_rect_={4,3,28,16}; overlay.has_clip_=true;
        for(int mark=0;mark<=7;++mark) {
            std::fill(overlay.marks_.begin(),overlay.marks_.end(),static_cast<solve::Mark>(mark));
            overlay.paint_surface();
            for(int y=-1;y<=17;++y) for(int x=-1;x<=31;++x) {
                const auto value=static_cast<solve::Mark>(mark);
                const bool expected=x>=4 && y>=3 && x<28 && y<16 && value!=solve::Mark::Safe &&
                    value!=solve::Mark::Guess && value!=solve::Mark::ChordReady;
                if(overlay.unsafe_at({x,y},overlay.visible_rect_)!=expected)
                    throw std::runtime_error("Mouse safety disagrees with visible clipped marks");
            }
        }
        overlay.marks_.clear(); overlay.paint_surface();
        if(!overlay.unsafe_at({5,5},overlay.visible_rect_)) throw std::runtime_error("Missing marks permitted an unsafe click");
        overlay.geom_={}; overlay.paint_surface();
        if(overlay.unsafe_at({5,5},overlay.visible_rect_)) throw std::runtime_error("Removed board retained old hit testing edges");
        overlay.bits_=nullptr;
    }
    static void lifecycle() {
        OverlayWindow overlay;
        if(!overlay.create()) throw std::runtime_error("Hidden overlay creation failed");
        const HWND original=overlay.hwnd_;
        if(IsWindowVisible(original) || !overlay.create() || overlay.hwnd_!=original)
            throw std::runtime_error("Repeated creation leaked or displayed an overlay window");
        overlay.destroy();
        if(IsWindow(original) || overlay.hwnd_ || overlay.safety_mode_ || overlay.surface_painted_)
            throw std::runtime_error("Destroy retained live overlay state");
        if(!overlay.create() || IsWindowVisible(overlay.hwnd_)) throw std::runtime_error("Hidden overlay recreation failed");
    }
    static void incremental() {
        OverlayWindow overlay;
        std::vector<unsigned char> pixels;
        std::mt19937 rng(53);
        for(int frame=0;frame<250;++frame) {
            if(frame%13==0) {
                overlay.surf_w_=1+rng()%600; overlay.surf_h_=1+rng()%320;
                overlay.stride_=overlay.surf_w_*4;
                pixels.assign(static_cast<size_t>(overlay.stride_)*overlay.surf_h_,0xcc);
                overlay.bits_=pixels.data(); overlay.surface_painted_=false;
                overlay.geom_.board_w=1+rng()%30; overlay.geom_.board_h=1+rng()%16;
                overlay.marks_.assign(static_cast<size_t>(overlay.geom_.board_w)*overlay.geom_.board_h,solve::Mark::None);
            }
            for(int i=0;i<1+frame%7;++i) overlay.marks_[rng()%overlay.marks_.size()]=static_cast<solve::Mark>(rng()%8);
            if(frame%17==0) overlay.mines_total_=static_cast<int>(rng()%100)-1;
            if(frame%19==0) overlay.excluded_from_capture_=!overlay.excluded_from_capture_;
            if(frame%23==0) overlay.safety_mode_=!overlay.safety_mode_;
            overlay.paint_surface();

            OverlayWindow reference;
            std::vector<unsigned char> expected(pixels.size(),0xcc);
            reference.bits_=expected.data(); reference.stride_=overlay.stride_;
            reference.surf_w_=overlay.surf_w_; reference.surf_h_=overlay.surf_h_;
            reference.geom_=overlay.geom_; reference.marks_=overlay.marks_;
            reference.mines_total_=overlay.mines_total_;
            reference.excluded_from_capture_=overlay.excluded_from_capture_;
            reference.safety_mode_=overlay.safety_mode_;
            reference.paint_surface();
            if(pixels!=expected) throw std::runtime_error("Incremental drawing differs from a fresh full paint");
            reference.bits_=nullptr;
        }
        overlay.bits_=nullptr;
    }
    static void cell_bounds() {
        for(int width : {3,7,17,31,64}) for(int height : {3,7,17,31,64}) {
            OverlayWindow overlay;
            std::vector<unsigned char> pixels(static_cast<size_t>(width)*height*4);
            overlay.bits_=pixels.data(); overlay.stride_=width*4;
            overlay.surf_w_=width; overlay.surf_h_=height;
            overlay.geom_.board_w=3; overlay.geom_.board_h=3;
            overlay.marks_.resize(9,solve::Mark::None);
            for(int cell=0;cell<9;++cell) for(int mark=1;mark<=7;++mark) {
                std::fill(overlay.marks_.begin(),overlay.marks_.end(),solve::Mark::None);
                overlay.marks_[cell]=static_cast<solve::Mark>(mark);
                overlay.paint_surface();
                const int left=overlay.cached_xEdge_[cell%3], right=overlay.cached_xEdge_[cell%3+1];
                const int top=overlay.cached_yEdge_[cell/3], bottom=overlay.cached_yEdge_[cell/3+1];
                for(int y=0;y<height;++y) for(int x=0;x<width;++x) {
                    const auto color=reinterpret_cast<const uint32_t*>(pixels.data())[y*width+x];
                    if(color && color!=0xFFFFFFFF && (x<left || x>=right || y<top || y>=bottom))
                        throw std::runtime_error("A mark painted into its neighboring cell");
                }
            }
            overlay.bits_=nullptr;
        }
    }
    static void benchmark() {
        OverlayWindow overlay;
        constexpr int width=600,height=320;
        std::vector<unsigned char> pixels(static_cast<size_t>(width)*height*4);
        overlay.bits_=pixels.data(); overlay.stride_=width*4;
        overlay.surf_w_=width; overlay.surf_h_=height;
        overlay.geom_.board_w=30; overlay.geom_.board_h=16;
        overlay.marks_.resize(480);
        for(size_t i=0;i<overlay.marks_.size();++i) overlay.marks_[i]=static_cast<solve::Mark>(i%8);
        overlay.paint_surface();
        const auto start=std::chrono::steady_clock::now();
        for(int i=0;i<5000;++i) {
            auto& mark=overlay.marks_[static_cast<size_t>(i)%overlay.marks_.size()];
            mark=static_cast<solve::Mark>((static_cast<int>(mark)+1)%8);
            overlay.paint_surface();
        }
        const auto ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
        std::cout<<"5000 single-cell raster updates: "<<ms<<" ms\n";
        overlay.bits_=nullptr;
    }
    static void paint(int width, int height, int cols, int rows) {
        OverlayWindow overlay;
        std::vector<unsigned char> pixels(static_cast<size_t>(width) * height * 4, 0xcc);
        overlay.bits_=pixels.data(); overlay.stride_=width*4;
        overlay.surf_w_=width; overlay.surf_h_=height;
        overlay.geom_.board_w=cols; overlay.geom_.board_h=rows;
        overlay.marks_.resize(cols*rows);
        for(int mark=0;mark<=7;++mark) {
            std::fill(overlay.marks_.begin(),overlay.marks_.end(),static_cast<solve::Mark>(mark));
            overlay.paint_surface(); // Exact-size heap buffer is checked by AddressSanitizer.
            if(overlay.cached_xEdge_.back()!=width || overlay.cached_yEdge_.back()!=height)
                throw std::runtime_error("Cell edges do not cover the surface");
            for(size_t i=0;i<pixels.size();i+=4) {
                if(pixels[i]>pixels[i+3] || pixels[i+1]>pixels[i+3] || pixels[i+2]>pixels[i+3])
                    throw std::runtime_error("Pixel is not premultiplied alpha");
            }
        }
        overlay.bits_=nullptr;
    }
    static void allocation() {
        OverlayWindow overlay;
        if(overlay.ensureSurface(0,1) || overlay.ensureSurface(100000,100000) || overlay.ensureSurface(8192,8192))
            throw std::runtime_error("Unbounded surface allocation");
        if(!overlay.ensureSurface(64,32) || !overlay.ensureSurface(32,64))
            throw std::runtime_error("Valid surface allocation failed");
    }
};

int main(int argc, char** argv) {
    try {
        if(argc>1 && std::string(argv[1])=="--benchmark") { OverlayTestAccess::benchmark(); return 0; }
        OverlayTestAccess::cell_bounds();
        OverlayTestAccess::incremental();
        OverlayTestAccess::geometry_and_hit_testing();
        OverlayTestAccess::dom_clipping();
        OverlayTestAccess::lifecycle();
        // Includes five-pixel cells (the old X divided by zero), and cells
        // narrower than one screen pixel (the old filled rectangles overflowed).
        for(int width=1;width<=24;++width) for(int height=1;height<=12;++height) {
            OverlayTestAccess::paint(width,height,1,1);
            OverlayTestAccess::paint(width,height,30,16);
        }
        OverlayTestAccess::allocation();
        std::cout<<"Overlay raster and allocation checks passed\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n'; return 1;
    }
}
