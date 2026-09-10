#include "overlay.hpp"
#include <algorithm>
#include <iostream>
#include <stdexcept>

struct OverlayTestAccess {
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

int main() {
    try {
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
