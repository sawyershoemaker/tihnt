#include "board.hpp"

namespace game {

Board::Board() {}

void Board::resize(int w, int h){
    w_ = w; h_ = h; cells_.assign(w_*h_, CellState::Unknown);
}

int Board::width() const { return w_; }
int Board::height() const { return h_; }

int Board::index(int x, int y) const { return y*w_ + x; }

CellState Board::at(int x, int y) const {
    if(x<0||y<0||x>=w_||y>=h_) return CellState::Unknown;
    return cells_[index(x,y)];
}

void Board::set(int x, int y, CellState s){
    if(x<0||y<0||x>=w_||y>=h_) return;
    cells_[index(x,y)] = s;
}

void Board::apply_updates(const std::vector<CellUpdate>& updates){
    for(const auto& u : updates){ set(u.x,u.y,u.state); }
}

void Board::apply_full(const std::vector<CellState>& all, int w, int h){
    w_ = w; h_ = h; cells_ = all; if((int)cells_.size()!=w_*h_) cells_.assign(w_*h_, CellState::Unknown);
}

std::vector<CellUpdate> Board::diff(const Board& other) const {
    std::vector<CellUpdate> d;
    if(w_!=other.w_||h_!=other.h_) return d;
    for(int y=0;y<h_;++y){
        for(int x=0;x<w_;++x){
            auto a = at(x,y); auto b = other.at(x,y);
            if(a!=b) d.push_back(CellUpdate{x,y,a});
        }
    }
    return d;
}

}
