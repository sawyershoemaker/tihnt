#include "board.hpp"
#include <algorithm>
#include <stdexcept>
#include <utility>

namespace game {

bool valid_dimensions(int w, int h) {
    return (w == 0 && h == 0) || (w > 0 && h > 0 && w <= MaxBoardDimension &&
        h <= MaxBoardDimension && w <= MaxBoardCells / h);
}

bool valid_cell_state(int state) {
    return state == 0 || state == 1 || state == 2 || (state >= 10 && state <= 18);
}

static void validate(const std::vector<CellState>& cells, int w, int h) {
    if (!valid_dimensions(w, h) || cells.size() != static_cast<size_t>(w) * h ||
        !std::all_of(cells.begin(), cells.end(), [](CellState s) { return valid_cell_state(static_cast<int>(s)); })) {
        throw std::invalid_argument("Invalid board dimensions or cells");
    }
}

Board::Board() {}

void Board::resize(int w, int h){
    if (!valid_dimensions(w, h)) throw std::invalid_argument("Invalid board dimensions");
    std::vector<CellState> cells(static_cast<size_t>(w) * h, CellState::Unknown);
    cells_ = std::move(cells); w_ = w; h_ = h;
}

int Board::width() const { return w_; }
int Board::height() const { return h_; }

int Board::index(int x, int y) const { return y*w_ + x; }

CellState Board::at(int x, int y) const {
    if(x<0||y<0||x>=w_||y>=h_) return CellState::Unknown;
    return cells_[index(x,y)];
}

void Board::set(int x, int y, CellState s){
    if(x<0||y<0||x>=w_||y>=h_||!valid_cell_state(static_cast<int>(s))) return;
    cells_[index(x,y)] = s;
}

void Board::apply_updates(const std::vector<CellUpdate>& updates){
    for(const auto& u : updates){ set(u.x,u.y,u.state); }
}

void Board::apply_full(const std::vector<CellState>& all, int w, int h){
    validate(all, w, h);
    cells_ = all; w_ = w; h_ = h;
}

void Board::apply_full(std::vector<CellState>&& all, int w, int h){
    validate(all, w, h);
    cells_ = std::move(all); w_ = w; h_ = h;
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
