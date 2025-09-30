#pragma once

#include <cstdint>
#include <vector>

namespace game {

enum class CellState : uint8_t {
    Unknown = 0,
    Clear = 1,
    Mine = 2,
    Number0 = 10,
    Number1,
    Number2,
    Number3,
    Number4,
    Number5,
    Number6,
    Number7,
    Number8
};

struct CellUpdate { int x; int y; CellState state; };

class Board {
public:
    Board();
    void resize(int w, int h);
    int width() const; 
    int height() const;

    CellState at(int x, int y) const;
    void set(int x, int y, CellState s);

    void apply_updates(const std::vector<CellUpdate>& updates);
    void apply_full(const std::vector<CellState>& all, int w, int h);

    std::vector<CellUpdate> diff(const Board& other) const;

    const std::vector<CellState>& data() const { return cells_; }

private:
    int w_ = 0;
    int h_ = 0;
    std::vector<CellState> cells_;
    int index(int x, int y) const;
};

}
