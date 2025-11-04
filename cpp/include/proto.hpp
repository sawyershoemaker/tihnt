#pragma once

#include <string>
#include <vector>
#include "board.hpp"

namespace proto {

struct FullMsg {
    int w=0;
    int h=0;
    std::vector<game::CellState> cells;
    int cell_px=0;
    int origin_x=0;
    int origin_y=0;
	int mines_total=-1;
    double rect_l=0.0;
    double rect_t=0.0;
    double rect_w=0.0;
    double rect_h=0.0;
    double vv_x=0.0;
    double vv_y=0.0;
    double vv_scale=1.0;
    double dpr=1.0;
};
struct DeltaMsg {
    std::vector<game::CellUpdate> updates;
    int cell_px=0;
    int origin_x=0;
    int origin_y=0;
	int mines_total=-1;
    double rect_l=0.0;
    double rect_t=0.0;
    double rect_w=0.0;
    double rect_h=0.0;
    double vv_x=0.0;
    double vv_y=0.0;
    double vv_scale=1.0;
    double dpr=1.0;
};

struct BindMsg {
    int pid = 0;
};

enum class MsgType { Full, Delta, Bind, Unknown };

struct ParsedMessage {
    MsgType type = MsgType::Unknown;
    FullMsg full;
    DeltaMsg delta;
    BindMsg bind;
};

bool parse_message(const std::string& json, ParsedMessage& out);

}
