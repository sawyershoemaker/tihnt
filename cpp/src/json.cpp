#include "proto.hpp"
#include "utf8.hpp"

#include <charconv>
#include <cmath>
#include <string_view>
#include <unordered_set>
#include <utility>

namespace proto {
namespace {

class Reader {
public:
    explicit Reader(std::string_view text) : text_(text) {}
    void whitespace() {
        while (pos_ < text_.size() && (text_[pos_] == ' ' || text_[pos_] == '\t' || text_[pos_] == '\r' || text_[pos_] == '\n')) ++pos_;
    }
    bool take(char c) {
        whitespace();
        if (pos_ == text_.size() || text_[pos_] != c) return false;
        ++pos_; return true;
    }
    bool finished() { whitespace(); return pos_ == text_.size(); }
    bool string(std::string& out) {
        if (!take('"')) return false;
        out.clear();
        while (pos_ < text_.size()) {
            char c = text_[pos_++];
            if (c == '"') return true;
            if (static_cast<unsigned char>(c) < 0x20) return false;
            if (c != '\\') { out.push_back(c); continue; }
            if (pos_ == text_.size()) return false;
            c = text_[pos_++];
            switch (c) {
            case '"': case '\\': case '/': out.push_back(c); break;
            case 'b': out.push_back('\b'); break;
            case 'f': out.push_back('\f'); break;
            case 'n': out.push_back('\n'); break;
            case 'r': out.push_back('\r'); break;
            case 't': out.push_back('\t'); break;
            case 'u': {
                uint32_t cp;
                if (!hex4(cp)) return false;
                if (cp >= 0xd800 && cp <= 0xdbff) {
                    if (text_.substr(pos_, 2) != "\\u") return false;
                    pos_ += 2;
                    uint32_t low;
                    if (!hex4(low) || low < 0xdc00 || low > 0xdfff) return false;
                    cp = 0x10000 + ((cp - 0xd800) << 10) + (low - 0xdc00);
                } else if (cp >= 0xdc00 && cp <= 0xdfff) return false;
                if (cp < 0x80) out.push_back(static_cast<char>(cp));
                else {
                    if (cp >= 0x10000) out.push_back(static_cast<char>(0xf0 | (cp >> 18)));
                    if (cp >= 0x800) out.push_back(static_cast<char>((cp >= 0x10000 ? 0x80 : 0xe0) | ((cp >> 12) & 0x3f)));
                    out.push_back(static_cast<char>((cp >= 0x800 ? 0x80 : 0xc0) | ((cp >> 6) & 0x3f)));
                    out.push_back(static_cast<char>(0x80 | (cp & 0x3f)));
                }
                break;
            }
            default: return false;
            }
        }
        return false;
    }
    bool integer(int& out) {
        std::string_view token;
        if (!number_token(token) || token.find_first_of(".eE") != std::string_view::npos) return false;
        const auto result = std::from_chars(token.data(), token.data() + token.size(), out);
        return result.ec == std::errc{} && result.ptr == token.data() + token.size();
    }
    bool number(double& out) {
        std::string_view token;
        if (!number_token(token)) return false;
        const auto result = std::from_chars(token.data(), token.data() + token.size(), out);
        return result.ec == std::errc{} && result.ptr == token.data() + token.size() && std::isfinite(out);
    }
    template<class Field> bool object(Field field) {
        if (!take('{')) return false;
        if (take('}')) return true;
        std::unordered_set<std::string> keys;
        do {
            std::string key;
            if (!string(key) || !keys.insert(key).second || !take(':') || !field(key)) return false;
            if (take('}')) return true;
        } while (take(','));
        return false;
    }
    template<class Item> bool array(Item item) {
        if (!take('[')) return false;
        if (take(']')) return true;
        do {
            if (!item()) return false;
            if (take(']')) return true;
        } while (take(','));
        return false;
    }
    bool skip(int depth = 0) {
        if (depth >= 32) return false;
        whitespace();
        if (pos_ == text_.size()) return false;
        if (text_[pos_] == '{') return object([&](const auto&) { return skip(depth + 1); });
        if (text_[pos_] == '[') return array([&] { return skip(depth + 1); });
        if (text_[pos_] == '"') { std::string unused; return string(unused); }
        for (std::string_view literal : {"true", "false", "null"}) {
            if (text_.substr(pos_, literal.size()) == literal) { pos_ += literal.size(); return true; }
        }
        double unused; return number(unused);
    }
private:
    bool hex4(uint32_t& value) {
        value = 0;
        for (int i = 0; i < 4; ++i) {
            if (pos_ == text_.size()) return false;
            const char c = text_[pos_++];
            int digit = c >= '0' && c <= '9' ? c - '0' : c >= 'a' && c <= 'f' ? c - 'a' + 10 : c >= 'A' && c <= 'F' ? c - 'A' + 10 : -1;
            if (digit < 0) return false;
            value = value * 16 + digit;
        }
        return true;
    }
    bool number_token(std::string_view& out) {
        whitespace(); const size_t start = pos_;
        if (pos_ < text_.size() && text_[pos_] == '-') ++pos_;
        auto digit = [&] { return pos_ < text_.size() && text_[pos_] >= '0' && text_[pos_] <= '9'; };
        if (!digit()) return false;
        if (text_[pos_] == '0') { ++pos_; if (digit()) return false; }
        else while (digit()) ++pos_;
        if (pos_ < text_.size() && text_[pos_] == '.') {
            ++pos_; if (!digit()) return false; while (digit()) ++pos_;
        }
        if (pos_ < text_.size() && (text_[pos_] == 'e' || text_[pos_] == 'E')) {
            ++pos_;
            if (pos_ < text_.size() && (text_[pos_] == '+' || text_[pos_] == '-')) ++pos_;
            if (!digit()) return false;
            while (digit()) ++pos_;
        }
        out = text_.substr(start, pos_ - start); return true;
    }
    std::string_view text_;
    size_t pos_ = 0;
};

} // namespace

bool parse_message(const std::string& json, ParsedMessage& out) {
    out = {};
    if (json.size() > MaxMessageBytes || !valid_utf8(json)) return false;
    Reader r(json);
    ParsedMessage parsed;
    GeometryMsg geometry;
    std::string type;
    bool gotW = false, gotH = false, gotCells = false, gotUpdates = false, gotPid = false;
    unsigned rectKeys = 0, clipKeys = 0;
    const bool ok = r.object([&](const std::string& key) {
        if (key == "type") return r.string(type);
        if (key == "w") { gotW = true; return r.integer(parsed.full.w); }
        if (key == "h") { gotH = true; return r.integer(parsed.full.h); }
        if (key == "pid") { gotPid = true; return r.integer(parsed.bind.pid) && parsed.bind.pid > 0; }
        if (key == "cells") {
            gotCells = true;
            return r.array([&] {
                int state;
                if (parsed.full.cells.size() >= game::MaxBoardCells || !r.integer(state) || !game::valid_cell_state(state)) return false;
                parsed.full.cells.push_back(static_cast<game::CellState>(state)); return true;
            });
        }
        if (key == "updates") {
            gotUpdates = true;
            return r.array([&] {
                if (parsed.delta.updates.size() >= game::MaxBoardCells) return false;
                game::CellUpdate update{};
                unsigned fields = 0;
                if (!r.object([&](const std::string& k) {
                    if (k == "x") { fields |= 1; return r.integer(update.x); }
                    if (k == "y") { fields |= 2; return r.integer(update.y); }
                    if (k == "s") {
                        int state; fields |= 4;
                        if (!r.integer(state) || !game::valid_cell_state(state)) return false;
                        update.state = static_cast<game::CellState>(state); return true;
                    }
                    return r.skip();
                }) || fields != 7 || update.x < 0 || update.y < 0 || update.x >= game::MaxBoardDimension || update.y >= game::MaxBoardDimension) return false;
                parsed.delta.updates.push_back(update); return true;
            });
        }
        if (key == "mines_total") {
            geometry.has_mines_total = true;
            return r.integer(geometry.mines_total) && geometry.mines_total >= -1 && geometry.mines_total <= game::MaxBoardCells;
        }
        if (key == "cell_px") return r.integer(geometry.cell_px) && geometry.cell_px >= 0 && geometry.cell_px <= 16384;
        if (key == "ox") return r.integer(geometry.origin_x);
        if (key == "oy") return r.integer(geometry.origin_y);
        double* value = nullptr;
        if (key == "rect_l") { value = &geometry.rect_l; rectKeys |= 1; }
        else if (key == "rect_t") { value = &geometry.rect_t; rectKeys |= 2; }
        else if (key == "rect_w") { value = &geometry.rect_w; rectKeys |= 4; }
        else if (key == "rect_h") { value = &geometry.rect_h; rectKeys |= 8; }
        else if (key == "clip_l") { value = &geometry.clip_l; clipKeys |= 1; }
        else if (key == "clip_t") { value = &geometry.clip_t; clipKeys |= 2; }
        else if (key == "clip_w") { value = &geometry.clip_w; clipKeys |= 4; }
        else if (key == "clip_h") { value = &geometry.clip_h; clipKeys |= 8; }
        else if (key == "vv_x") value = &geometry.vv_x;
        else if (key == "vv_y") value = &geometry.vv_y;
        else if (key == "vv_scale") value = &geometry.vv_scale;
        else if (key == "dpr") value = &geometry.dpr;
        if (!value) return r.skip();
        return r.number(*value) && std::abs(*value) <= 1000000;
    });
    if (!ok || !r.finished() || (rectKeys != 0 && rectKeys != 15) ||
        (clipKeys != 0 && (clipKeys != 15 || rectKeys != 15)) ||
        geometry.rect_w < 0 || geometry.rect_h < 0 || geometry.rect_w > 16384 || geometry.rect_h > 16384 ||
        geometry.clip_w < 0 || geometry.clip_h < 0 || geometry.clip_w > 16384 || geometry.clip_h > 16384 ||
        geometry.vv_scale <= 0 || geometry.vv_scale > 16 || geometry.dpr <= 0 || geometry.dpr > 16) return false;
    geometry.has_geometry = rectKeys == 15;
    geometry.has_clip = clipKeys == 15;
    if (type == "full") {
        if (!gotW || !gotH || !gotCells || !game::valid_dimensions(parsed.full.w, parsed.full.h) ||
            parsed.full.cells.size() != static_cast<size_t>(parsed.full.w) * parsed.full.h) return false;
        parsed.type = MsgType::Full;
        static_cast<GeometryMsg&>(parsed.full) = geometry;
    } else if (type == "delta") {
        if (!gotUpdates) return false;
        parsed.type = MsgType::Delta;
        static_cast<GeometryMsg&>(parsed.delta) = geometry;
    } else if (type == "bind") {
        if (!gotPid) return false;
        parsed.type = MsgType::Bind;
    } else return false;
    out = std::move(parsed);
    return true;
}

} // namespace proto
