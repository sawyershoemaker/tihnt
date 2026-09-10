#pragma once
#include <cstdint>
#include <string_view>

inline bool valid_utf8(std::string_view s) {
    for (size_t i = 0; i < s.size();) {
        const auto c = static_cast<uint8_t>(s[i++]);
        if (c < 0x80) continue;
        int extra = 0;
        uint32_t value = 0, minimum = 0;
        if (c >= 0xc2 && c <= 0xdf) { extra = 1; value = c & 0x1f; minimum = 0x80; }
        else if (c >= 0xe0 && c <= 0xef) { extra = 2; value = c & 0x0f; minimum = 0x800; }
        else if (c >= 0xf0 && c <= 0xf4) { extra = 3; value = c & 7; minimum = 0x10000; }
        else return false;
        if (s.size() - i < static_cast<size_t>(extra)) return false;
        while (extra--) {
            const auto next = static_cast<uint8_t>(s[i++]);
            if ((next & 0xc0) != 0x80) return false;
            value = (value << 6) | (next & 0x3f);
        }
        if (value < minimum || value > 0x10ffff || (value >= 0xd800 && value <= 0xdfff)) return false;
    }
    return true;
}
