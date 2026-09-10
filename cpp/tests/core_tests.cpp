#include "board.hpp"
#include "proto.hpp"
#include "solver.hpp"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <future>
#include <iostream>
#include <random>
#include <stdexcept>
#include <string>

using game::CellState;
using solve::Mark;

static void require(bool ok, const std::string& message) {
    if (!ok) throw std::runtime_error(message);
}

static game::Board make_board(int w, int h, const std::vector<int>& cells) {
    std::vector<CellState> states;
    for (int s : cells) states.push_back(static_cast<CellState>(s));
    game::Board b;
    b.apply_full(std::move(states), w, h);
    return b;
}

// Independent oracle: enumerate complete boards, without using solver constraints.
static void verify_oracle(const game::Board& b, int total, bool chords, int threads = 1) {
    const int n = b.width() * b.height();
    std::vector<int> unknown;
    std::vector<int> mines(n, 0);
    for (int i = 0; i < n; ++i) {
        if (b.data()[i] == CellState::Unknown) unknown.push_back(i);
        mines[i] = b.data()[i] == CellState::Mine;
    }
    require(unknown.size() <= 20, "oracle board too large");
    std::vector<int> hits(n, 0);
    int solutions = 0;
    for (unsigned mask = 0; mask < (1u << unknown.size()); ++mask) {
        for (size_t j = 0; j < unknown.size(); ++j) mines[unknown[j]] = (mask >> j) & 1;
        if (total >= 0 && std::count(mines.begin(), mines.end(), 1) != total) continue;
        bool valid = true;
        for (int i = 0; i < n && valid; ++i) {
            int value = static_cast<int>(b.data()[i]) - 10;
            if (value < 0 || value > 8) continue;
            int count = 0;
            for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
                int x = i % b.width() + dx, y = i / b.width() + dy;
                if ((dx || dy) && x >= 0 && y >= 0 && x < b.width() && y < b.height()) count += mines[y * b.width() + x];
            }
            valid = count == value;
        }
        if (!valid) continue;
        ++solutions;
        for (int i : unknown) hits[i] += mines[i];
    }
    const auto out = solve::compute_overlay(b, total, chords, threads);
    require(out.marks.size() == size_t(n), "wrong mark dimensions");
    bool safe = false;
    for (int i : unknown) {
        const auto m = out.marks[i];
        if (m == Mark::Safe) {
            safe = true;
            require(solutions > 0 && hits[i] == 0, "unsafe Safe at " + std::to_string(i));
        }
        if (m == Mark::Mine || m == Mark::FlagForChord) {
            require(solutions > 0 && hits[i] == solutions, "false Mine/Flag at " + std::to_string(i));
        }
        if (m == Mark::Guess) require(solutions > 0 && hits[i] != solutions, "guess on proven mine");
        if (solutions > 0 && total >= 0) {
            require(out.mineProbability.size() == size_t(n), "missing probabilities");
            double expected = double(hits[i]) / solutions;
            double actual = out.mineProbability[i];
            require(std::isfinite(actual) && std::abs(actual - expected) < 1e-9,
                "probability at " + std::to_string(i) + ": " + std::to_string(actual) + " expected " + std::to_string(expected));
        }
    }
    require(out.hasGuaranteedSafe == safe, "safe summary mismatch");
    if (safe) require(std::find(out.marks.begin(), out.marks.end(), Mark::Guess) == out.marks.end(), "guess alongside a safe move");
    for (int i=0;i<n;++i) if(out.marks[i]==Mark::Chord || out.marks[i]==Mark::ChordReady) {
        for(int dy=-1;dy<=1;++dy) for(int dx=-1;dx<=1;++dx) {
            const int x=i%b.width()+dx, y=i/b.width()+dy;
            if(!(dx||dy) || x<0 || y<0 || x>=b.width() || y>=b.height()) continue;
            const int neighbor=y*b.width()+x;
            if(b.data()[neighbor]!=CellState::Unknown || out.marks[neighbor]==Mark::FlagForChord) continue;
            require(solutions>0 && hits[neighbor]==0, "chord would expose a possible mine");
        }
    }
}

static game::Board random_board(std::mt19937& rng, int w, int h) {
    std::vector<int> mines(w * h), cells(w * h);
    for (int& m : mines) m = rng() % 5 == 0;
    for (int i = 0; i < w * h; ++i) {
        if (mines[i]) { cells[i] = rng() % 4 == 0 ? 2 : 0; continue; }
        if (rng() % 2) continue;
        cells[i] = 10;
        for (int dy = -1; dy <= 1; ++dy) for (int dx = -1; dx <= 1; ++dx) {
            const int x = i % w + dx, y = i / w + dy;
            if ((dx || dy) && x >= 0 && y >= 0 && x < w && y < h) cells[i] += mines[y * w + x];
        }
    }
    return make_board(w, h, cells);
}

static void solver_tests() {
    verify_oracle(make_board(3, 1, {0, 11, 0}), 1, true);
    verify_oracle(make_board(3, 1, {11, 0, 0}), 1, true);
    verify_oracle(make_board(3, 1, {0, 0, 0}), 0, true);
    verify_oracle(make_board(3, 1, {0, 0, 0}), 3, false);
    verify_oracle(make_board(1, 1, {11}), -1, true);
    auto large = make_board(128, 32, std::vector<int>(4096, 0));
    large.set(0, 0, CellState::Number1);
    auto global = solve::compute_overlay(large, 100, false, 1);
    require(std::isfinite(global.mineProbability.back()), "large board probability overflow");
    require(std::abs(global.mineProbability.back() - 99.0 / 4092.0) < 1e-9, "incorrect outside probability");
    require(std::abs(global.mineProbability[1] - 1.0 / 3.0) < 1e-9, "incorrect frontier probability");
    const auto cancelled = solve::compute_overlay(large, 100, true, 4, [] { return true; });
    require(std::all_of(cancelled.marks.begin(), cancelled.marks.end(), [](Mark m) { return m == Mark::None; }), "cancelled work exposed hints");
    std::mt19937 rng(0x71a17);
    for (int i = 0; i < 800; ++i) {
        auto b = random_board(rng, 4, 3);
        // Choose a consistent mine count by constructing the same partial board's
        // valid totals through the oracle; unknown total still checks every guarantee.
        try {
            verify_oracle(b, -1, i % 2 == 0, i % 3 == 0 ? 4 : 1);
            for (int total = 0; total <= 6; ++total) verify_oracle(b, total, i % 2 == 0);
        } catch (const std::exception&) {
            std::cerr << "case " << i << ":";
            for (auto s : b.data()) std::cerr << ' ' << int(s);
            std::cerr << '\n';
            throw;
        }
    }
    // Rapid consecutive and simultaneous calls used to miss worker wakeups.
    std::vector<std::future<void>> tasks;
    for (int t = 0; t < 4; ++t) tasks.push_back(std::async(std::launch::async, [] {
        auto b = make_board(20, 1, {0,11,0,1,1,0,11,0,1,1,0,11,0,1,1,0,11,0,1,1});
        for (int i = 0; i < 100; ++i) verify_oracle(b, 4, true, i % 3 == 0 ? 2 : 8);
    }));
    for (auto& task : tasks) task.get();
    // Wide connected frontiers exercise sampling, bounded SAT, and belief
    // propagation. Check every guarantee against the planted ground truth.
    for(int width : {24, 65}) for(int round=0;round<12;++round) {
        const int n=width*3;
        std::vector<int> truth(n,0), cells(n,0);
        for(int i=0;i<n;++i) if(i/width!=1 || i%width%2) truth[i]=rng()%4==0;
        for(int x=0;x<width;x+=2) {
            int count=0;
            for(int dy=-1;dy<=1;++dy) for(int dx=-1;dx<=1;++dx)
                if((dx || dy) && x+dx>=0 && x+dx<width) count+=truth[(1+dy)*width+x+dx];
            cells[width+x]=10+count;
        }
        const auto out=solve::compute_overlay(make_board(width,3,cells),-1,true,4);
        for(int i=0;i<n;++i) {
            if(out.marks[i]==Mark::Safe) require(!truth[i], "approximation declared a mine safe");
            if(out.marks[i]==Mark::Mine || out.marks[i]==Mark::FlagForChord)
                require(truth[i], "approximation declared a safe cell mined");
            require(std::isfinite(out.mineProbability[i]), "nonfinite approximate probability");
        }
    }
}

static void protocol_tests() {
    proto::ParsedMessage out;
    require(proto::parse_message(R"({"type":"full","w":2,"h":1,"cells":[0,11]})", out), "valid full rejected");
    const char* invalid[] = {
        R"({"type":"full","w":-1,"h":1,"cells":[]})",
        R"({"type":"full","w":2147483647,"h":2,"cells":[]})",
        R"({"type":"full","w":2,"h":1,"cells":[0]})",
        R"({"type":"full","w":1,"h":1,"cells":[256]})",
        R"({"type":"full","w":1,"h":1,"cells":[3]})",
        R"({"type":"delta","updates":[]} trailing)",
        R"({"type":"delta","updates":[],"dpr":1e})",
        R"({"type":"delta","updates":[],"dpr":0})",
        R"({"type":"delta","updates":[],"rect_w":1e300})",
        R"({"type":"delta","updates":[],"type":"full"})",
        R"({"type":"delta","updates":[{"x":-1,"y":0,"s":0}]})",
        R"({"type":"delta","updates":[],"extra":"\q"})"
    };
    for (const char* json : invalid) require(!proto::parse_message(json, out), std::string("accepted malformed message: ") + json);
    require(proto::parse_message("{\"type\" : \"delta\",\"updates\":[],\"extra\":{\"type\":\"full\"}}", out) && out.type == proto::MsgType::Delta, "nested type overrides root");
    std::string deep = "{\"type\":\"delta\",\"updates\":[],\"extra\":" + std::string(100, '[') + "0" + std::string(100, ']') + "}";
    require(!proto::parse_message(deep, out), "unbounded JSON depth");
    require(out.type == proto::MsgType::Unknown, "failed parse retains previous type");
    require(proto::parse_message(R"({"type":"delta","updates":[],"extra":"\uD83D\uDE00"})", out), "valid Unicode escape rejected");
    require(!proto::parse_message(R"({"type":"delta","updates":[],"extra":"\uD83D"})", out), "unpaired surrogate accepted");
    std::mt19937 rng(71);
    const std::string valid = R"({"type":"full","w":2,"h":1,"cells":[0,11]})";
    for (size_t n = 0; n < valid.size(); ++n) require(!proto::parse_message(valid.substr(0, n), out), "truncated message accepted");
    for (int i = 0; i < 10000; ++i) {
        std::string input(rng() % 256, ' ');
        for (char& c : input) c = static_cast<char>(rng() & 255);
        proto::parse_message(input, out); // Must remain bounded for arbitrary bytes.
    }
    auto board = make_board(2, 1, {0, 11});
    try { board.resize(-1, 3); require(false, "negative resize accepted"); } catch (const std::invalid_argument&) {}
    require(board.width() == 2 && board.at(1, 0) == CellState::Number1, "invalid resize corrupts board");
    try { board.apply_full(std::vector<CellState>{}, 2, 1); require(false, "short board accepted"); } catch (const std::invalid_argument&) {}
    require(board.width() == 2 && board.at(1, 0) == CellState::Number1, "invalid full corrupts board");
}

int main(int argc, char** argv) {
    try {
        if (argc > 1 && std::string(argv[1]) == "--benchmark") {
            std::mt19937 rng(17);
            std::vector<game::Board> boards;
            for (int i = 0; i < 100; ++i) boards.push_back(random_board(rng, 30, 16));
            const auto start = std::chrono::steady_clock::now();
            size_t marks = 0;
            for (int rep = 0; rep < 3; ++rep) for (const auto& b : boards) marks += solve::compute_overlay(b, -1, true, 1).marks.size();
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
            std::cout << "300 expert-size snapshots: " << ms << " ms; " << marks << " output cells\n";
            return 0;
        }
        if (argc > 1 && std::string(argv[1]) == "--protocol") protocol_tests();
        else { solver_tests(); protocol_tests(); }
        std::cout << "Core regression checks passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
