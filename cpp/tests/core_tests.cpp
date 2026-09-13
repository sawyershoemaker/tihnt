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

static game::Board symmetric_frontier(int width) {
    std::vector<int> cells(width * 3, 0);
    cells[width] = cells[width * 2 - 1] = 1;
    for(int x=1;x<width-1;++x) cells[width+x]=11;
    return make_board(width, 3, cells);
}

static game::Board separated_frontiers(int count) {
    const int width=512, perRow=100;
    const int height=3*((count+perRow-1)/perRow);
    std::vector<int> cells(width*height, 1);
    for(int i=0;i<count;++i) {
        const int base=(i/perRow)*3*width+(i%perRow)*5;
        cells[base]=cells[base+2]=0;
        cells[base+1]=11;
    }
    cells.back()=0;
    return make_board(width,height,cells);
}

static void verify_column_frontier(const std::vector<int>& clue, int total) {
    const int width=static_cast<int>(clue.size());
    std::vector<int> cells(width*3,0);
    cells[width]=cells[width*2-1]=1;
    for(int x=1;x<width-1;++x) cells[width+x]=10+clue[x];
    std::vector<long double> mined(width,0);
    long double ways=0;
    // A column has two cells. Once its first two mine counts are chosen,
    // each remaining count follows directly from the preceding clue.
    for(int a=0;a<=2;++a) for(int b=0;b<=2;++b) {
        std::vector<int> counts(width,0);
        counts[0]=a; counts[1]=b;
        bool valid=true;
        int sum=a+b;
        for(int x=2;x<width;++x) {
            counts[x]=clue[x-1]-counts[x-2]-counts[x-1];
            if(counts[x]<0 || counts[x]>2) { valid=false; break; }
            sum+=counts[x];
        }
        if(!valid || (total>=0 && sum!=total)) continue;
        long double weight=1;
        for(int count : counts) if(count==1) weight*=2;
        ways+=weight;
        for(int x=0;x<width;++x) mined[x]+=weight*counts[x]/2;
    }
    const auto out=solve::compute_overlay(make_board(width,3,cells),total,false,1);
    for(int y : {0,2}) for(int x=0;x<width;++x) {
        const int i=y*width+x;
        if(ways==0) {
            require(out.marks[i]==Mark::None, "inconsistent column frontier exposed hints");
            continue;
        }
        const double expected=static_cast<double>(mined[x]/ways);
        require(std::abs(out.mineProbability[i]-expected)<1e-9,
            "column frontier w="+std::to_string(width)+" total="+std::to_string(total)+" x="+std::to_string(x)+
            " probability="+std::to_string(out.mineProbability[i])+" expected="+std::to_string(expected));
        if(expected==0) require(out.marks[i]==Mark::Safe, "column oracle safe cell not proved");
        if(expected==1) require(out.marks[i]==Mark::Mine, "column oracle mine not proved");
    }
}

static void large_frontier_tests() {
    // Every three adjacent columns contain one mine. The two cells in each
    // column are interchangeable, so there are three weighted periodic layouts.
    auto board=symmetric_frontier(59);
    const auto local=solve::compute_overlay(board,-1,false,1);
    const auto conditioned=solve::compute_overlay(board,19,false,1);
    for(int y : {0,2}) for(int x=0;x<59;++x) {
        const int i=y*59+x;
        const double expected=x%3==2 ? 0.1 : 0.2;
        require(std::abs(local.mineProbability[i]-expected)<1e-9, "symmetric frontier lost exact weights");
        require(std::abs(conditioned.mineProbability[i]-(x%3==2 ? 0.5 : 0.0))<1e-9,
            "large frontier ignores the total mine count");
        if(x%3!=2) require(conditioned.marks[i]==Mark::Safe, "global count did not prove the safe columns");
    }
    auto separated=separated_frontiers(300);
    const auto saturated=solve::compute_overlay(separated,301,false,4);
    require(saturated.marks.back()==Mark::Mine, "large disconnected frontier lost outside certainty");
    const auto exhausted=solve::compute_overlay(separated,300,false,4);
    require(exhausted.marks.back()==Mark::Safe, "large disconnected frontier lost outside safe move");
    const auto impossible=solve::compute_overlay(separated,299,false,4);
    require(std::all_of(impossible.marks.begin(),impossible.marks.end(),[](Mark m){return m==Mark::None;}),
        "globally inconsistent frontier exposed hints");

    const int count=200,width=512,height=9;
    std::vector<int> cells(width*height,1);
    for(int c=0;c<count;++c) {
        const int base=(c/70)*3*width+(c%70)*7;
        cells[base]=cells[base+2]=cells[base+4]=0;
        cells[base+1]=cells[base+3]=11;
    }
    cells.back()=0;
    const auto varied=solve::compute_overlay(make_board(width,height,cells),225,false,4);
    const double p=25.0/201.0;
    for(int c=0;c<count;++c) {
        const int base=(c/70)*3*width+(c%70)*7;
        require(std::abs(varied.mineProbability[base]-p)<1e-9 &&
            std::abs(varied.mineProbability[base+2]-(1-p))<1e-9 &&
            std::abs(varied.mineProbability[base+4]-p)<1e-9,
            "backward conditioning lost component count weights");
    }
    require(std::abs(varied.mineProbability.back()-p)<1e-9, "backward conditioning lost outside count weights");

    int checks=0;
    const auto stopped=solve::compute_overlay(separated,300,false,1,[&]{return ++checks>650;});
    require(checks>650 && std::all_of(stopped.marks.begin(),stopped.marks.end(),[](Mark m){return m==Mark::None;}),
        "cancellation during global conditioning exposed partial hints");

    std::mt19937 rng(0x19a51);
    for(int w : {17,33,59,64}) for(int round=0;round<16;++round) {
        std::vector<int> truth(w),clue(w);
        int total=0;
        for(int& value : truth) { value=round<2 ? round+1 : rng()%3; total+=value; }
        for(int x=1;x<w-1;++x) clue[x]=truth[x-1]+truth[x]+truth[x+1];
        for(int mines : {-1,total,total-1}) verify_column_frontier(clue,mines);
    }
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
    const std::string geometry=R"(,"rect_l":-20.25,"rect_t":10.5,"rect_w":80,"rect_h":40)";
    const std::string clipFields[]={R"("clip_l":-12.25)",R"("clip_t":15.5)",R"("clip_w":60.5)",R"("clip_h":20.25)"};
    for(unsigned mask=0;mask<16;++mask) {
        std::string message=R"({"type":"delta","updates":[])"+geometry;
        for(int i=0;i<4;++i) if(mask&(1u<<i)) message+=","+clipFields[i];
        const bool valid=mask==0 || mask==15;
        require(proto::parse_message(message+"}",out)==valid,"partial clipping rectangle was accepted");
        if(valid) require(out.delta.has_geometry && out.delta.has_clip==(mask==15),"optional clipping presence was lost");
        if(mask==15) require(out.delta.clip_l==-12.25 && out.delta.clip_t==15.5 &&
            out.delta.clip_w==60.5 && out.delta.clip_h==20.25,"fractional clipping values were changed");
    }
    const std::string clip=R"(,"clip_l":-1000000,"clip_t":1000000,"clip_w":16384,"clip_h":0)";
    require(proto::parse_message(R"({"type":"full","w":2,"h":1,"cells":[0,11])"+geometry+clip+"}",out) &&
        out.full.has_clip && out.full.clip_h==0,"bounded or empty full clipping was rejected");
    require(!proto::parse_message(R"({"type":"delta","updates":[])"+clip+"}",out),"clip without board geometry was accepted");
    const std::string fields[]={"rect_l","rect_t","rect_w","rect_h"};
    for(const auto& missing:fields) {
        std::string message=R"({"type":"delta","updates":[])"+clip;
        for(const auto& field:fields) if(field!=missing) message+=",\""+field+"\":10";
        require(!proto::parse_message(message+"}",out),"clip accepted incomplete board geometry");
    }
    const std::string names[]={"clip_l","clip_t","clip_w","clip_h"};
    for(int axis=0;axis<4;++axis) for(const char* invalidValue : {"null","true","\"1\"","1e309","1000000.25","-1000000.25","NaN"}) {
        std::string message=R"({"type":"delta","updates":[])"+geometry;
        for(int i=0;i<4;++i) message+=",\""+names[i]+"\":"+(i==axis ? invalidValue : "1");
        require(!proto::parse_message(message+"}",out),"invalid clipping coordinate was accepted");
    }
    for(int axis=2;axis<4;++axis) for(const char* invalidValue : {"-0.1","16384.1"}) {
        std::string message=R"({"type":"delta","updates":[])"+geometry;
        for(int i=0;i<4;++i) message+=",\""+names[i]+"\":"+(i==axis ? invalidValue : "1");
        require(!proto::parse_message(message+"}",out),"invalid clipping dimension was accepted");
    }
    require(proto::parse_message(R"({"type":"delta","updates":[]})",out) && !out.delta.has_clip &&
        !out.delta.has_geometry,"legacy delta retained previous clipping metadata");
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
        if(argc>1 && std::string(argv[1])=="--benchmark-frontier") {
            const auto board=symmetric_frontier(59);
            const auto start=std::chrono::steady_clock::now();
            double sum=0;
            for(int i=0;i<100;++i) sum+=solve::compute_overlay(board,19,false,1).mineProbability[0];
            const double ms=std::chrono::duration<double,std::milli>(std::chrono::steady_clock::now()-start).count();
            std::cout << "100 connected 118-cell frontiers: " << ms << " ms; probability sum " << sum << '\n';
            return 0;
        }
        if(argc>1 && std::string(argv[1])=="--large-frontier") {
            large_frontier_tests();
            std::cout << "Large frontier checks passed\n";
            return 0;
        }
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
        else { solver_tests(); large_frontier_tests(); protocol_tests(); }
        std::cout << "Core regression checks passed\n";
        return 0;
    } catch (const std::exception& e) {
        std::cerr << e.what() << '\n';
        return 1;
    }
}
