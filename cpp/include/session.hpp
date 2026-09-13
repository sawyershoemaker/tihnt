#pragma once

#include <atomic>
#include <cstdint>
#include <optional>
#include <vector>

#include "board.hpp"
#include "proto.hpp"
#include "solver.hpp"

namespace app {

struct SolverJob {
    game::Board board;
    int minesTotal=-1;
    bool enableChords=true;
    uint64_t revision=0;
};

struct SolverResult {
    std::vector<solve::Mark> marks;
    uint64_t revision=0;
};

struct Presentation {
    std::vector<solve::Mark> marks;
    proto::GeometryMsg geometry;
    int width=0, height=0;
    int minesTotal=-1;
};

// The caller serializes state access; current_revision() is safe for solver workers.
class Session {
public:
    void reset();
    bool apply(proto::ParsedMessage&& message);
    bool set_chords(bool enabled);
    bool chords_enabled() const { return enableChords_; }
    uint64_t current_revision() const { return revision_.load(); }

    std::optional<SolverJob> take_job();
    bool accept(SolverResult&& result);
    bool take_presentation_changed();
    const Presentation& presentation() const { return presentation_; }
    std::optional<uint32_t> take_binding();

private:
    void invalidate();
    bool set_geometry(const proto::GeometryMsg& geometry);

    game::Board board_;
    Presentation presentation_;
    std::atomic<uint64_t> revision_{0};
    bool enableChords_=true;
    bool jobPending_=false, presentationChanged_=false, bindingChanged_=false;
    uint32_t targetPid_=0;
};

}
