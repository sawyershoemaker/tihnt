#include "session.hpp"

#include <algorithm>
#include <iostream>
#include <stdexcept>
#include <string>
#include <utility>

using solve::Mark;

static void require(bool ok,const char* message) {
    if(!ok) throw std::runtime_error(message);
}

static bool deliver(app::Session& session,const std::string& json) {
    proto::ParsedMessage message;
    require(proto::parse_message(json,message),"invalid test message");
    return session.apply(std::move(message));
}

static const char* initial=R"({"type":"full","w":3,"h":1,"cells":[10,0,0],"mines_total":0,"rect_l":40,"rect_t":50,"rect_w":60,"rect_h":20,"dpr":1})";

static app::SolverJob take_job(app::Session& session) {
    auto job=session.take_job();
    require(job.has_value(),"expected a solver job");
    return std::move(*job);
}

static app::SolverResult solve_job(const app::SolverJob& job) {
    auto result=solve::compute_overlay(job.board,job.minesTotal,job.enableChords,1);
    return {std::move(result.marks),job.revision};
}

static void baseline_and_duplicate_messages() {
    app::Session session;
    require(!deliver(session,R"({"type":"delta","updates":[{"x":0,"y":0,"s":10}]})"),"a delta established a board");
    require(!session.take_job() && !session.take_presentation_changed(),"ignored delta produced work");
    require(deliver(session,initial),"initial snapshot was ignored");
    auto job=take_job(session);
    require(session.take_presentation_changed(),"new board did not clear previous marks");
    require(session.accept(solve_job(job)),"current solver result was rejected");
    require(session.presentation().marks[1]==Mark::Safe && session.presentation().marks[2]==Mark::Safe,"safe marks were not presented");
    require(session.take_presentation_changed(),"solved marks were not published");
    const auto revision=session.current_revision();
    require(!deliver(session,initial),"identical snapshot scheduled redundant work");
    require(!deliver(session,R"({"type":"delta","updates":[{"x":1,"y":0,"s":0}],"mines_total":0})"),"no-op delta scheduled redundant work");
    require(!deliver(session,R"({"type":"delta","updates":[{"x":1,"y":0,"s":2},{"x":1,"y":0,"s":0}]})"),"net-no-op batch cancelled valid work");
    require(!session.take_job() && !session.take_presentation_changed(),"unchanged board recomputed or repainted");
    require(session.current_revision()==revision,"unchanged board cancelled a solve");
}

static void latest_geometry_preserves_solves() {
    app::Session session;
    deliver(session,initial);
    auto job=take_job(session);
    session.take_presentation_changed();
    require(deliver(session,R"({"type":"delta","updates":[],"rect_l":15,"rect_t":-12,"rect_w":90,"rect_h":30,"vv_x":2,"vv_y":3,"vv_scale":1.5,"dpr":2})"),"viewport movement was ignored");
    require(!session.take_job(),"scrolling recomputed the board");
    require(session.current_revision()==job.revision,"scrolling cancelled a valid solve");
    require(session.accept(solve_job(job)),"geometry change rejected a valid result");
    const auto& view=session.presentation();
    require(view.geometry.rect_l==15 && view.geometry.rect_t==-12 && view.geometry.rect_w==90,"a solver result rewound current geometry");
    require(view.geometry.vv_x==2 && view.geometry.vv_y==3 && view.geometry.vv_scale==1.5 && view.geometry.dpr==2,"viewport metadata was lost");
    require(view.minesTotal==0,"geometry-only delta reset the mine total");
    require(session.take_presentation_changed(),"latest geometry was not published");
}

static void changed_boards_retire_results() {
    app::Session session;
    deliver(session,initial);
    const auto oldJob=take_job(session);
    auto oldResult=solve_job(oldJob);
    require(session.accept(solve_job(oldJob)),"initial result rejected");
    session.take_presentation_changed();
    deliver(session,R"({"type":"full","w":3,"h":1,"cells":[0,11,0],"mines_total":1,"rect_l":40,"rect_t":50,"rect_w":60,"rect_h":20})");
    require(std::all_of(session.presentation().marks.begin(),session.presentation().marks.end(),[](Mark mark){return mark==Mark::None;}),"same-size new game retained old hints");
    require(session.take_presentation_changed(),"new game did not publish cleared marks");
    require(!session.accept(std::move(oldResult)),"obsolete game result was accepted");
    auto newJob=take_job(session);
    require(newJob.board.at(1,0)==game::CellState::Number1 && newJob.minesTotal==1,"new game job used old state");
    require(oldJob.board.at(0,0)==game::CellState::Number0 && oldJob.minesTotal==0,"pending job snapshot mutated");
    require(!session.accept({{Mark::Safe},newJob.revision}),"wrong-sized result was accepted");
    require(session.accept(solve_job(newJob)),"current new-game result rejected");
}

static void clipping_preserves_solves() {
    app::Session session;
    deliver(session,initial);
    session.take_presentation_changed();
    const auto revision=session.current_revision();
    const char* clipped=R"({"type":"delta","updates":[],"rect_l":40,"rect_t":50,"rect_w":60,"rect_h":20,"clip_l":50.25,"clip_t":55.5,"clip_w":30.5,"clip_h":10})";
    require(deliver(session,clipped) && session.take_presentation_changed(),"clip-only change was not presented");
    auto job=take_job(session);
    require(job.revision==revision && !session.take_job(),"clipping replaced or duplicated a pending solve");
    require(!deliver(session,clipped) && !session.take_presentation_changed(),"identical clipping republished the presentation");
    require(deliver(session,R"({"type":"delta","updates":[],"rect_l":40,"rect_t":50,"rect_w":60,"rect_h":20,"clip_l":50.25,"clip_t":55.5,"clip_w":0,"clip_h":10})"),"empty clipping was ignored");
    require(session.current_revision()==revision && !session.take_job(),"empty clipping invalidated a solve");
    require(session.accept(solve_job(job)),"clipping rejected a valid in-flight result");
    const auto marks=session.presentation().marks;
    require(session.presentation().geometry.has_clip && session.presentation().geometry.clip_w==0 && marks[1]==Mark::Safe,
        "a completed solve lost the latest clip or safe marks");
    session.take_presentation_changed();
    require(!deliver(session,R"({"type":"delta","updates":[]})") && !session.take_presentation_changed() &&
        session.presentation().geometry.has_clip,"geometry-less delta reset clipping");
    require(deliver(session,clipped) && session.take_presentation_changed(),"restored clipping did not update presentation");
    require(session.presentation().marks==marks && !session.take_job(),"clip restoration discarded completed marks");
    require(!deliver(session,R"({"type":"delta","updates":[{"x":3,"y":0,"s":2}],"rect_l":40,"rect_t":50,"rect_w":60,"rect_h":20})") &&
        session.presentation().geometry.has_clip,"invalid delta partially reset clipping");
    require(deliver(session,R"({"type":"delta","updates":[],"rect_l":40,"rect_t":50,"rect_w":60,"rect_h":20})") &&
        !session.presentation().geometry.has_clip,"fresh legacy geometry did not remove clipping");
    require(deliver(session,clipped),"reapplying clipping was ignored");
    require(deliver(session,initial) && !session.presentation().geometry.has_clip,
        "fresh legacy snapshot did not remove clipping");
    require(session.current_revision()==revision && session.presentation().marks==marks && !session.take_job(),
        "legacy clip removal invalidated completed work");
}

static void deltas_are_atomic_and_keep_optional_state() {
    app::Session session;
    deliver(session,initial);
    const auto original=take_job(session);
    session.take_presentation_changed();
    require(!deliver(session,R"({"type":"delta","updates":[{"x":1,"y":0,"s":2},{"x":3,"y":0,"s":2}],"mines_total":2,"rect_l":99,"rect_t":99,"rect_w":60,"rect_h":20})"),"out-of-board update was accepted");
    require(!session.take_job() && !session.take_presentation_changed(),"invalid delta partially changed state");
    require(session.current_revision()==original.revision && session.presentation().minesTotal==0 && session.presentation().geometry.rect_l==40,"invalid delta changed metadata");
    deliver(session,R"({"type":"delta","updates":[{"x":2,"y":0,"s":10}]})");
    const auto changed=take_job(session);
    require(changed.board.at(1,0)==game::CellState::Unknown && changed.board.at(2,0)==game::CellState::Number0,"delta failed to update exactly its cells");
    require(changed.minesTotal==0 && session.presentation().geometry.rect_l==40,"omitted metadata was reset");
    deliver(session,R"({"type":"delta","updates":[{"x":2,"y":0,"s":0},{"x":1,"y":0,"s":2},{"x":2,"y":0,"s":11}]})");
    const auto unordered=take_job(session);
    require(unordered.board.at(1,0)==game::CellState::Mine && unordered.board.at(2,0)==game::CellState::Number1,"unordered updates did not retain the final state per cell");
}

static void settings_and_latest_work() {
    app::Session session;
    deliver(session,initial);
    auto oldJob=take_job(session);
    require(!session.set_chords(true),"unchanged chord setting created work");
    require(session.set_chords(false),"chord setting did not invalidate work");
    deliver(session,R"({"type":"delta","updates":[{"x":1,"y":0,"s":10}]})");
    deliver(session,R"({"type":"delta","updates":[],"mines_total":1})");
    auto latest=take_job(session);
    require(!latest.enableChords && latest.minesTotal==1 && latest.board.at(1,0)==game::CellState::Number0,"coalesced job did not contain the latest board and settings");
    require(!session.take_job(),"coalesced job remained pending twice");
    require(!session.accept(solve_job(oldJob)),"settings change accepted a previous result");
}

static void connection_resets_and_empty_boards() {
    app::Session session;
    deliver(session,initial);
    auto oldJob=take_job(session);
    session.set_chords(false);
    deliver(session,R"({"type":"bind","pid":123})");
    require(session.take_binding()==123u,"target binding was not delivered");
    require(!session.take_binding(),"target binding was delivered twice");
    session.reset();
    require(session.take_binding()==0u,"connection reset kept an old process binding");
    const auto& empty=session.presentation();
    require(empty.width==0 && empty.height==0 && empty.marks.empty() && empty.minesTotal==-1 && empty.geometry.rect_w==0,"connection reset retained previous board state");
    require(session.take_presentation_changed(),"connection reset did not clear the overlay");
    require(!session.take_job(),"empty board invoked the solver");
    require(!session.accept(solve_job(oldJob)),"connection reset accepted an old result");
    require(!deliver(session,R"({"type":"delta","updates":[]})"),"delta after reconnect established a baseline");
    deliver(session,initial);
    auto newJob=take_job(session);
    require(!newJob.enableChords,"reconnect reset the user's chord preference");
    require(newJob.revision!=oldJob.revision,"identical board after reconnect reused a previous revision");
    deliver(session,R"({"type":"full","w":0,"h":0,"cells":[]})");
    require(session.presentation().marks.empty() && !session.take_job(),"empty full snapshot did not retire the board");
    require(!session.accept(solve_job(newJob)),"empty full snapshot accepted a delayed result");
}

int main() {
    try {
        baseline_and_duplicate_messages();
        latest_geometry_preserves_solves();
        clipping_preserves_solves();
        changed_boards_retire_results();
        deltas_are_atomic_and_keep_optional_state();
        settings_and_latest_work();
        connection_resets_and_empty_boards();
        std::cout<<"Session checks passed: baseline, snapshots, deltas, resets, settings, immutable jobs, stale results, latest geometry.\n";
        return 0;
    } catch(const std::exception& error) {
        std::cerr<<error.what()<<'\n';
        return 1;
    }
}
