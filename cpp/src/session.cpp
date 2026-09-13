#include "session.hpp"

#include <utility>

namespace app {

void Session::invalidate() {
    ++revision_;
    jobPending_=!board_.data().empty();
    presentation_.width=board_.width();
    presentation_.height=board_.height();
    presentation_.marks.assign(board_.data().size(),solve::Mark::None);
    presentationChanged_=true;
}

void Session::reset() {
    board_.resize(0,0);
    presentation_={};
    targetPid_=0;
    bindingChanged_=true;
    invalidate();
}

bool Session::set_geometry(const proto::GeometryMsg& geometry) {
    const auto& old=presentation_.geometry;
    const bool changed=old.rect_l!=geometry.rect_l || old.rect_t!=geometry.rect_t ||
        old.rect_w!=geometry.rect_w || old.rect_h!=geometry.rect_h ||
        old.vv_x!=geometry.vv_x || old.vv_y!=geometry.vv_y ||
        old.vv_scale!=geometry.vv_scale || old.dpr!=geometry.dpr || old.has_clip!=geometry.has_clip ||
        (geometry.has_clip && (old.clip_l!=geometry.clip_l || old.clip_t!=geometry.clip_t ||
            old.clip_w!=geometry.clip_w || old.clip_h!=geometry.clip_h));
    presentation_.geometry=geometry;
    presentationChanged_|=changed;
    return changed;
}

bool Session::apply(proto::ParsedMessage&& message) {
    bool changed=false, moved=false;
    if(message.type==proto::MsgType::Full) {
        auto& full=message.full;
        const bool boardChanged=board_.width()!=full.w || board_.height()!=full.h || board_.data()!=full.cells;
        changed=boardChanged || presentation_.minesTotal!=full.mines_total;
        if(boardChanged) board_.apply_full(std::move(full.cells),full.w,full.h);
        presentation_.minesTotal=full.mines_total;
        moved=set_geometry(full);
    } else if(message.type==proto::MsgType::Delta) {
        if(board_.width()==0) return false;
        const auto& delta=message.delta;
        bool ordered=true;
        int previous=-1;
        for(const auto& update:delta.updates) {
            if(update.x<0 || update.y<0 || update.x>=board_.width() || update.y>=board_.height() ||
                !game::valid_cell_state(static_cast<int>(update.state))) return false;
            const int index=update.y*board_.width()+update.x;
            ordered&=index>previous;
            previous=index;
        }
        if(ordered) {
            for(const auto& update:delta.updates) {
                if(board_.at(update.x,update.y)!=update.state) {
                    board_.set(update.x,update.y,update.state);
                    changed=true;
                }
            }
        } else {
            // Repeated coordinates use the final update; transient values must not cancel valid work.
            auto cells=board_.data();
            for(const auto& update:delta.updates) cells[update.y*board_.width()+update.x]=update.state;
            changed=cells!=board_.data();
            if(changed) board_.apply_full(std::move(cells),board_.width(),board_.height());
        }
        if(delta.has_mines_total && presentation_.minesTotal!=delta.mines_total) {
            presentation_.minesTotal=delta.mines_total;
            changed=true;
        }
        if(delta.has_geometry) moved=set_geometry(delta);
    } else if(message.type==proto::MsgType::Bind) {
        const auto pid=static_cast<uint32_t>(message.bind.pid);
        if(targetPid_==pid) return false;
        targetPid_=pid;
        bindingChanged_=true;
        return true;
    }
    if(changed) invalidate();
    return changed || moved;
}

bool Session::set_chords(bool enabled) {
    if(enableChords_==enabled) return false;
    enableChords_=enabled;
    invalidate();
    return true;
}

std::optional<SolverJob> Session::take_job() {
    if(!jobPending_) return std::nullopt;
    SolverJob job{board_,presentation_.minesTotal,enableChords_,current_revision()};
    jobPending_=false;
    return job;
}

bool Session::accept(SolverResult&& result) {
    if(result.revision!=current_revision() || result.marks.size()!=board_.data().size()) return false;
    if(presentation_.marks!=result.marks) {
        presentation_.marks=std::move(result.marks);
        presentationChanged_=true;
    }
    return true;
}

bool Session::take_presentation_changed() {
    return std::exchange(presentationChanged_,false);
}

std::optional<uint32_t> Session::take_binding() {
    if(!std::exchange(bindingChanged_,false)) return std::nullopt;
    return targetPid_;
}

}
