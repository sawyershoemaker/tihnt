#include "solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <unordered_map>
#include <unordered_set>
#include <thread>
#include <random>
#include <atomic>
#include <condition_variable>
#include <mutex>

namespace {

enum class Sat { No, Yes, Unknown };

class WorkerPool {
public:
    static WorkerPool& instance() { static WorkerPool pool; return pool; }
    void run(int requested, int count, const std::function<void(int)>& fn) {
        if (count <= 0) return;
        std::lock_guard<std::mutex> submission(runMutex_);
        const int active = std::min({std::max(1, requested), count, 8});
        std::unique_lock<std::mutex> lock(mtx_);
        while (static_cast<int>(workers_.size()) < active) {
            const int index = static_cast<int>(workers_.size());
            const uint64_t seen = generation_;
            workers_.emplace_back([this, index, seen] { worker_loop(index, seen); });
        }
        currentTask_ = fn;
        taskCount_ = count;
        activeThreads_ = remaining_ = active;
        error_ = nullptr;
        nextTask_.store(0, std::memory_order_relaxed);
        ++generation_;
        cv_.notify_all();
        doneCv_.wait(lock, [&] { return remaining_ == 0; });
        currentTask_ = {};
        if (error_) std::rethrow_exception(error_);
    }
private:
    ~WorkerPool() {
        {
            std::lock_guard<std::mutex> lock(mtx_);
            stopping_ = true;
        }
        cv_.notify_all();
        for (auto& worker : workers_) worker.join();
    }
    void worker_loop(int index, uint64_t seen) {
        for (;;) {
            std::unique_lock<std::mutex> lock(mtx_);
            cv_.wait(lock, [&] { return stopping_ || generation_ != seen; });
            if (stopping_) return;
            seen = generation_;
            if (index >= activeThreads_) continue;
            auto task = currentTask_;
            const int count = taskCount_;
            lock.unlock();
            try {
                for (;;) {
                    int id = nextTask_.fetch_add(1, std::memory_order_relaxed);
                    if (id >= count) break;
                    task(id);
                }
            } catch (...) {
                std::lock_guard<std::mutex> guard(mtx_);
                if (!error_) error_ = std::current_exception();
            }
            lock.lock();
            if (--remaining_ == 0) doneCv_.notify_one();
        }
    }
    std::mutex runMutex_, mtx_;
    std::condition_variable cv_, doneCv_;
    std::vector<std::thread> workers_;
    std::function<void(int)> currentTask_;
    std::exception_ptr error_;
    std::atomic<int> nextTask_{0};
    uint64_t generation_ = 0;
    int taskCount_ = 0, activeThreads_ = 0, remaining_ = 0;
    bool stopping_ = false;
};

}

namespace solve {

using game::Board;
using game::CellState;

static inline int to_index(int x, int y, int w){ return y*w + x; }

static inline bool in_bounds(int x, int y, int w, int h){ return x>=0 && y>=0 && x<w && y<h; }

static inline int cell_number(CellState s){
	int si = static_cast<int>(s);
	int base = static_cast<int>(CellState::Number0);
	if(si >= base && si <= base + 8){ return si - base; }
	return -1;
}

static thread_local std::mt19937 rng(std::random_device{}());

Overlay compute_overlay(const Board& board, int totalMines, bool enableChords, int threads, const std::function<bool()>& isCancelled){
	constexpr int UNKNOWN_TOTAL = -1;
	Overlay ov{};
	if(totalMines < UNKNOWN_TOTAL) totalMines = UNKNOWN_TOTAL;
	const int w = board.width();
	const int h = board.height();
	const int N = w * h;
	ov.marks.assign(N, Mark::None);
    auto cancelled = [&] { return isCancelled && isCancelled(); };
    std::atomic<int64_t> searchWork{8000000};
    auto searchBlock = [&] {
        // Account in blocks to keep synchronization out of the hot recursion.
        return !cancelled() && searchWork.fetch_sub(1024, std::memory_order_relaxed) > 0;
    };
    auto invalid = [&] { Overlay empty; empty.marks.assign(N, Mark::None); empty.mineProbability.assign(N, -1.0); return empty; };
    if(cancelled()) return invalid();

	if(w==0 || h==0) return ov;

	struct NeighborCache { int w=-1, h=-1; std::vector<std::array<int,8>> nbr; std::vector<int> cnt; };
	thread_local NeighborCache cache;
    if(cache.w!=w || cache.h!=h){
        cache.w=w; cache.h=h; cache.nbr.assign(N, {}); cache.cnt.assign(N, 0);
        for(int y=0;y<h;++y){
            for(int x=0;x<w;++x){
                int idx = to_index(x,y,w);
                int c=0;
                for(int dy=-1; dy<=1; ++dy){
                    for(int dx=-1; dx<=1; ++dx){
                        if(dx==0 && dy==0) continue;
                        int nx=x+dx, ny=y+dy;
                        if(!in_bounds(nx,ny,w,h)) continue;
                        cache.nbr[idx][c++] = to_index(nx,ny,w);
                    }
                }
                cache.cnt[idx]=c;
            }
        }
    }
    const auto& neighbors = cache.nbr;
    const auto& neighborCounts = cache.cnt;
	const auto& cells = board.data();

	// precompute numeric values for all cells to avoid repeated decoding
	std::vector<int> numbers(N, -1);
	for(int i=0;i<N;++i){ numbers[i] = cell_number(cells[i]); }

	std::vector<int> numberCells;
    numberCells.reserve(N);
	for(int i=0;i<N;++i){ if(numbers[i] >= 0) numberCells.push_back(i); }

    std::vector<uint8_t> inQueue(N,0);
    std::vector<int> queue;
    queue.reserve(numberCells.size());
    for(int i : numberCells){ queue.push_back(i); inQueue[i]=1; }

	auto enqueueNbrNumbers = [&](int cellIdx){
        for(int k=0;k<neighborCounts[cellIdx];++k){
            int nb = neighbors[cellIdx][k];
			if(numbers[nb] >= 0 && !inQueue[nb]){ queue.push_back(nb); inQueue[nb]=1; }
        }
    };

    while(!queue.empty()){
        int idxCenter = queue.back(); queue.pop_back(); inQueue[idxCenter]=0;
		int num = numbers[idxCenter];
        if(num < 0) continue;
        int knownMines = 0;
        int unknownIdx[8]; int ucount = 0;
		for(int k=0;k<neighborCounts[idxCenter];++k){
            int nb = neighbors[idxCenter][k];
			CellState s = cells[nb];
            if(s == CellState::Mine || ov.marks[nb] == Mark::Mine){ knownMines++; continue; }
            if(s == CellState::Unknown && ov.marks[nb] != Mark::Safe){ unknownIdx[ucount++] = nb; }
        }
        int remaining = num - knownMines;
        if(remaining < 0 || remaining > ucount){
            return invalid();
        }
        bool any=false;
        if(remaining == 0 && ucount > 0){
            for(int i=0;i<ucount;++i){ int u=unknownIdx[i]; if(ov.marks[u] != Mark::Safe){ ov.marks[u]=Mark::Safe; any=true; enqueueNbrNumbers(u);} }
        } else if(remaining == ucount && ucount>0){
            for(int i=0;i<ucount;++i){ int u=unknownIdx[i]; if(ov.marks[u] != Mark::Mine){ ov.marks[u]=Mark::Mine; any=true; enqueueNbrNumbers(u);} }
        }
		// Rebuild both sets after deductions; mixing old and new sets is unsound.
        if(any) continue;
		if(ucount>0){
            for(int k=0;k<neighborCounts[idxCenter];++k){
                int nbNumIdx = neighbors[idxCenter][k];
				int num2 = numbers[nbNumIdx]; if(num2<0) continue;
                int known2=0; int u2c=0; int u2[8];
                for(int t=0;t<neighborCounts[nbNumIdx];++t){
                    int q = neighbors[nbNumIdx][t];
					CellState sq = cells[q];
                    if(sq==CellState::Mine || ov.marks[q]==Mark::Mine){ known2++; continue; }
                    if(sq==CellState::Unknown && ov.marks[q]!=Mark::Safe){ u2[u2c++]=q; }
                }
                int rem1 = remaining;
                int rem2 = num2 - known2;
                if(rem2<0 || rem2>u2c) return invalid();
                auto isIn = [&](int v, const int* arr, int n){ for(int i=0;i<n;++i) if(arr[i]==v) return true; return false; };
                int common=0;
                for(int i=0;i<ucount;++i){ if(isIn(unknownIdx[i], u2, u2c)) common++; }
                if(common==ucount && ucount<u2c){
                    int diff = rem2 - rem1;
                    if(diff==0){
                        for(int i=0;i<u2c;++i){ int q=u2[i]; if(!isIn(q, unknownIdx, ucount)){ if(ov.marks[q]!=Mark::Safe){ ov.marks[q]=Mark::Safe; any=true; enqueueNbrNumbers(q);} }
                        }
                    } else if(diff == (u2c-common)){
                        for(int i=0;i<u2c;++i){ int q=u2[i]; if(!isIn(q, unknownIdx, ucount)){ if(ov.marks[q]!=Mark::Mine){ ov.marks[q]=Mark::Mine; any=true; enqueueNbrNumbers(q);} }
                        }
                    }
                } else if(common==u2c && u2c<ucount){
                    int diff = rem1 - rem2;
                    if(diff==0){
                        for(int i=0;i<ucount;++i){ int q=unknownIdx[i]; if(!isIn(q, u2, u2c)){ if(ov.marks[q]!=Mark::Safe){ ov.marks[q]=Mark::Safe; any=true; enqueueNbrNumbers(q);} }
                        }
                    } else if(diff == (ucount-common)){
                        for(int i=0;i<ucount;++i){ int q=unknownIdx[i]; if(!isIn(q, u2, u2c)){ if(ov.marks[q]!=Mark::Mine){ ov.marks[q]=Mark::Mine; any=true; enqueueNbrNumbers(q);} }
                        }
                    }
                } else {
                    // general difference rule for overlapping sets (neither subset of the other)
                    int only1 = ucount - common;
                    int only2 = u2c - common;
                    if(only1 > 0 && (rem1 - rem2) == only1){
                        for(int i=0;i<ucount;++i){
                            int q = unknownIdx[i];
                            if(!isIn(q, u2, u2c)){
                                if(ov.marks[q] != Mark::Mine){ ov.marks[q] = Mark::Mine; any=true; enqueueNbrNumbers(q); }
                            }
                        }
                    }
                    if(only2 > 0 && (rem2 - rem1) == only2){
                        for(int i=0;i<u2c;++i){
                            int q = u2[i];
                            if(!isIn(q, unknownIdx, ucount)){
                                if(ov.marks[q] != Mark::Mine){ ov.marks[q] = Mark::Mine; any=true; enqueueNbrNumbers(q); }
                            }
                        }
                    }
                }
                if(any) { if(!inQueue[idxCenter]) { queue.push_back(idxCenter); inQueue[idxCenter]=1; } break; }
            }
        }
        if(cancelled()) return invalid();
	}

    int fixedMines = 0, unresolved = 0;
    for(int i=0; i<N; ++i) {
        if(cells[i] == CellState::Mine || ov.marks[i] == Mark::Mine) ++fixedMines;
        else if(cells[i] == CellState::Unknown && ov.marks[i] != Mark::Safe) ++unresolved;
    }
    if(totalMines >= 0) {
        const int remaining = totalMines - fixedMines;
        if(remaining < 0 || remaining > unresolved) return invalid();
        if(remaining == 0 || remaining == unresolved) {
            for(int i=0; i<N; ++i) {
                if(cells[i] == CellState::Unknown && ov.marks[i] == Mark::None)
                    ov.marks[i] = remaining == 0 ? Mark::Safe : Mark::Mine;
            }
        }
    }
    for(int idx : numberCells) {
        int known=0, unknown=0;
        for(int k=0; k<neighborCounts[idx]; ++k) {
            int nb=neighbors[idx][k];
            if(cells[nb]==CellState::Mine || ov.marks[nb]==Mark::Mine) ++known;
            else if(cells[nb]==CellState::Unknown && ov.marks[nb]!=Mark::Safe) ++unknown;
        }
        if(numbers[idx] < known || numbers[idx] > known+unknown) return invalid();
    }

	// determine if we have any guaranteed safe so far
	for(const auto m : ov.marks){ if(m == Mark::Safe){ ov.hasGuaranteedSafe = true; break; } }

        std::vector<uint8_t> isFrontier(N,0);
        for(int i=0;i<N;++i){ if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Safe && ov.marks[i]!=Mark::Mine){ for(int k=0;k<neighborCounts[i];++k){ int nb=neighbors[i][k]; if(cell_number(board.data()[nb])>=0){ isFrontier[i]=1; break; } } } }

        std::vector<int> compId(N,-1);
        std::vector<int> uf(N, -1);
        auto findp = [&](int v){ int r=v; while(uf[r]!=r) r=uf[r]; while(uf[v]!=v){ int t=uf[v]; uf[v]=r; v=t; } return r; };
        auto unite = [&](int a, int b){ int ra=findp(a), rb=findp(b); if(ra!=rb) uf[rb]=ra; };
        for(int i=0;i<N;++i){ if(isFrontier[i]) uf[i]=i; }
        for(int idx=0; idx<N; ++idx){ if(cell_number(cells[idx])>=0){
            int ulist[8]; int uc=0;
            for(int k=0;k<neighborCounts[idx];++k){ int nb=neighbors[idx][k]; if(isFrontier[nb]) ulist[uc++]=nb; }
            for(int a=1;a<uc;++a) unite(ulist[0], ulist[a]);
        } }
        std::vector<int> rootToCid(N, -1);
        int compCount=0;
        for(int i=0;i<N;++i){
            if(isFrontier[i]){
                int r = findp(i);
                int mapped = rootToCid[r];
                if(mapped==-1){ mapped = compCount; rootToCid[r] = compCount; compCount++; }
                compId[i] = mapped;
            }
        }

        ov.mineProbability.assign(N, -1.0);
        // track which probabilities are logically certain (exactly 0 or 1) via exact enumeration
        std::vector<uint8_t> probCertain(N, 0);

        struct Con { int num; int knownMines; std::vector<int> uidx; };
        struct CompData {
            std::vector<int> U;
            std::vector<Con> cons;
            int m = 0;
            std::vector<long double> waysK;
            std::vector<std::vector<long double>> mineWaysPerCellK;
            std::vector<std::vector<long double>> safeWaysPerCellK;
            long double totalSolutions = 0.0L;
            bool enumerated = false;
        };
        std::vector<CompData> comps(compCount);
        std::vector<int> localIndex(N, -1);
        for(int i=0;i<N;++i){ if(compId[i]>=0){
            auto& U = comps[compId[i]].U;
            localIndex[i] = static_cast<int>(U.size());
            U.push_back(i);
        } }
        int fixedBeforeSearch=0, freeBeforeSearch=0;
        for(int i=0;i<N;++i) {
            if(cells[i]==CellState::Mine || ov.marks[i]==Mark::Mine) ++fixedBeforeSearch;
            else if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Safe && !isFrontier[i]) ++freeBeforeSearch;
        }
        for(int idx=0; idx<N; ++idx){ if(cell_number(cells[idx])>=0){
            int uc=0; int tmp[8]; int known=0;
            for(int k=0;k<neighborCounts[idx];++k){ int nb=neighbors[idx][k]; if(cells[nb]==CellState::Mine || ov.marks[nb]==Mark::Mine){ known++; continue; } if(isFrontier[nb]) tmp[uc++]=nb; }
            if(uc==0) continue;
            int cid = compId[tmp[0]];
            Con con; con.num = cell_number(cells[idx]); con.knownMines = known;
            for(int t=0;t<uc;++t) con.uidx.push_back(localIndex[tmp[t]]);
            comps[cid].cons.push_back(std::move(con));
        } }

	const int ENUM_MAX_VARS = 128;
	const size_t ENUM_NODE_BUDGET = 8000000;

        struct CompResult {
            bool enumerated = false;
            bool inconsistent = false;
            std::vector<long double> waysK;
            std::vector<std::vector<long double>> mineWaysPerCellK;
            std::vector<std::vector<long double>> safeWaysPerCellK;
            long double totalSolutions = 0.0L;
            std::vector<double> beliefs; // size m, [-1 if unknown]
            std::vector<uint8_t> probCertainBits; // size m
            std::vector<int> forcedSafe; // global indices
            std::vector<int> forcedMine; // global indices
        };

        std::vector<CompResult> compResults(compCount);

        auto process_component = [&](int cid){
            CompResult R;
            if(cancelled()) return;
            auto& C = comps[cid];
            C.m = (int)C.U.size();
            if(C.m==0 || (int)C.cons.size()==0){ compResults[cid] = std::move(R); return; }

            // degree order once
            std::vector<int> degree(C.m, 0);
            std::vector<std::vector<int>> varToCons(C.m);
            for(size_t ci=0; ci<C.cons.size(); ++ci){ const auto& con = C.cons[ci];
                for(int ui : con.uidx){ if(ui>=0 && ui<C.m){ degree[ui]++; varToCons[ui].push_back((int)ci); } }
            }
            std::vector<int> order(C.m); for(int i=0;i<C.m;++i) order[i]=i;
            std::sort(order.begin(), order.end(), [&](int a, int b){ if(degree[a]!=degree[b]) return degree[a]>degree[b]; return a<b; });

			bool needApprox = true;
			if(C.m <= ENUM_MAX_VARS){
                // Cells touching identical clues are interchangeable. Search the
                // number of mines in each group, weighted by binomial choices.
                struct Group { std::vector<int> vars; std::vector<int> cons; std::vector<long double> choose; };
                std::vector<Group> groups;
                std::vector<int> groupedOrder(C.m);
                for(int t=0;t<C.m;++t) groupedOrder[t]=t;
                // Keep bounded search order identical across standard libraries.
                std::sort(groupedOrder.begin(),groupedOrder.end(),[&](int a,int b){
                    return varToCons[a]!=varToCons[b] ? varToCons[a]<varToCons[b] : a<b;
                });
                for(int t : groupedOrder) {
                    if(groups.empty() || groups.back().cons!=varToCons[t]) groups.push_back({{},varToCons[t],{}});
                    groups.back().vars.push_back(t);
                }
                std::sort(groups.begin(),groups.end(),[](const Group& a,const Group& b){
                    if(a.cons.size()!=b.cons.size()) return a.cons.size()>b.cons.size();
                    return a.vars.front()<b.vars.front();
                });
                for(auto& group : groups) {
                    const int size=static_cast<int>(group.vars.size());
                    group.choose.assign(size+1,1);
                    for(int k=1;k<size;++k) group.choose[k]=group.choose[k-1]*(size-k+1)/k;
                }
                R.totalSolutions=0;
                R.waysK.assign(C.m+1,0);
                R.mineWaysPerCellK.assign(C.m,std::vector<long double>(C.m+1,0));
                R.safeWaysPerCellK.assign(C.m,std::vector<long double>(C.m+1,0));
                std::vector<int> need(C.cons.size()), left(C.cons.size()), assignment(groups.size());
                for(size_t ci=0;ci<C.cons.size();++ci) {
                    need[ci]=C.cons[ci].num-C.cons[ci].knownMines;
                    left[ci]=static_cast<int>(C.cons[ci].uidx.size());
                }
                // Give larger frontiers a bounded exact attempt before sampling.
                size_t nodeBudget=C.m<=32 ? ENUM_NODE_BUDGET : 100000;
                bool aborted=false;
                std::function<void(size_t,int,long double)> dfs=[&](size_t depth,int mines,long double weight) {
                    if(aborted) return;
                    if(nodeBudget==0 || ((nodeBudget&1023)==0 && !searchBlock())) { aborted=true; return; }
                    --nodeBudget;
                    if(depth==groups.size()) {
                        R.totalSolutions+=weight;
                        R.waysK[mines]+=weight;
                        for(size_t g=0;g<groups.size();++g) {
                            const int size=static_cast<int>(groups[g].vars.size()), count=assignment[g];
                            const long double mined=weight*count/size, safe=weight*(size-count)/size;
                            for(int t : groups[g].vars) {
                                R.mineWaysPerCellK[t][mines]+=mined;
                                R.safeWaysPerCellK[t][mines]+=safe;
                            }
                        }
                        return;
                    }
                    const auto& group=groups[depth];
                    const int size=static_cast<int>(group.vars.size());
                    int low=0,high=size;
                    for(int ci : group.cons) {
                        low=std::max(low,need[ci]-(left[ci]-size));
                        high=std::min(high,need[ci]);
                    }
                    for(int count=low;count<=high && !aborted;++count) {
                        assignment[depth]=count;
                        for(int ci : group.cons) { need[ci]-=count; left[ci]-=size; }
                        dfs(depth+1,mines+count,weight*group.choose[count]);
                        for(int ci : group.cons) { need[ci]+=count; left[ci]+=size; }
                    }
                };
                dfs(0,0,1);
				if(!aborted && R.totalSolutions == 0.0L){ R.inconsistent = true; compResults[cid] = std::move(R); return; }
                if(!aborted && R.totalSolutions > 0.0L){
					R.enumerated = true;
					needApprox = false;
					R.beliefs.assign(C.m, 0.5);
					R.probCertainBits.assign(C.m, 0);
					for(int t=0;t<C.m;++t){
						long double mw=0.0L, sw=0.0L;
                        for(int k=0;k<=C.m;++k) { mw+=R.mineWaysPerCellK[t][k]; sw+=R.safeWaysPerCellK[t][k]; }
						if(mw <= 0.0L){ R.beliefs[t] = 0.0; R.probCertainBits[t]=1; }
						else if(sw <= 0.0L){ R.beliefs[t] = 1.0; R.probCertainBits[t]=1; }
						else { R.beliefs[t] = (double)(mw / (mw+sw)); }
					}
				} else {
					R.enumerated = false;
					R.totalSolutions = 0.0L;
					R.waysK.clear();
					R.mineWaysPerCellK.clear();
					R.safeWaysPerCellK.clear();
					R.probCertainBits.clear();
				}
			}

		if(needApprox){
			const int MC_MAX_VARS = 64;
			const int MC_MAX_SAMPLES = 768;
			const int MC_MIN_VALID = 64;
			const int MC_MAX_ATTEMPTS = MC_MAX_SAMPLES * 4;
			const size_t MC_NODE_BUDGET = 200000;
			if(C.m <= MC_MAX_VARS){
				std::vector<int> sampleOrder = order;
				std::vector<uint8_t> assignMC(C.m, 255);
				std::vector<int> mineHits(C.m, 0);
				int samples = 0;
				int attempts = 0;
				bool unsatDetected = false;
				bool budgetExceeded = false;
                size_t sampleBudget = MC_NODE_BUDGET;
				auto constraint_ok = [&](int conIdx)->bool{
					const auto& con = C.cons[conIdx];
					int need = con.num - con.knownMines;
					int placed = 0, unk = 0;
					for(int ui : con.uidx){ uint8_t v = assignMC[ui]; if(v==1) placed++; else if(v==255) unk++; }
					return !(need < placed || need > placed + unk);
				};
				auto consistent_after = [&](int varIdx)->bool{
					for(int conIdx : varToCons[varIdx]){ if(!constraint_ok(conIdx)) return false; }
					return true;
				};
				std::uniform_int_distribution<int> dist01(0,1);
				std::function<bool(int,size_t&)> dfsSample = [&](int depth, size_t& budget)->bool{
					if(depth == C.m) return true;
					if(budget == 0 || ((budget & 1023) == 0 && !searchBlock())){ budgetExceeded = true; return false; }
					--budget;
					int var = sampleOrder[depth];
					std::array<int,2> choices{}; int choiceCount = 0;
					for(int val=0; val<=1; ++val){
						assignMC[var] = static_cast<uint8_t>(val);
						if(consistent_after(var)) choices[choiceCount++] = val;
					}
					assignMC[var] = 255;
					if(choiceCount == 0) return false;
					if(choiceCount == 2 && dist01(rng)) std::swap(choices[0], choices[1]);
					for(int i=0;i<choiceCount;++i){
						assignMC[var] = static_cast<uint8_t>(choices[i]);
						if(dfsSample(depth+1, budget)) return true;
					}
					assignMC[var] = 255;
					return false;
				};
				while(samples < MC_MAX_SAMPLES && attempts < MC_MAX_ATTEMPTS && !budgetExceeded){
					attempts++;
						std::fill(assignMC.begin(), assignMC.end(), uint8_t{255});
					std::shuffle(sampleOrder.begin(), sampleOrder.end(), rng);
					if(dfsSample(0, sampleBudget)){
						samples++;
						for(int i=0;i<C.m;++i){ if(assignMC[i]==1) mineHits[i]++; }
					} else {
						if(!budgetExceeded) unsatDetected = true;
						break;
					}
				}
				if(samples >= MC_MIN_VALID){
					R.beliefs.assign(C.m, 0.5);
					R.probCertainBits.assign(C.m, 0);
					for(int i=0;i<C.m;++i){
						// Samples never prove a cell; smooth endpoints away from certainty.
                        R.beliefs[i] = (mineHits[i] + 1.0) / (samples + 2.0);
					}
					needApprox = false;
				}
				if(unsatDetected){ R.inconsistent = true; compResults[cid] = std::move(R); return; }
                if(budgetExceeded){
					needApprox = true;
				}
			} else {
				needApprox = true;
			}
		}
		if(needApprox){
				struct EdgeRef { int conIndex; int posInCon; };
				std::vector<std::vector<EdgeRef>> varEdges(C.m);
				for(int j=0;j<(int)C.cons.size();++j){ const auto& con=C.cons[j]; for(int p=0;p<(int)con.uidx.size();++p){ int v=con.uidx[p]; varEdges[v].push_back({j,p}); } }
                // A clue has at most eight edges: avoid the dense variables x clues table.
                std::vector<std::vector<int>> edgeIndexLookup(C.cons.size());
                for(size_t j=0;j<C.cons.size();++j) edgeIndexLookup[j].resize(C.cons[j].uidx.size());
                for(int v=0;v<C.m;++v) for(int e=0;e<(int)varEdges[v].size();++e) {
                    const auto& ref=varEdges[v][e];
                    edgeIndexLookup[ref.conIndex][ref.posInCon]=e;
                }
				std::vector<std::vector<double>> vToF(C.m), fToV(C.cons.size());
				for(int v=0; v<C.m; ++v){ vToF[v].assign(varEdges[v].size(), 0.5); }
				for(size_t j=0;j<C.cons.size();++j){ fToV[j].assign(C.cons[j].uidx.size(), 0.5); }
				const int BP_MAX_ITERS = 20; const double BP_DAMP = 0.5; const double BP_EPS = 1e-6;
				int maxConSize = 0; for(const auto& con : C.cons){ maxConSize = std::max(maxConSize, (int)con.uidx.size()); }
				std::vector<long double> dpBuf(maxConSize+1), ndpBuf(maxConSize+1);
				int maxDeg = 0; for(int v=0;v<C.m;++v){ maxDeg = std::max(maxDeg, (int)varEdges[v].size()); }
				std::vector<double> pref1Buf(maxDeg+1), pref0Buf(maxDeg+1), suf1Buf(maxDeg+1), suf0Buf(maxDeg+1);
				for(int it=0; it<BP_MAX_ITERS; ++it){
                    if(cancelled()) return;
					for(int j=0;j<(int)C.cons.size(); ++j){
						const auto& con = C.cons[j];
						const int sz = (int)con.uidx.size();
						const int need = con.num - con.knownMines;
						for(int p=0; p<sz; ++p){
							for(int k=0;k<sz;++k) dpBuf[k]=0.0L;
							dpBuf[0]=1.0L;
							for(int q=0; q<sz; ++q){ if(q==p) continue; int vv=con.uidx[q]; int eidx=edgeIndexLookup[j][q]; double prob=0.5; if(eidx>=0) prob=vToF[vv][eidx];
								for(int k=0;k<sz;++k) ndpBuf[k]=0.0L;
								for(int k=0;k<sz-1;++k){ if(dpBuf[k]==0.0L) continue; ndpBuf[k] += dpBuf[k]*(1.0L-(long double)prob); ndpBuf[k+1]+=dpBuf[k]*(long double)prob; }
								std::swap(dpBuf,ndpBuf);
							}
							long double A=0.0L,B=0.0L; if(need-1>=0 && need-1<=sz-1) A=dpBuf[need-1]; if(need>=0 && need<=sz-1) B=dpBuf[need];
							double msg=0.5; if(A==0.0L && B==0.0L) msg=0.5; else if(A==0.0L) msg=0.0; else if(B==0.0L) msg=1.0; else { double a=(double)A,b=(double)B; msg=a/(a+b);}
							msg = BP_DAMP * fToV[j][p] + (1.0 - BP_DAMP) * std::min(1.0-1e-9, std::max(1e-9, msg));
							fToV[j][p] = msg;
						}
					}
					for(int v=0; v<C.m; ++v){ int deg=(int)varEdges[v].size(); if(deg==0) continue;
						pref1Buf[0]=1.0; pref0Buf[0]=1.0;
						for(int i=0;i<deg;++i){ const auto&e=varEdges[v][i]; double m1=fToV[e.conIndex][e.posInCon]; pref1Buf[i+1]=pref1Buf[i]*m1; pref0Buf[i+1]=pref0Buf[i]*(1.0-m1);}
						suf1Buf[deg]=1.0; suf0Buf[deg]=1.0;
						for(int i=deg-1;i>=0;--i){ const auto&e=varEdges[v][i]; double m1=fToV[e.conIndex][e.posInCon]; suf1Buf[i]=suf1Buf[i+1]*m1; suf0Buf[i]=suf0Buf[i+1]*(1.0-m1);}
						for(int i=0;i<deg;++i){ double p1=pref1Buf[i]*suf1Buf[i+1]; double p0=pref0Buf[i]*suf0Buf[i+1]; double msg=(p1==0.0 && p0==0.0)?0.5:(p1/(p1+p0)); msg=std::min(1.0-BP_EPS, std::max(BP_EPS, msg)); vToF[v][i] = BP_DAMP * vToF[v][i] + (1.0 - BP_DAMP) * msg; } }
				}
				R.beliefs.assign(C.m, 0.5);
				for(int v=0; v<C.m; ++v){ int deg=(int)varEdges[v].size(); double prod1=1.0,prod0=1.0; for(int i=0;i<deg;++i){ const auto&e=varEdges[v][i]; double m1=fToV[e.conIndex][e.posInCon]; prod1*=m1; prod0*=(1.0-m1);} double p = (prod1==0.0&&prod0==0.0)?0.5:(prod1/(prod1+prod0)); p=std::min(1.0-BP_EPS, std::max(BP_EPS, p)); R.beliefs[v]=p; }
			}

            if(C.m>0 && C.m<=128 && !C.cons.empty() && !R.enumerated && !cancelled()){
                size_t budget = 500000; // Shared by all probes in this component.
                std::vector<uint8_t> assignment(C.m, 255);
                auto consistent = [&](int var) {
                    for(int ci : varToCons[var]) {
                        const auto& con=C.cons[ci];
                        int placed=0, unknown=0;
                        for(int v : con.uidx) {
                            placed += assignment[v]==1;
                            unknown += assignment[v]==255;
                        }
                        const int need=con.num-con.knownMines;
                        if(need<placed || need>placed+unknown) return false;
                    }
                    return true;
                };
                std::function<Sat(int)> search = [&](int depth) -> Sat {
                    if(depth==C.m) return Sat::Yes;
                    if(budget==0 || ((budget & 1023)==0 && !searchBlock())) return Sat::Unknown;
                    --budget;
                    const int var=order[depth];
                    if(assignment[var]!=255) return search(depth+1);
                    bool incomplete=false;
                    for(uint8_t value : {uint8_t(0), uint8_t(1)}) {
                        assignment[var]=value;
                        const Sat result=consistent(var) ? search(depth+1) : Sat::No;
                        assignment[var]=255;
                        if(result==Sat::Yes) return result;
                        incomplete |= result==Sat::Unknown;
                        if(budget==0) return Sat::Unknown;
                    }
                    return incomplete ? Sat::Unknown : Sat::No;
                };
                for(int t=0;t<C.m && budget>0 && !cancelled();++t) {
                    assignment[t]=0;
                    const Sat safe=consistent(t) ? search(0) : Sat::No;
                    assignment[t]=1;
                    const Sat mine=consistent(t) ? search(0) : Sat::No;
                    assignment[t]=255;
                    if(safe==Sat::No && mine==Sat::No) { R.inconsistent=true; break; }
                    if(safe==Sat::Yes && mine==Sat::No) R.forcedSafe.push_back(C.U[t]);
                    if(safe==Sat::No && mine==Sat::Yes) R.forcedMine.push_back(C.U[t]);
                }
            }

            compResults[cid] = std::move(R);
        };

        int threadCount = threads;
        if(threadCount <= 0){ threadCount = (int)std::thread::hardware_concurrency(); if(threadCount<=0) threadCount = 1; }
        threadCount = std::min({threadCount, std::max(1, compCount), 8});

        if(threadCount == 1){
            for(int cid=0; cid<compCount; ++cid){ process_component(cid); }
        } else {
            WorkerPool::instance().run(threadCount, compCount, [&](int cid){ process_component(cid); });
        }

        if(cancelled()) return invalid();
        for(const auto& result : compResults) if(result.inconsistent) return invalid();

        // apply forced marks from SAT and enumerated 0/1 certainty
        bool anyForced=false;
		for(int cid=0; cid<compCount; ++cid){
			const auto& R = compResults[cid];
			for(int g : R.forcedSafe){ if(g>=0 && g<N){ if(ov.marks[g] != Mark::Safe){ ov.marks[g]=Mark::Safe; ov.hasGuaranteedSafe = true; anyForced=true; enqueueNbrNumbers(g);} } }
			for(int g : R.forcedMine){ if(g>=0 && g<N){ if(ov.marks[g] != Mark::Mine){ ov.marks[g]=Mark::Mine; anyForced=true; enqueueNbrNumbers(g);} } }
        }
        if(anyForced){
            while(!queue.empty()){
                int idxCenter = queue.back(); queue.pop_back(); inQueue[idxCenter]=0;
                int num2 = numbers[idxCenter]; if(num2 < 0) continue;
                int knownMines2 = 0; int unknownIdx2[8]; int ucount2 = 0;
                for(int k=0;k<neighborCounts[idxCenter];++k){ int nb2=neighbors[idxCenter][k]; CellState s2=cells[nb2]; if(s2==CellState::Mine || ov.marks[nb2]==Mark::Mine){ knownMines2++; continue; } if(s2==CellState::Unknown && ov.marks[nb2]!=Mark::Safe){ unknownIdx2[ucount2++]=nb2; } }
                int remaining2 = num2 - knownMines2; if(remaining2 < 0 || remaining2 > ucount2) return invalid();
			if(remaining2 == 0 && ucount2 > 0){ for(int i=0;i<ucount2;++i){ int u=unknownIdx2[i]; if(ov.marks[u] != Mark::Safe){ ov.marks[u]=Mark::Safe; ov.hasGuaranteedSafe = true; enqueueNbrNumbers(u);} } }
			else if(remaining2 == ucount2 && ucount2>0){ for(int i=0;i<ucount2;++i){ int u=unknownIdx2[i]; if(ov.marks[u] != Mark::Mine){ ov.marks[u]=Mark::Mine; enqueueNbrNumbers(u);} } }
            }
        }

        // expose enumerated results to comps for global combination step
        for(int cid=0; cid<compCount; ++cid){
            auto& C = comps[cid];
            auto& R = compResults[cid];
            if(R.enumerated){
                C.enumerated = true;
                C.waysK = std::move(R.waysK);
                C.mineWaysPerCellK = std::move(R.mineWaysPerCellK);
                C.safeWaysPerCellK = std::move(R.safeWaysPerCellK);
                C.totalSolutions = R.totalSolutions;
            }
        }

		for(int cid=0; cid<compCount; ++cid){
			const auto& C = comps[cid];
			const auto& R = compResults[cid];
			if(C.m==0 || R.beliefs.empty()) continue;
			for(int t=0; t<C.m; ++t){
				int g = C.U[t];
				if(g>=0 && g<N){
					ov.mineProbability[g] = R.beliefs[t];
					if(t < (int)R.probCertainBits.size() && R.probCertainBits[t]) probCertain[g] = 1;
				}
			}
		}

        for(int i=0;i<N;++i) {
            if(cells[i]==CellState::Mine || ov.marks[i]==Mark::Mine) ov.mineProbability[i]=1;
            else if(cells[i]!=CellState::Unknown || ov.marks[i]==Mark::Safe) ov.mineProbability[i]=0;
        }

        int enumeratedVars = 0;
        bool allEnumerated = true;
        for(const auto& C : comps) { allEnumerated &= C.enumerated; enumeratedVars += C.m; }
        // Only a completely enumerated frontier supports exact conditioning.
        if(totalMines>=0 && allEnumerated) {
            const int remaining=totalMines-fixedBeforeSearch;
            const int freeCount=freeBeforeSearch;
            if(remaining<0 || remaining>enumeratedVars+freeCount) return invalid();
            const long double negInf=-std::numeric_limits<long double>::infinity();
            auto logAdd = [&](long double a, long double b) {
                if(a==negInf) return b;
                if(b==negInf) return a;
                if(a<b) std::swap(a,b);
                return a+std::log1p(std::exp(b-a));
            };
            // Compact away impossible counts. A component with a fixed mine
            // count occupies one slot, irrespective of its number of cells.
            struct Distribution { int offset=0; std::vector<long double> values; };
            std::vector<Distribution> logs(compCount), prefix(compCount+1);
            for(int c=0;c<compCount;++c) {
                const auto& ways=comps[c].waysK;
                int first=0,last=comps[c].m;
                while(first<=last && ways[first]==0) ++first;
                while(last>=first && ways[last]==0) --last;
                logs[c].offset=first;
                for(int k=first;k<=last;++k) logs[c].values.push_back(ways[k]>0 ? std::log(ways[k]) : negInf);
            }
            // Bound actual support width and work, rather than rejecting every
            // frontier above an arbitrary number of cells.
            size_t entries=1,work=0,width=1;
            bool withinBudget=true;
            for(const auto& distribution : logs) {
                work+=width*distribution.values.size();
                width+=distribution.values.size()-1;
                entries+=width;
                if(entries>1000000 || work>8000000) { withinBudget=false; break; }
            }
            if(withinBudget) {
            prefix[0].values={0};
            for(int c=0;c<compCount;++c) {
                if(cancelled()) return invalid();
                const auto& a=prefix[c]; const auto& b=logs[c];
                auto& out=prefix[c+1]; out.offset=a.offset+b.offset;
                out.values.assign(a.values.size()+b.values.size()-1,negInf);
                for(size_t i=0;i<a.values.size();++i) {
                    if((i&127)==0 && cancelled()) return invalid();
                    if(a.values[i]!=negInf) for(size_t j=0;j<b.values.size();++j) if(b.values[j]!=negInf)
                        out.values[i+j]=logAdd(out.values[i+j],a.values[i]+b.values[j]);
                }
            }
            std::vector<long double> freeLog(freeCount+1, 0);
            for(int k=1;k<=freeCount;++k)
                freeLog[k]=freeLog[k-1]+std::log(static_cast<long double>(freeCount-k+1))-std::log(static_cast<long double>(k));
            auto outside = [&](int k) { return k>=0 && k<=freeCount ? freeLog[k] : negInf; };
            long double total=negInf;
            std::vector<long double> tail(prefix.back().values.size(),negInf);
            for(size_t k=0;k<tail.size();++k) {
                tail[k]=outside(remaining-prefix.back().offset-static_cast<int>(k));
                if(tail[k]!=negInf && prefix.back().values[k]!=negInf)
                    total=logAdd(total,prefix.back().values[k]+tail[k]);
            }
            if(total==negInf) return invalid();
            // Backward messages reuse the prefix convolution. Re-convolving all
            // other components for every cell group grows quadratically in the
            // number of components and is unnecessary.
            for(int c=compCount-1;c>=0;--c) {
                if(cancelled()) return invalid();
                const auto& C=comps[c];
                std::vector<long double> weights(C.m+1,negInf);
                std::vector<long double> previous(prefix[c].values.size(),negInf);
                for(size_t p=0;p<previous.size();++p) {
                    if((p&127)==0 && cancelled()) return invalid();
                    for(size_t k=0;k<logs[c].values.size();++k) {
                        if(tail[p+k]==negInf || logs[c].values[k]==negInf) continue;
                        previous[p]=logAdd(previous[p],logs[c].values[k]+tail[p+k]);
                        if(prefix[c].values[p]!=negInf) {
                            const int count=logs[c].offset+static_cast<int>(k);
                            weights[count]=logAdd(weights[count],prefix[c].values[p]+tail[p+k]);
                        }
                    }
                }
                tail=std::move(previous);
                for(int t=0;t<C.m;++t) {
                    long double mine=negInf, safe=negInf;
                    for(int k=0;k<=C.m;++k) if(weights[k]!=negInf) {
                        long double m=C.mineWaysPerCellK[t][k], f=C.safeWaysPerCellK[t][k];
                        if(m>0) mine=logAdd(mine,std::log(m)+weights[k]);
                        if(f>0) safe=logAdd(safe,std::log(f)+weights[k]);
                    }
                    const int g=C.U[t];
                    probCertain[g]=(mine==negInf || safe==negInf);
                    ov.mineProbability[g]=mine==negInf ? 0 : safe==negInf ? 1 :
                        static_cast<double>(std::exp(mine-logAdd(mine,safe)));
                }
            }
            if(freeCount>0) {
                long double mine=negInf, safe=negInf;
                for(size_t k=0;k<prefix.back().values.size();++k) {
                    const int out=remaining-prefix.back().offset-static_cast<int>(k);
                    if(out<0 || out>freeCount || prefix.back().values[k]==negInf) continue;
                    const long double weight=prefix.back().values[k]+freeLog[out];
                    if(out>0) mine=logAdd(mine,weight+std::log(static_cast<long double>(out)/freeCount));
                    if(out<freeCount) safe=logAdd(safe,weight+std::log(static_cast<long double>(freeCount-out)/freeCount));
                }
                const double p=mine==negInf ? 0 : safe==negInf ? 1 : static_cast<double>(std::exp(mine-logAdd(mine,safe)));
                for(int i=0;i<N;++i) if(cells[i]==CellState::Unknown && !isFrontier[i] && ov.marks[i]==Mark::None) {
                    ov.mineProbability[i]=p; probCertain[i]=(mine==negInf || safe==negInf);
                }
            }
            }
        }

        std::vector<double> sumRisk(N, 0.0); std::vector<int> countRisk(N, 0);
        for(int i=0;i<N;++i){
            if(cell_number(cells[i])>=0){
                int knownMines=0; int ucount=0; int uidx[8];
                for(int k=0;k<neighborCounts[i];++k){
                    int nb=neighbors[i][k]; CellState s=cells[nb];
                    if(s==CellState::Mine||ov.marks[nb]==Mark::Mine){ knownMines++; }
                    else if(s==CellState::Unknown && ov.marks[nb]!=Mark::Safe){ uidx[ucount++]=nb; }
                }
                int remaining = cell_number(cells[i]) - knownMines;
                if(ucount>0 && remaining>0){
                    double contrib=(double)remaining/(double)ucount;
                    for(int t=0;t<ucount;++t){ sumRisk[uidx[t]]+=contrib; countRisk[uidx[t]]++; }
                }
            }
        }

		// complete probabilities for all unknown cells if totalMines is known
		if(totalMines>=0){
            int knownGlobal=0; for(int i=0;i<N;++i){ if(cells[i]==CellState::Mine || ov.marks[i]==Mark::Mine) knownGlobal++; }
            int remainingGlobal = std::max(0, totalMines - knownGlobal);

            int freeUnknown=0; long double expectedEnumerated=0.0L;
            for(int i=0;i<N;++i){
                if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Safe && ov.marks[i]!=Mark::Mine){
                    if(ov.mineProbability[i] >= 0.0){ expectedEnumerated += ov.mineProbability[i]; }
                    else { freeUnknown++; }
                }
            }
            if(freeUnknown>0){
                double expectedOutside = (double)remainingGlobal - (double)expectedEnumerated;
                if(expectedOutside < 0.0) expectedOutside = 0.0;
                if(expectedOutside > (double)freeUnknown) expectedOutside = (double)freeUnknown;
                double p_free = expectedOutside / (double)freeUnknown;
                for(int i=0;i<N;++i){
                    if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Safe && ov.marks[i]!=Mark::Mine && ov.mineProbability[i] < 0.0){
                        ov.mineProbability[i] = p_free;
                    }
                }
            }
        } else {
            for(int i=0;i<N;++i){
                if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Safe && ov.marks[i]!=Mark::Mine && ov.mineProbability[i] < 0.0){
                    double r = countRisk[i]>0 ? (sumRisk[i]/(double)countRisk[i]) : -1.0;
                    ov.mineProbability[i] = r;
                }
            }
        }

        // clamp near-certain probabilities exactly to 0/1 to harden certainty detection
        if(ov.mineProbability.size()==(size_t)N){
            for(int i=0;i<N;++i){ if(ov.mineProbability[i] >= 0.0){ if(ov.mineProbability[i] <= 1e-12) ov.mineProbability[i]=0.0; else if(ov.mineProbability[i] >= 1.0 - 1e-12) ov.mineProbability[i]=1.0; } }
        }
		const double kCascadeDampen = 0.5;
		auto probOf = [&](int idx)->double{
			double p = 1.0;
			if(idx>=0 && idx<N){
				if(cells[idx]==CellState::Unknown){
					p = ov.mineProbability[idx];
					if(!(p>=0.0 && p<=1.0)) p = 1.0;
				} else if(cells[idx]==CellState::Mine || ov.marks[idx]==Mark::Mine){
					p = 1.0;
				} else {
					p = 0.0;
				}
			}
			if(p < 0.0) p = 0.0;
			if(p > 1.0) p = 1.0;
			return p;
		};
		bool hasProb=false;
		if(ov.mineProbability.size()==(size_t)N){
			for(int i=0;i<N;++i){ if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Mine && ov.mineProbability[i]>=0.0){ hasProb=true; break; } }
		}
		std::vector<int> changedByProb; changedByProb.reserve(N/8+1);
		if(hasProb){
            // only convert to Safe/Mine when probabilities are logically certain
            for(int i=0;i<N;++i){
                if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Mine && ov.mineProbability[i]>=0.0 && probCertain[i]){
                    if(ov.mineProbability[i] == 0.0){ if(ov.marks[i]!=Mark::Safe){ ov.marks[i]=Mark::Safe; ov.hasGuaranteedSafe = true; changedByProb.push_back(i); } }
                    else if(ov.mineProbability[i] == 1.0){ if(ov.marks[i]!=Mark::Mine){ ov.marks[i]=Mark::Mine; changedByProb.push_back(i); } }
                }
            }
        // if new certain marks were added via probabilities, re-run local propagation around them
        if(!changedByProb.empty()){
            for(int idxChanged : changedByProb){ enqueueNbrNumbers(idxChanged); }
            while(!queue.empty()){
                int idxCenter = queue.back(); queue.pop_back(); inQueue[idxCenter]=0;
                int y2 = idxCenter / w; (void)y2; int x2 = idxCenter % w; (void)x2;
                int num2 = numbers[idxCenter];
                if(num2 < 0) continue;
                int knownMines2 = 0;
                int unknownIdx2[8]; int ucount2 = 0;
                for(int k=0;k<neighborCounts[idxCenter];++k){
                    int nb2 = neighbors[idxCenter][k];
                    CellState s2 = cells[nb2];
                    // only count preexisting flags or deduced mines here; do not treat FlagForChord as flagged
                    if(s2 == CellState::Mine || ov.marks[nb2] == Mark::Mine){ knownMines2++; continue; }
                    if(s2 == CellState::Unknown && ov.marks[nb2] != Mark::Safe){ unknownIdx2[ucount2++] = nb2; }
                }
                int remaining2 = num2 - knownMines2;
                if(remaining2 < 0 || remaining2 > ucount2){ return invalid(); }
                if(remaining2 == 0 && ucount2 > 0){
                    for(int i=0;i<ucount2;++i){ int u=unknownIdx2[i]; if(ov.marks[u] != Mark::Safe){ ov.marks[u]=Mark::Safe; enqueueNbrNumbers(u);} }
                } else if(remaining2 == ucount2 && ucount2>0){
                    for(int i=0;i<ucount2;++i){
                        int u=unknownIdx2[i];
                        if(ov.marks[u] == Mark::FlagForChord || ov.marks[u] == Mark::FlagForChordReady){ continue; }
                        if(ov.marks[u] != Mark::Mine){ ov.marks[u] = Mark::Mine; enqueueNbrNumbers(u);}
                    }
                }
            }
            // refresh guaranteed safe after propagation
            ov.hasGuaranteedSafe = false; for(const auto m : ov.marks){ if(m==Mark::Safe){ ov.hasGuaranteedSafe=true; break; } }
        }
        for(int i=0;i<N;++i) {
            if(cells[i]==CellState::Mine || ov.marks[i]==Mark::Mine) ov.mineProbability[i]=1;
            else if(cells[i]!=CellState::Unknown || ov.marks[i]==Mark::Safe) ov.mineProbability[i]=0;
        }

			// only set Guess when no guaranteed safe exists
			if(!ov.hasGuaranteedSafe){
				const double eps = 1e-9;
				double bestScore = -std::numeric_limits<double>::infinity();
				double bestRisk = 1.0;
				std::vector<int> bestGuess;
				bestGuess.reserve(8);
				for(int i=0;i<N;++i){
					if(cells[i]!=CellState::Unknown) continue;
					if(ov.marks[i]==Mark::Mine) continue;
					double p = probOf(i);
					if(!(p>=0.0 && p<1.0)) continue;
					double safeProb = 1.0 - p;
					if(safeProb < 0.0) safeProb = 0.0;
					long double zeroProb = 1.0L;
					int unknownCount = 0;
					for(int k=0;k<neighborCounts[i];++k){
						int nb = neighbors[i][k];
						double pn = probOf(nb);
						zeroProb *= (1.0 - pn);
						if(cells[nb]==CellState::Unknown && ov.marks[nb]!=Mark::Mine){ unknownCount++; }
					}
					double cascadeBonus = (double)zeroProb * (double)unknownCount * kCascadeDampen * safeProb;
					double score = safeProb + cascadeBonus;
					if(score > bestScore + eps){
						bestScore = score;
						bestRisk = p;
						bestGuess.clear();
						bestGuess.push_back(i);
					} else if(std::abs(score - bestScore) <= eps){
						if(p < bestRisk - eps){
							bestRisk = p;
							bestGuess.clear();
							bestGuess.push_back(i);
						} else if(std::abs(p - bestRisk) <= eps){
							bestGuess.push_back(i);
						}
					}
				}
				for(int idxGuess : bestGuess){ ov.marks[idxGuess] = Mark::Guess; }
			}

    // chord valuation using probabilities and cascace heuristic
    if(enableChords && ov.mineProbability.size()==(size_t)N){
        struct ChordCandidate {
            int centerIdx;
            int missing;
            std::vector<int> unknowns;
            std::vector<int> flagList;        // all required flags (pending + already placed)
            std::vector<int> pendingFlags;    // subset that still need to be placed
            std::vector<int> readyFlags;      // subset already flagged on the live board
            int clicks;                       // flags to place (excluding already-flagged board cells) + 1 chord click
            double expectedValue;             // reveals + cascade
            double bestSingle;                // best single-click alternative near this center
            bool inProgress;                  // true if user just placed a flag enabling a break-even chord
        };

        auto computeExpected = [&](const std::vector<int>& unknowns, const std::vector<int>& plannedFlags, int centerA, int centerB, double dampen)->std::pair<double,double>{
            double expectedReveals = 0.0; double cascadeBonus = 0.0;
            for(int u : unknowns){
                bool isPlannedFlag = false; for(int f : plannedFlags){ if(f==u){ isPlannedFlag=true; break; } }
                if(isPlannedFlag) continue;
                double pu = probOf(u);
                expectedReveals += (1.0 - pu);
                int uc2=0; int list2[8];
                for(int kk=0; kk<neighborCounts[u]; ++kk){ int m = neighbors[u][kk]; if(m==centerA || m==centerB) continue; if(cells[m]==CellState::Unknown){ bool pf=false; for(int f : plannedFlags){ if(f==m){ pf=true; break; } } if(pf) continue; list2[uc2++] = m; } }
                if(uc2>0){ long double zeroProb = 1.0L; for(int t=0;t<uc2;++t){ double pm = probOf(list2[t]); zeroProb *= (long double)(1.0 - pm); } cascadeBonus += (double)zeroProb * (double)uc2 * dampen * (1.0 - pu); }
            }
            return { expectedReveals, cascadeBonus };
        };

        auto computeSingleValue = [&](int u, const std::vector<int>& plannedFlags, int centerA, int centerB, double dampen)->double{
            double pu = probOf(u);
            double singleVal = (1.0 - pu);
            int uc2=0; int list2[8];
            for(int kk=0; kk<neighborCounts[u]; ++kk){ int m = neighbors[u][kk]; if(m==centerA || m==centerB) continue; if(cells[m]==CellState::Unknown){ bool pf=false; for(int f : plannedFlags){ if(f==m){ pf=true; break; } } if(pf) continue; list2[uc2++] = m; } }
            if(uc2>0){ long double zeroProb = 1.0L; for(int t=0;t<uc2;++t){ double pm = probOf(list2[t]); zeroProb *= (long double)(1.0 - pm); } singleVal += (double)zeroProb * (double)uc2 * dampen; }
            return singleVal;
        };

        auto computeBestKSingles = [&](const std::vector<int>& unknowns, const std::vector<int>& plannedFlags, int centerA, int centerB, double dampen, int K)->double{
            if(K <= 0) return 0.0;
            struct ValIdx { double v; int idx; };
            std::vector<ValIdx> vals; vals.reserve(unknowns.size());
            for(int u : unknowns){
                bool isPlannedFlag = false; for(int f : plannedFlags){ if(f==u){ isPlannedFlag=true; break; } }
                if(isPlannedFlag) continue;
                double v = computeSingleValue(u, plannedFlags, centerA, centerB, dampen);
                vals.push_back({v, u});
            }
            if(vals.empty()) return 0.0;
            std::nth_element(vals.begin(), vals.begin() + std::min<int>(K, (int)vals.size()) - 1, vals.end(), [](const ValIdx& a, const ValIdx& b){ return a.v > b.v; });
            // nth_element already places the K best values in the prefix; their order is irrelevant.
            double sum = 0.0; int take = std::min<int>(K, (int)vals.size());
            for(int i=0;i<take;++i) sum += vals[i].v;
            return sum;
        };

        std::vector<ChordCandidate> candidates; candidates.reserve(w*h/2);
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                int idx = to_index(x,y,w);
                int num = cell_number(cells[idx]); if(num<0) continue;

                int flagged=0; std::vector<int> unknowns; unknowns.reserve(8);
                std::vector<int> needFlags; needFlags.reserve(8);
                std::vector<int> placedFlags; placedFlags.reserve(8);
                int userFlags = 0; // flags present on board but not deduced by solver
                for(int k=0;k<neighborCounts[idx];++k){
                    int nb = neighbors[idx][k];
                    CellState s = cells[nb];
                    if(s==CellState::Mine){
                        flagged++;
                        placedFlags.push_back(nb);
                        if(ov.marks[nb] != Mark::Mine) userFlags++;
                        continue;
                    }
                    if(s==CellState::Unknown){
                        unknowns.push_back(nb);
                        if(ov.marks[nb]==Mark::Mine){ needFlags.push_back(nb); }
                    }
                }
                int missing = num - flagged;
                if(missing < 0 || missing > (int)unknowns.size()) continue;

                // derive flags required (certain mines) to satisfy missing, or none if missing==0
                std::vector<int> certainFlags = needFlags;
                if(missing>0){
                    if((int)certainFlags.size() != missing) continue; // cannot chord safely
                } else {
                    // if user just placed at least one neighbor flag, allow a break-even chord with 1 unknown
                    // to persist planned action; otherwise require 2+ unknowns to ensure savings.
                    if(!((int)unknowns.size() >= 1 && userFlags > 0)){
                        if((int)unknowns.size() < 2) continue; // require savings when not in-progress
                    }
                }

                for(int fPlaced : placedFlags){
                    bool exists=false;
                    for(int v : certainFlags){ if(v==fPlaced){ exists=true; break; } }
                    if(!exists){ certainFlags.push_back(fPlaced); }
                }

                std::vector<int> pendingFlags; pendingFlags.reserve(certainFlags.size());
                std::vector<int> readyFlags; readyFlags.reserve(certainFlags.size());
                for(int f : certainFlags){
                    if(cells[f] == CellState::Mine){
                        readyFlags.push_back(f);
                    } else {
                        pendingFlags.push_back(f);
                    }
                }

                int safeReveals = (int)unknowns.size() - (int)pendingFlags.size();
                if(safeReveals <= 0){
                    continue;
                }

                int clicks = (int)pendingFlags.size() + 1;
                // Keep break-even candidates for in-progress chords and shared-flag pairs.

                auto ev = computeExpected(unknowns, certainFlags, idx, -1, kCascadeDampen);
                double totalValue = ev.first + ev.second;
                int budgetK = clicks; // alternatively spend the same number of clicks on best singles
                double bestAltSingles = computeBestKSingles(unknowns, certainFlags, idx, -1, kCascadeDampen, budgetK);

                bool inProgress = (missing==0 && (int)unknowns.size()==1 && userFlags>0);
                candidates.push_back(ChordCandidate{
                    idx,
                    missing,
                    std::move(unknowns),
                    std::move(certainFlags),
                    std::move(pendingFlags),
                    std::move(readyFlags),
                    clicks,
                    totalValue,
                    bestAltSingles,
                    inProgress
                });
            }
        }

        // find synergistic pairs that can share flags and prefer pairs with strong positive net improvement
        struct PairPick { int a; int b; int clicks; double value; double bestTwoSingles; double improvement; std::vector<int> unionFlags; std::vector<int> unionUnknowns; };
        std::vector<PairPick> pairPicks; pairPicks.reserve(candidates.size());

        auto unionUnique = [&](const std::vector<int>& A, const std::vector<int>& B){
            std::vector<int> out = A; out.reserve(A.size()+B.size());
            for(int v : B){ bool found=false; for(int u : out){ if(u==v){ found=true; break; } } if(!found) out.push_back(v); }
            return out;
        };

        const double kMargin = 1e-3;

        std::vector<int> candGrid(N, -1);
        for(int ci=0;ci<(int)candidates.size();++ci) candGrid[candidates[ci].centerIdx]=ci;

        auto tryPair = [&](int i, int j){
            if(i >= j) return;
            const auto& A = candidates[i]; const auto& B = candidates[j];
            std::vector<int> flagUnion = unionUnique(A.flagList, B.flagList);
            int flagsToPlace=0; for(int f : flagUnion){ if(cells[f] != CellState::Mine) flagsToPlace++; }
            int pairClicks = flagsToPlace + 2;
            std::vector<int> unknownUnion = unionUnique(A.unknowns, B.unknowns);
            auto evPair = computeExpected(unknownUnion, flagUnion, A.centerIdx, B.centerIdx, kCascadeDampen);
            double pairValue = evPair.first + evPair.second;
            double bestKSinglesPair = computeBestKSingles(unknownUnion, flagUnion, A.centerIdx, B.centerIdx, kCascadeDampen, pairClicks);
            double improvement = pairValue - bestKSinglesPair;
            if(improvement > kMargin){
                pairPicks.push_back(PairPick{ i, j, pairClicks, pairValue, bestKSinglesPair, improvement, std::move(flagUnion), std::move(unknownUnion) });
            }
        };

        // Any shared adjacent cell places the two centers within two cells.
        // Each center occupies one grid slot, so neither hashing nor deduplication is needed.
        for(int ci=0;ci<(int)candidates.size();++ci) {
            int cx=candidates[ci].centerIdx%w, cy=candidates[ci].centerIdx/w;
            for(int dy=-2;dy<=2;++dy) for(int dx=-2;dx<=2;++dx) {
                int nx=cx+dx, ny=cy+dy;
                if(!in_bounds(nx,ny,w,h)) continue;
                const int cj=candGrid[ny*w+nx];
                if(cj>ci) tryPair(ci,cj);
            }
        }

        // greedy select non-overlapping best pairs
        std::vector<uint8_t> picked(candidates.size(), 0);
        std::vector<uint8_t> selected(candidates.size(), 0);
        std::sort(pairPicks.begin(), pairPicks.end(), [](const PairPick& a, const PairPick& b){ return a.improvement > b.improvement; });
        for(const auto& p : pairPicks){
            if(picked[p.a] || picked[p.b]) continue;
            picked[p.a]=picked[p.b]=1;
            selected[p.a]=selected[p.b]=1;
        }

        // always keep in-progress single chords visible (user recently placed enabling flag)
        for(int i=0;i<(int)candidates.size();++i){
            if(picked[i]) continue;
            const auto& C = candidates[i];
            if(C.inProgress){
                picked[i]=1;
                selected[i]=1;
            }
        }

        // pick remaining profitable singles
        for(int i=0;i<(int)candidates.size();++i){
            if(picked[i]) continue;
            const auto& C = candidates[i];
            // baseline: same click budget on best singles
            double improvement = C.expectedValue - C.bestSingle;
            if(improvement > kMargin){
                picked[i]=1;
                selected[i]=1;
            }
        }

        if(!selected.empty()){
            const uint8_t kReadyBit = 0x1;
            const uint8_t kPendingBit = 0x2;
            std::vector<uint8_t> chordMask(N, 0);
            std::vector<uint8_t> flagMask(N, 0);

            for(size_t i=0; i<candidates.size(); ++i){
                if(!selected[i]) continue;
                const auto& C = candidates[i];
                bool chordReady = !C.flagList.empty() && C.pendingFlags.empty();
                if(chordReady){
                    chordMask[C.centerIdx] |= kReadyBit;
                } else {
                    chordMask[C.centerIdx] |= kPendingBit;
                }
                for(int f : C.readyFlags){
                    flagMask[f] |= chordReady ? kReadyBit : kPendingBit;
                }
                for(int f : C.pendingFlags){
                    flagMask[f] |= kPendingBit;
                }
            }

            for(int idx=0; idx<N; ++idx){
                uint8_t mask = chordMask[idx];
                if(mask == 0) continue;
                bool pending = (mask & kPendingBit) != 0;
                bool ready = (mask & kReadyBit) != 0;
                if(pending){
                    if(ov.marks[idx]==Mark::None || ov.marks[idx]==Mark::Chord || ov.marks[idx]==Mark::ChordReady){
                        ov.marks[idx] = Mark::Chord;
                    }
                } else if(ready){
                    if(ov.marks[idx]==Mark::None || ov.marks[idx]==Mark::Chord || ov.marks[idx]==Mark::ChordReady){
                        ov.marks[idx] = Mark::ChordReady;
                    }
                }
            }

            for(int idx=0; idx<N; ++idx){
                uint8_t mask = flagMask[idx];
                if(mask == 0) continue;
                bool pending = (mask & kPendingBit) != 0;
                bool ready = (mask & kReadyBit) != 0;
                if(pending){
                    if(cells[idx] != CellState::Mine){
                        ov.marks[idx] = Mark::FlagForChord;
                    }
                } else if(ready){
                    if(cells[idx] == CellState::Mine){
                        ov.marks[idx] = Mark::FlagForChordReady;
                    }
                }
            }
        }
    }

        }

		// if we still have neither safe moves nor probabilities to guide, fall back to a region-weighted guess
		if(!ov.hasGuaranteedSafe){
			bool anyGuess=false; for(const auto m : ov.marks){ if(m==Mark::Guess){ anyGuess=true; break; } }
			if(!anyGuess){
				std::vector<uint8_t> visited(N, 0);
				struct RegionCandidate {
					std::vector<int> cells;
					int size = 0;
					int parityCount[2] = {0,0};
					bool touchesEdge = false;
					int minCenterDist = std::numeric_limits<int>::max();
				};
				RegionCandidate bestRegion;
				int cx = w/2;
				int cy = h/2;
				auto enqueueNeighbor = [&](int idx, std::vector<int>& stack){ if(!visited[idx]){ visited[idx]=1; stack.push_back(idx); } };
				const int dx4[4] = {1,-1,0,0};
				const int dy4[4] = {0,0,1,-1};
				for(int idx=0; idx<N; ++idx){
					if(visited[idx]) continue;
					if(cells[idx] != CellState::Unknown) continue;
					if(ov.marks[idx] == Mark::Mine) continue;
					visited[idx] = 1;
					RegionCandidate region;
					std::vector<int> stack; stack.push_back(idx);
					while(!stack.empty()){
						int cur = stack.back(); stack.pop_back();
						region.cells.push_back(cur);
						region.size++;
						int x = cur % w;
						int y = cur / w;
						region.parityCount[(x + y) & 1]++;
						if(x==0 || x==w-1 || y==0 || y==h-1) region.touchesEdge = true;
						int dx = x - cx; int dy = y - cy;
						int centerDist = dx*dx + dy*dy;
						if(centerDist < region.minCenterDist) region.minCenterDist = centerDist;
						for(int dir=0; dir<4; ++dir){
							int nx = x + dx4[dir];
							int ny = y + dy4[dir];
							if(!in_bounds(nx, ny, w, h)) continue;
							int nidx = to_index(nx, ny, w);
							if(visited[nidx]) continue;
							if(cells[nidx] != CellState::Unknown) continue;
							if(ov.marks[nidx] == Mark::Mine) continue;
							enqueueNeighbor(nidx, stack);
						}
					}
					int majority = std::max(region.parityCount[0], region.parityCount[1]);
					int bestMajority = std::max(bestRegion.parityCount[0], bestRegion.parityCount[1]);
					if(region.size > bestRegion.size
						|| (region.size == bestRegion.size && majority > bestMajority)
						|| (region.size == bestRegion.size && majority == bestMajority && region.touchesEdge && !bestRegion.touchesEdge)
						|| (region.size == bestRegion.size && majority == bestMajority && region.touchesEdge == bestRegion.touchesEdge && region.minCenterDist < bestRegion.minCenterDist)){
						bestRegion = std::move(region);
					}
				}
				if(!bestRegion.cells.empty()){
					std::vector<uint8_t> inRegion(N, 0);
					for(int cellIdx : bestRegion.cells) inRegion[cellIdx] = 1;
					std::vector<int> regionDist(N, -1);
std::vector<int> regionQueue;
					regionQueue.reserve(bestRegion.cells.size());
					for(int cellIdx : bestRegion.cells){
						int x = cellIdx % w;
						int y = cellIdx / w;
						bool boundary = false;
						for(int dir=0; dir<4; ++dir){
							int nx = x + dx4[dir];
							int ny = y + dy4[dir];
							if(!in_bounds(nx, ny, w, h) || !inRegion[to_index(nx, ny, w)]){ boundary = true; break; }
						}
						if(boundary){
							regionDist[cellIdx] = 0;
							regionQueue.push_back(cellIdx);
						}
					}
					for(size_t qi=0; qi<regionQueue.size(); ++qi){
						int cur = regionQueue[qi];
						int distCur = regionDist[cur];
						int x = cur % w;
						int y = cur / w;
						for(int dir=0; dir<4; ++dir){
							int nx = x + dx4[dir];
							int ny = y + dy4[dir];
							if(!in_bounds(nx, ny, w, h)) continue;
							int nidx = to_index(nx, ny, w);
							if(!inRegion[nidx]) continue;
							if(regionDist[nidx] != -1) continue;
							regionDist[nidx] = distCur + 1;
							regionQueue.push_back(nidx);
						}
					}
					int bestParity = (bestRegion.parityCount[0] >= bestRegion.parityCount[1]) ? 0 : 1;
					bool preferEdge = bestRegion.touchesEdge;
					int bestIdx = -1;
					int bestDepth = -1;
					bool bestOnEdge = false;
					int bestCenterScore = std::numeric_limits<int>::max();
					for(int cellIdx : bestRegion.cells){
						int x = cellIdx % w;
						int y = cellIdx / w;
						int parity = (x + y) & 1;
						if(parity != bestParity) continue;
						int depth = regionDist[cellIdx]; if(depth < 0) depth = 0;
						bool onEdge = (x==0 || x==w-1 || y==0 || y==h-1);
						int dx = x - cx; int dy = y - cy;
						int centerDist = dx*dx + dy*dy;
						bool take = false;
						if(depth > bestDepth + 0){
							take = true;
						} else if(depth == bestDepth){
							if(preferEdge && onEdge && !bestOnEdge) take = true;
							else if((!preferEdge || onEdge == bestOnEdge) && centerDist < bestCenterScore) take = true;
							else if((!preferEdge || onEdge == bestOnEdge) && centerDist == bestCenterScore && cellIdx < bestIdx) take = true;
						}
						if(take){
							bestIdx = cellIdx;
							bestDepth = depth;
							bestOnEdge = onEdge;
							bestCenterScore = centerDist;
						}
					}
					if(bestIdx < 0){
						bestDepth = -1; bestOnEdge = false; bestCenterScore = std::numeric_limits<int>::max();
						for(int cellIdx : bestRegion.cells){
							int x = cellIdx % w;
							int y = cellIdx / w;
							int depth = regionDist[cellIdx]; if(depth < 0) depth = 0;
							bool onEdge = (x==0 || x==w-1 || y==0 || y==h-1);
							int dx = x - cx; int dy = y - cy;
							int centerDist = dx*dx + dy*dy;
							bool take = false;
							if(depth > bestDepth) take = true;
							else if(depth == bestDepth){
								if(preferEdge && onEdge && !bestOnEdge) take = true;
								else if((!preferEdge || onEdge == bestOnEdge) && centerDist < bestCenterScore) take = true;
								else if((!preferEdge || onEdge == bestOnEdge) && centerDist == bestCenterScore && cellIdx < bestIdx) take = true;
							}
							if(take){
								bestIdx = cellIdx;
								bestDepth = depth;
								bestOnEdge = onEdge;
								bestCenterScore = centerDist;
							}
						}
					}
					if(bestIdx >= 0){
						ov.marks[bestIdx] = Mark::Guess;
					}
				} else {
					int bestIdx=-1; int bestDist=std::numeric_limits<int>::max();
					for(int y=0;y<h;++y){ for(int x=0;x<w;++x){ if(board.at(x,y)==CellState::Unknown && ov.marks[to_index(x,y,w)]!=Mark::Mine && ov.marks[to_index(x,y,w)]!=Mark::FlagForChord){ int dx=x-cx, dy=y-cy; int d=dx*dx+dy*dy; if(d<bestDist){ bestDist=d; bestIdx=to_index(x,y,w);} } } }
					if(bestIdx>=0) ov.marks[bestIdx]=Mark::Guess;
				}
			}
		}

    if(cancelled()) return invalid();
    for(int i=0;i<N;++i) {
        if(cells[i]==CellState::Mine || ov.marks[i]==Mark::Mine || ov.marks[i]==Mark::FlagForChord) ov.mineProbability[i]=1;
        else if(cells[i]!=CellState::Unknown || ov.marks[i]==Mark::Safe) ov.mineProbability[i]=0;
    }
	return ov;
}

}
