#include "solver.hpp"

#include <algorithm>
#include <array>
#include <cmath>
#include <functional>
#include <limits>
#include <unordered_map>

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

Overlay compute_overlay(const Board& board, int totalMines, bool enableChords){
	Overlay ov{};
	const int w = board.width();
	const int h = board.height();
	const int N = w * h;
	ov.marks.assign(N, Mark::None);

	if(w==0 || h==0) return ov;

    struct NeighborCache { int w=-1, h=-1; std::vector<std::array<int,8>> nbr; std::vector<int> cnt; };
    static NeighborCache cache;
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
		int y = idxCenter / w; int x = idxCenter % w;
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
            // inconsistent; skip
            continue;
        }
        bool any=false;
        if(remaining == 0 && ucount > 0){
            for(int i=0;i<ucount;++i){ int u=unknownIdx[i]; if(ov.marks[u] != Mark::Safe){ ov.marks[u]=Mark::Safe; any=true; enqueueNbrNumbers(u);} }
        } else if(remaining == ucount && ucount>0){
            for(int i=0;i<ucount;++i){ int u=unknownIdx[i]; if(ov.marks[u] != Mark::Mine){ ov.marks[u]=Mark::Mine; any=true; enqueueNbrNumbers(u);} }
        }
		if(ucount>0){
            for(int k=0;k<neighborCounts[idxCenter];++k){
                int nbNumIdx = neighbors[idxCenter][k];
				int ny = nbNumIdx / w, nx = nbNumIdx % w;
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
                if(rem2<0 || rem2>u2c) continue;
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
                }
            }
        }
        (void)any; // already enqueued impacted numbers when we changed marks
	}

	if(enableChords){
		for(int y=0; y<h; ++y){
			for(int x=0; x<w; ++x){
				int idx = to_index(x,y,w);
				int num = cell_number(board.at(x,y));
                if(num < 0) continue;
				int flagged = 0;
				int unknownCount = 0;
				int neededMineIdx[8]; int neededCnt = 0;
				for(int k=0; k<neighborCounts[idx]; ++k){
					int nb = neighbors[idx][k];
					CellState s = cells[nb];
					if(s == CellState::Mine){ flagged++; continue; }
					if(s == CellState::Unknown){
						unknownCount++;
						if(ov.marks[nb] == Mark::Mine){ neededMineIdx[neededCnt++] = nb; }
					}
				}
				int missing = num - flagged;
				if(missing < 0 || missing > unknownCount) continue;
                if(missing == 0){
                    int safeClicks = unknownCount;
                    if(safeClicks > 1){
                        if(ov.marks[idx] == Mark::None){ ov.marks[idx] = Mark::Chord; }
                    }
                } else {
					if(neededCnt == missing){
                        int chordCost = missing + 1;
                        int safeClicks = unknownCount - missing;
                        if(chordCost < safeClicks){
							if(ov.marks[idx] == Mark::None){ ov.marks[idx] = Mark::Chord; }
                            for(int t=0; t<neededCnt; ++t){ int nb = neededMineIdx[t]; if(cells[nb] != CellState::Mine){ ov.marks[nb] = Mark::FlagForChord; } }
						}
					}
				}
			}
		}
	}

	// determine if we have any guaranteed safe so far
	for(const auto m : ov.marks){ if(m == Mark::Safe){ ov.hasGuaranteedSafe = true; break; } }

	// Build frontier components and compute probabilities using exact enumeration where feasible,
	// then fall back to BP. Use results to mark certain Safe/Mine and best Guess if needed.
        std::vector<uint8_t> isFrontier(N,0);
        for(int i=0;i<N;++i){ if(board.data()[i]==CellState::Unknown){ for(int k=0;k<neighborCounts[i];++k){ int nb=neighbors[i][k]; if(cell_number(board.data()[nb])>=0){ isFrontier[i]=1; break; } } } }

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
        std::unordered_map<int,int> rootToCid; rootToCid.reserve(N/4+1);
        int compCount=0;
        for(int i=0;i<N;++i){ if(isFrontier[i]){ int r=findp(i); auto it=rootToCid.find(r); if(it==rootToCid.end()){ rootToCid[r]=compCount; compId[i]=compCount; compCount++; } else { compId[i]=it->second; } } }

        ov.mineProbability.assign(N, -1.0);

        struct Con { int num; int knownMines; std::vector<int> uidx; };
        struct CompData {
            std::vector<int> U;
            std::vector<Con> cons;
            int m = 0;
            std::vector<long double> waysK;
            std::vector<std::vector<long double>> mineWaysPerCellK;
            long double totalSolutions = 0.0L;
            bool enumerated = false;
        };
        std::vector<CompData> comps(compCount);
        for(int i=0;i<N;++i){ if(compId[i]>=0){ comps[compId[i]].U.push_back(i); } }
        for(int idx=0; idx<N; ++idx){ if(cell_number(cells[idx])>=0){
            int uc=0; int tmp[8]; int known=0;
            for(int k=0;k<neighborCounts[idx];++k){ int nb=neighbors[idx][k]; if(cells[nb]==CellState::Mine || ov.marks[nb]==Mark::Mine){ known++; continue; } if(isFrontier[nb]) tmp[uc++]=nb; }
            if(uc==0) continue;
            int cid = compId[tmp[0]];
            Con con; con.num = cell_number(cells[idx]); con.knownMines = known;
            for(int t=0;t<uc;++t){ int g=tmp[t]; for(int j=0;j<(int)comps[cid].U.size(); ++j){ if(comps[cid].U[j]==g){ con.uidx.push_back(j); break; } } }
            comps[cid].cons.push_back(std::move(con));
        } }

        const int ENUM_CAP = 22;
        for(int cid=0; cid<compCount; ++cid){
            auto& C = comps[cid];
            C.m = (int)C.U.size();
            if(C.m==0 || (int)C.cons.size()==0){ C.enumerated=false; continue; }
            if(C.m > ENUM_CAP){ C.enumerated=false; continue; }
            C.waysK.assign(C.m+1, 0.0L);
            C.mineWaysPerCellK.assign(C.m, std::vector<long double>(C.m+1, 0.0L));
            std::vector<uint8_t> assign(C.m, 255);
            std::function<void(int,int)> dfs = [&](int idx, int minesSoFar){
                // local feasibility against all constraints
                for(const auto& con : C.cons){
                    int need = con.num - con.knownMines;
                    int placed=0, unk=0;
                    for(int ui: con.uidx){ if(assign[ui]==1) placed++; else if(assign[ui]==255) unk++; }
                    if(need < placed) return; if(need > placed + unk) return;
                }
                if(idx==C.m){
                    C.totalSolutions += 1.0L;
                    C.waysK[minesSoFar] += 1.0L;
                    for(int t=0;t<C.m;++t){ if(assign[t]==1) C.mineWaysPerCellK[t][minesSoFar] += 1.0L; }
                    return;
                }
                assign[idx]=0; dfs(idx+1, minesSoFar);
                assign[idx]=1; dfs(idx+1, minesSoFar+1);
                assign[idx]=255;
            };
            dfs(0,0);
            C.enumerated = (C.totalSolutions > 0.0L);
            if(C.enumerated && totalMines==0){
                for(int t=0;t<C.m;++t){ long double mw=0.0L; for(int k=0;k<=C.m;++k) mw += C.mineWaysPerCellK[t][k]; ov.mineProbability[C.U[t]] = (double)(mw / C.totalSolutions); }
            }
        }

        const int BP_MAX_ITERS = 20;
        const double BP_DAMP = 0.5;
        const double BP_EPS = 1e-6;
        for(int cid=0; cid<compCount; ++cid){
            auto& C = comps[cid];
            if(C.m==0 || (int)C.cons.size()==0) continue;
            if(C.enumerated) continue; // already handled exactly

            // build variable-to-constraint adjacency
            struct EdgeRef { int conIndex; int posInCon; };
            std::vector<std::vector<EdgeRef>> varEdges(C.m);
            for(int j=0;j<(int)C.cons.size();++j){
                const auto& con = C.cons[j];
                for(int p=0;p<(int)con.uidx.size();++p){ int v=con.uidx[p]; varEdges[v].push_back({j,p}); }
            }
            auto findEdgeIndex = [&](int v, int conIndex)->int{
                const auto& edges = varEdges[v];
                for(int i=0;i<(int)edges.size();++i){ if(edges[i].conIndex==conIndex) return i; }
                return -1;
            };

            std::vector<std::vector<double>> vToF(C.m);
            for(int v=0; v<C.m; ++v){
                vToF[v].assign(varEdges[v].size(), 0.5);
                double prior = 0.5;
                int g = C.U[v];
                if(g>=0 && g<(int)ov.mineProbability.size() && ov.mineProbability[g] >= 0.0) prior = std::min(1.0-1e-6, std::max(1e-6, ov.mineProbability[g]));
                for(double& x : vToF[v]) x = prior;
            }
            std::vector<std::vector<double>> fToV(C.cons.size());
            for(size_t j=0;j<C.cons.size();++j){ fToV[j].assign(C.cons[j].uidx.size(), 0.5); }

            for(int it=0; it<BP_MAX_ITERS; ++it){
                for(int j=0;j<(int)C.cons.size(); ++j){
                    const auto& con = C.cons[j];
                    const int s = (int)con.uidx.size();
                    const int need = con.num - con.knownMines;
                    for(int p=0; p<s; ++p){
                        std::vector<long double> dp(s, 0.0L);
                        dp[0] = 1.0L;
                        for(int q=0; q<s; ++q){ if(q==p) continue; int vv = con.uidx[q]; int eidx = findEdgeIndex(vv, j); double prob = 0.5; if(eidx>=0) prob = vToF[vv][eidx];
                            std::vector<long double> ndp(s, 0.0L);
                            for(int k=0;k<s-1;++k){ if(dp[k]==0.0L) continue; ndp[k] += dp[k] * (1.0L - (long double)prob); ndp[k+1] += dp[k] * (long double)prob; }
                            dp.swap(ndp);
                        }
                        long double A = 0.0L, B = 0.0L;
                        if(need-1 >= 0 && need-1 <= s-1) A = dp[need-1];
                        if(need >= 0 && need <= s-1) B = dp[need];
                        double msg = 0.5;
                        if(A==0.0L && B==0.0L){ msg = 0.5; }
                        else if(A==0.0L){ msg = 0.0; }
                        else if(B==0.0L){ msg = 1.0; }
                        else { double a=(double)A, b=(double)B; msg = a / (a + b); }
                        msg = BP_DAMP * fToV[j][p] + (1.0 - BP_DAMP) * std::min(1.0-1e-9, std::max(1e-9, msg));
                        fToV[j][p] = msg;
                    }
                }

                for(int v=0; v<C.m; ++v){
                    int deg = (int)varEdges[v].size();
                    if(deg==0) continue;
                    std::vector<double> pref1(deg+1, 1.0), pref0(deg+1, 1.0);
                    for(int i=0;i<deg;++i){ const auto& e = varEdges[v][i]; double m1 = fToV[e.conIndex][e.posInCon]; pref1[i+1] = pref1[i]*m1; pref0[i+1] = pref0[i]*(1.0 - m1); }
                    std::vector<double> suf1(deg+1, 1.0), suf0(deg+1, 1.0);
                    for(int i=deg-1;i>=0;--i){ const auto& e = varEdges[v][i]; double m1 = fToV[e.conIndex][e.posInCon]; suf1[i] = suf1[i+1]*m1; suf0[i] = suf0[i+1]*(1.0 - m1); }
                    for(int i=0;i<deg;++i){
                        double p1 = pref1[i]*suf1[i+1];
                        double p0 = pref0[i]*suf0[i+1];
                        double msg = 0.5;
                        if(p1==0.0 && p0==0.0) msg = 0.5; else msg = p1 / (p1 + p0);
                        msg = std::min(1.0-BP_EPS, std::max(BP_EPS, msg));
                        vToF[v][i] = BP_DAMP * vToF[v][i] + (1.0 - BP_DAMP) * msg;
                    }
                }
            }

            // compute final beliefs per variable
            for(int v=0; v<C.m; ++v){
                int deg = (int)varEdges[v].size();
                double prod1=1.0, prod0=1.0;
                for(int i=0;i<deg;++i){ const auto& e = varEdges[v][i]; double m1 = fToV[e.conIndex][e.posInCon]; prod1 *= m1; prod0 *= (1.0 - m1); }
                double p = 0.5; if(!(prod1==0.0 && prod0==0.0)) p = prod1 / (prod1 + prod0);
                p = std::min(1.0-BP_EPS, std::max(BP_EPS, p));
                int g = C.U[v]; if(g>=0 && g<(int)ov.mineProbability.size()) ov.mineProbability[g] = p;
            }
        }

        bool anyEnumerated=false; for(const auto& C : comps){ if(C.enumerated){ anyEnumerated=true; break; } }
        if(totalMines>0 && anyEnumerated){
            int knownGlobal=0, unknownGlobal=0; for(int i=0;i<N;++i){ if(cells[i]==CellState::Mine || ov.marks[i]==Mark::Mine) knownGlobal++; else if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Safe) unknownGlobal++; }
            int unknownInEnum=0; for(const auto& C : comps){ if(C.enumerated) unknownInEnum += C.m; }
            int freeUnknown = std::max(0, unknownGlobal - unknownInEnum);
            int remainingGlobal = std::max(0, totalMines - knownGlobal);
            // comb(freeUnknown, j)
            std::vector<long double> combFree(freeUnknown+1, 0.0L); if(freeUnknown>=0){ combFree[0]=1.0L; for(int j=1;j<=freeUnknown;++j){ combFree[j] = combFree[j-1] * (long double)(freeUnknown - (j-1)) / (long double)j; } }
            // gather enumerated indices
            std::vector<int> enumIdx; for(int i=0;i<(int)comps.size();++i){ if(comps[i].enumerated) enumIdx.push_back(i); }
            int E=(int)enumIdx.size();
            std::vector<std::vector<long double>> prefix(E+1), suffix(E+1);
            prefix[0] = std::vector<long double>(1, 1.0L);
            for(int i=0;i<E;++i){ const auto& C = comps[enumIdx[i]]; std::vector<long double> next(prefix[i].size()+C.m, 0.0L); for(size_t a=0;a<prefix[i].size();++a){ for(int k=0;k<=C.m;++k){ next[a+k] += prefix[i][a]*C.waysK[k]; } } prefix[i+1].swap(next); }
            suffix[E] = std::vector<long double>(1, 1.0L);
            for(int i=E-1;i>=0;--i){ const auto& C = comps[enumIdx[i]]; std::vector<long double> next(suffix[i+1].size()+C.m, 0.0L); for(size_t b=0;b<suffix[i+1].size();++b){ for(int k=0;k<=C.m;++k){ next[b+k] += suffix[i+1][b]*C.waysK[k]; } } suffix[i].swap(next); }
            long double totalWeight=0.0L; for(size_t r=0;r<prefix[E].size();++r){ int needOut = remainingGlobal - (int)r; if(needOut>=0 && needOut<=freeUnknown) totalWeight += prefix[E][r]*combFree[needOut]; }
            if(totalWeight>0.0L){
                for(int ii=0; ii<E; ++ii){ int cid = enumIdx[ii]; auto& C = comps[cid];
                    std::vector<long double> other(prefix[ii].size()+suffix[ii+1].size()-1, 0.0L); for(size_t a=0;a<prefix[ii].size();++a){ for(size_t b=0;b<suffix[ii+1].size();++b){ other[a+b] += prefix[ii][a]*suffix[ii+1][b]; } }
                    for(int t=0;t<C.m;++t){ long double numer=0.0L; for(int k=0;k<=C.m;++k){ if(C.mineWaysPerCellK[t][k]==0.0L) continue; for(size_t rOther=0;rOther<other.size();++rOther){ int needOut = remainingGlobal - ((int)rOther + k); if(needOut<0 || needOut>freeUnknown) continue; numer += C.mineWaysPerCellK[t][k] * other[rOther] * combFree[needOut]; } } ov.mineProbability[C.U[t]] = (double)(numer / totalWeight); }
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
        if(totalMines>0){
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
                    double r = countRisk[i]>0 ? (sumRisk[i]/(double)countRisk[i]) : 1.0;
                    ov.mineProbability[i] = r;
                }
            }
        }

        bool hasProb=false; double bestProb=std::numeric_limits<double>::infinity();
        if(ov.mineProbability.size()==(size_t)N){
            for(int i=0;i<N;++i){ if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Mine && ov.mineProbability[i]>=0.0){ hasProb=true; if(ov.mineProbability[i] < bestProb) bestProb = ov.mineProbability[i]; } }
        }
        std::vector<int> changedByProb; changedByProb.reserve(N/8+1);
        if(hasProb && std::isfinite(bestProb)){
            const double eps=1e-9;
            for(int i=0;i<N;++i){
                if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Mine && ov.mineProbability[i]>=0.0){
                    if(ov.mineProbability[i] <= eps){ if(ov.marks[i]!=Mark::Safe){ ov.marks[i]=Mark::Safe; ov.hasGuaranteedSafe = true; changedByProb.push_back(i); } }
                    else if(ov.mineProbability[i] >= 1.0 - eps){ if(ov.marks[i]!=Mark::Mine){ ov.marks[i]=Mark::Mine; changedByProb.push_back(i); } }
                }
            }
            // Only set Guess when no guaranteed safe exists
            if(!ov.hasGuaranteedSafe){
                for(int i=0;i<N;++i){
                    if(cells[i]==CellState::Unknown && ov.marks[i]!=Mark::Mine && ov.mineProbability[i]>=0.0){
                        if(std::abs(ov.mineProbability[i] - bestProb) <= eps){ ov.marks[i]=Mark::Guess; }
                    }
                }
            }

    // chord valuation using probabilities and a simple cascade heuristic
    if(enableChords && ov.mineProbability.size()==(size_t)N){
        // clear prior chord/flag marks while preserving Safe/Mine/Guess
        for(int i=0;i<N;++i){ if(ov.marks[i]==Mark::Chord || ov.marks[i]==Mark::FlagForChord) ov.marks[i]=Mark::None; }

        auto probOf = [&](int idx)->double{
            double p = 1.0; // default to risky if unknown
            if(idx>=0 && idx<N && cells[idx]==CellState::Unknown){
                p = ov.mineProbability[idx];
                if(!(p>=0.0 && p<=1.0)) p = 1.0; // clamp invalids to risky
            } else if(idx>=0 && idx<N && (cells[idx]==CellState::Mine || ov.marks[idx]==Mark::Mine)){
                p = 1.0;
            } else {
                p = 0.0;
            }
            if(p<0.0) p=0.0; if(p>1.0) p=1.0;
            return p;
        };

        auto willFlagForChord = [&](int idx)->bool{ return ov.marks[idx]==Mark::Mine && cells[idx]!=CellState::Mine; };
        const double kProbEps = 1e-9;

        const double kCascadeDampen = 0.5; // dampen cascade bonus to avoid double counting

        struct ChordCandidate {
            int centerIdx;
            int missing;
            std::vector<int> unknowns;
            std::vector<int> flagList; // flags to place for this chord
            int clicks; // flags to place (excluding already-flagged board cells) + 1 chord click
            double expectedValue; // reveals + cascade
            double bestSingle; // best single-click alternative near this center
            bool inProgress; // true if user just placed a flag enabling a break-even chord
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

        auto computeBestSingle = [&](const std::vector<int>& unknowns, const std::vector<int>& plannedFlags, int centerA, int centerB, double dampen)->double{
            double bestSingle = 0.0;
            for(int u : unknowns){
                bool isPlannedFlag = false; for(int f : plannedFlags){ if(f==u){ isPlannedFlag=true; break; } }
                if(isPlannedFlag) continue;
                double pu = probOf(u);
                double singleVal = (1.0 - pu);
                int uc2=0; int list2[8];
                for(int kk=0; kk<neighborCounts[u]; ++kk){ int m = neighbors[u][kk]; if(m==centerA || m==centerB) continue; if(cells[m]==CellState::Unknown){ bool pf=false; for(int f : plannedFlags){ if(f==m){ pf=true; break; } } if(pf) continue; list2[uc2++] = m; } }
                if(uc2>0){ long double zeroProb = 1.0L; for(int t=0;t<uc2;++t){ double pm = probOf(list2[t]); zeroProb *= (long double)(1.0 - pm); } singleVal += (double)zeroProb * (double)uc2 * dampen; }
                if(singleVal > bestSingle) bestSingle = singleVal;
            }
            return bestSingle;
        };

        std::vector<ChordCandidate> candidates; candidates.reserve(w*h/2);
        for(int y=0; y<h; ++y){
            for(int x=0; x<w; ++x){
                int idx = to_index(x,y,w);
                int num = cell_number(cells[idx]); if(num<0) continue;

                int flagged=0; std::vector<int> unknowns; unknowns.reserve(8);
                std::vector<int> needFlags; needFlags.reserve(8);
                int userFlags = 0; // flags present on board but not deduced by solver
                for(int k=0;k<neighborCounts[idx];++k){
                    int nb = neighbors[idx][k];
                    CellState s = cells[nb];
                    if(s==CellState::Mine){
                        flagged++;
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
                    for(int u : unknowns){
                        if(ov.marks[u]==Mark::Mine) continue;
                        double pu = probOf(u);
                        if(pu >= 1.0 - kProbEps) certainFlags.push_back(u);
                    }
                    if((int)certainFlags.size() != missing) continue; // cannot chord safely
                } else {
                    // If user just placed at least one neighbor flag, allow a break-even chord with 1 unknown
                    // to persist planned action; otherwise require 2+ unknowns to ensure savings.
                    if(!((int)unknowns.size() >= 1 && userFlags > 0)){
                        if((int)unknowns.size() < 2) continue; // require savings when not in-progress
                    }
                }

                int numFlagsToPlace = 0; for(int f : certainFlags){ if(cells[f] != CellState::Mine) numFlagsToPlace++; }
                int clicks = (missing==0) ? 1 : (numFlagsToPlace + 1);

                auto ev = computeExpected(unknowns, certainFlags, idx, -1, kCascadeDampen);
                double totalValue = ev.first + ev.second;
                double bestSingle = computeBestSingle(unknowns, certainFlags, idx, -1, kCascadeDampen);

                bool inProgress = (missing==0 && (int)unknowns.size()==1 && userFlags>0);
                candidates.push_back(ChordCandidate{ idx, missing, std::move(unknowns), std::move(certainFlags), clicks, totalValue, bestSingle, inProgress });
            }
        }

        // find synergistic pairs that can share flags; prefer pairs with strong positive net improvement
        struct PairPick { int a; int b; int clicks; double value; double bestTwoSingles; double improvement; std::vector<int> unionFlags; std::vector<int> unionUnknowns; };
        std::vector<PairPick> pairPicks; pairPicks.reserve(candidates.size());

        auto unionUnique = [&](const std::vector<int>& A, const std::vector<int>& B){
            std::vector<int> out = A; out.reserve(A.size()+B.size());
            for(int v : B){ bool found=false; for(int u : out){ if(u==v){ found=true; break; } } if(!found) out.push_back(v); }
            return out;
        };

        const double kMargin = 1e-3;

        for(int i=0;i<(int)candidates.size();++i){
            for(int j=i+1;j<(int)candidates.size();++j){
                const auto& A = candidates[i]; const auto& B = candidates[j];
                // spatial locality: only consider if centers are neighbors or share unknown/flags
                int ax=A.centerIdx%w, ay=A.centerIdx/w; int bx=B.centerIdx%w, by=B.centerIdx/w;
                int dx = std::abs(ax-bx), dy = std::abs(ay-by);
                bool close = (dx<=2 && dy<=2);
                if(!close){
                    bool overlap=false;
                    for(int u1 : A.unknowns){ for(int u2 : B.unknowns){ if(u1==u2){ overlap=true; break; } } if(overlap) break; }
                    if(!overlap){
                        for(int f1 : A.flagList){ for(int f2 : B.flagList){ if(f1==f2){ overlap=true; break; } } if(overlap) break; }
                    }
                    if(!overlap) continue;
                }

                // to share flags effectively, require that flags do not conflict, and both are safe chordable
                std::vector<int> flagUnion = unionUnique(A.flagList, B.flagList);
                // clicks for pair: unique flags to place + 2 chord clicks
                int flagsToPlace=0; for(int f : flagUnion){ if(cells[f] != CellState::Mine) flagsToPlace++; }
                int pairClicks = flagsToPlace + 2;

                // unknown union used for valuation (avoid double counting)
                std::vector<int> unknownUnion = unionUnique(A.unknowns, B.unknowns);

                auto evPair = computeExpected(unknownUnion, flagUnion, A.centerIdx, B.centerIdx, kCascadeDampen);
                double pairValue = evPair.first + evPair.second;

                // approximate best two singles alternative from union unknowns
                // compute top two single values greedily without replacement
                struct ValIdx { double v; int idx; };
                std::vector<ValIdx> vals; vals.reserve(unknownUnion.size());
                for(int u : unknownUnion){ bool isFlag=false; for(int f:flagUnion){ if(f==u){ isFlag=true; break; } } if(isFlag) continue; double pu = probOf(u); double v = (1.0 - pu);
                    int uc2=0; int list2[8]; for(int kk=0; kk<neighborCounts[u]; ++kk){ int m = neighbors[u][kk]; if(m==A.centerIdx || m==B.centerIdx) continue; if(cells[m]==CellState::Unknown){ bool pf=false; for(int f:flagUnion){ if(f==m){ pf=true; break; } } if(pf) continue; list2[uc2++] = m; } }
                    if(uc2>0){ long double zeroProb = 1.0L; for(int t=0;t<uc2;++t){ double pm = probOf(list2[t]); zeroProb *= (long double)(1.0 - pm); } v += (double)zeroProb * (double)uc2 * kCascadeDampen; }
                    vals.push_back({v, u}); }
                std::sort(vals.begin(), vals.end(), [](const ValIdx& a, const ValIdx& b){ return a.v>b.v; });
                double bestTwoSingles = 0.0; if(!vals.empty()) bestTwoSingles += vals[0].v; if(vals.size()>1) bestTwoSingles += vals[1].v;

                double improvement = pairValue - std::max((double)pairClicks, bestTwoSingles);
                if(improvement > kMargin){
                    pairPicks.push_back(PairPick{ i, j, pairClicks, pairValue, bestTwoSingles, improvement, std::move(flagUnion), std::move(unknownUnion) });
                }
            }
        }

        // greedy select non-overlapping best pairs
        std::vector<uint8_t> picked(candidates.size(), 0);
        std::sort(pairPicks.begin(), pairPicks.end(), [](const PairPick& a, const PairPick& b){ return a.improvement > b.improvement; });
        for(const auto& p : pairPicks){
            if(picked[p.a] || picked[p.b]) continue;
            // mark both chords and required flags
            const auto& A = candidates[p.a]; const auto& B = candidates[p.b];
            if(ov.marks[A.centerIdx]==Mark::None) ov.marks[A.centerIdx]=Mark::Chord;
            if(ov.marks[B.centerIdx]==Mark::None) ov.marks[B.centerIdx]=Mark::Chord;
            for(int f : p.unionFlags){ if(cells[f] != CellState::Mine){ ov.marks[f]=Mark::FlagForChord; } }
            picked[p.a]=picked[p.b]=1;
        }

        // Always keep in-progress single chords visible (user recently placed enabling flag)
        for(int i=0;i<(int)candidates.size();++i){ if(picked[i]) continue; const auto& C = candidates[i];
            if(C.inProgress){ if(ov.marks[C.centerIdx]==Mark::None) ov.marks[C.centerIdx]=Mark::Chord; picked[i]=1; }
        }

        // pick remaining profitable singles
        for(int i=0;i<(int)candidates.size();++i){ if(picked[i]) continue; const auto& C = candidates[i];
            double improvement = C.expectedValue - std::max((double)C.clicks, C.bestSingle);
            if(improvement > kMargin){
                if(ov.marks[C.centerIdx]==Mark::None) ov.marks[C.centerIdx]=Mark::Chord;
                for(int f : C.flagList){ if(cells[f] != CellState::Mine){ ov.marks[f]=Mark::FlagForChord; } }
                picked[i]=1;
            }
        }
    }

        // If new certain marks were added via probabilities, re-run local propagation around them
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
                    // Only count preexisting flags or deduced mines here; do not treat FlagForChord as flagged
                    if(s2 == CellState::Mine || ov.marks[nb2] == Mark::Mine){ knownMines2++; continue; }
                    if(s2 == CellState::Unknown && ov.marks[nb2] != Mark::Safe){ unknownIdx2[ucount2++] = nb2; }
                }
                int remaining2 = num2 - knownMines2;
                if(remaining2 < 0 || remaining2 > ucount2){ continue; }
                if(remaining2 == 0 && ucount2 > 0){
                    for(int i=0;i<ucount2;++i){ int u=unknownIdx2[i]; if(ov.marks[u] != Mark::Safe){ ov.marks[u]=Mark::Safe; enqueueNbrNumbers(u);} }
                } else if(remaining2 == ucount2 && ucount2>0){
                    for(int i=0;i<ucount2;++i){ int u=unknownIdx2[i]; if(ov.marks[u] == Mark::FlagForChord){ continue; } if(ov.marks[u] != Mark::Mine){ ov.marks[u] = Mark::Mine; enqueueNbrNumbers(u);} }
                }
            }
            // refresh guaranteed safe after propagation
            ov.hasGuaranteedSafe = false; for(const auto m : ov.marks){ if(m==Mark::Safe){ ov.hasGuaranteedSafe=true; break; } }
        }
        }

        // If we still have neither safe moves nor probabilities to guide, fall back to a center guess
        if(!ov.hasGuaranteedSafe){
            bool anyGuess=false; for(const auto m : ov.marks){ if(m==Mark::Guess){ anyGuess=true; break; } }
            if(!anyGuess){
                int cx=w/2, cy=h/2; int bestIdx=-1; int bestDist=std::numeric_limits<int>::max();
                for(int y=0;y<h;++y){ for(int x=0;x<w;++x){ if(board.at(x,y)==CellState::Unknown){ int dx=x-cx, dy=y-cy; int d=dx*dx+dy*dy; if(d<bestDist){ bestDist=d; bestIdx=to_index(x,y,w);} } } }
                if(bestIdx>=0) ov.marks[bestIdx]=Mark::Guess;
            }
        }

	return ov;
}

}



