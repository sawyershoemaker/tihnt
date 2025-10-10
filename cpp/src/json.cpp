#include <string>
#include <vector>
#include <cctype>

#include "proto.hpp"

namespace proto {

static void skip_ws(const std::string& s, size_t& i){ while(i<s.size() && std::isspace((unsigned char)s[i])) ++i; }
static bool match(const std::string& s, size_t& i, char c){ skip_ws(s,i); if(i<s.size() && s[i]==c){ ++i; return true; } return false; }
static bool parse_string(const std::string& s, size_t& i, std::string& out){ skip_ws(s,i); if(i>=s.size()||s[i]!='"') return false; ++i; out.clear(); while(i<s.size()&&s[i]!='"'){ if(s[i]=='\\' && i+1<s.size()){ ++i; out.push_back(s[i++]); } else { out.push_back(s[i++]); } } if(i>=s.size()) return false; ++i; return true; }
static bool parse_int(const std::string& s, size_t& i, int& out){ skip_ws(s,i); bool neg=false; if(i<s.size()&&(s[i]=='-'||s[i]=='+')){ neg=s[i]=='-'; ++i; } int v=0; bool any=false; while(i<s.size() && std::isdigit((unsigned char)s[i])){ any=true; v = v*10 + (s[i++]-'0'); } if(!any) return false; out = neg?-v:v; return true; }
static bool parse_double(const std::string& s, size_t& i, double& out){
    skip_ws(s,i);
    size_t start=i;
    if(i<s.size() && (s[i]=='-'||s[i]=='+')) ++i;
    bool any=false; while(i<s.size() && std::isdigit((unsigned char)s[i])){ any=true; ++i; }
    if(i<s.size() && s[i]=='.'){ ++i; while(i<s.size() && std::isdigit((unsigned char)s[i])){ any=true; ++i; } }
    if(!any) return false;
    try { out = std::stod(s.substr(start, i-start)); } catch(...) { return false; }
    return true;
}

static bool parse_object_start(const std::string& s, size_t& i){ return match(s,i,'{'); }
static bool parse_object_end(const std::string& s, size_t& i){ return match(s,i,'}'); }
static bool parse_array_start(const std::string& s, size_t& i){ return match(s,i,'['); }
static bool parse_array_end(const std::string& s, size_t& i){ return match(s,i,']'); }
static bool parse_colon(const std::string& s, size_t& i){ return match(s,i,':'); }
static bool parse_comma(const std::string& s, size_t& i){ return match(s,i,','); }

using proto::FullMsg; using proto::DeltaMsg; using proto::ParsedMessage; using proto::MsgType; using proto::BindMsg;

static bool parse_cells_array(const std::string& s, size_t& i, std::vector<game::CellState>& out){
    if(!parse_array_start(s,i)) return false; out.clear();
    size_t j=i; skip_ws(s,j); if(j<s.size() && s[j]==']'){ i=j+1; return true; }
    for(;;){ int v; if(!parse_int(s,i,v)) return false; out.push_back(static_cast<game::CellState>(v)); size_t k=i; if(parse_comma(s,k)){ i=k; continue; } if(parse_array_end(s,i)) break; else return false; }
    return true;
}

static bool parse_updates_array(const std::string& s, size_t& i, std::vector<game::CellUpdate>& out){
    if(!parse_array_start(s,i)) return false; out.clear();
    size_t j=i; skip_ws(s,j); if(j<s.size() && s[j]==']'){ i=j+1; return true; }
    for(;;){ if(!parse_object_start(s,i)) return false; int x=0,y=0,sv=0; bool gotx=false, goty=false, gots=false; for(;;){ std::string key; if(!parse_string(s,i,key)) return false; if(!parse_colon(s,i)) return false; if(key=="x"){ if(!parse_int(s,i,x)) return false; gotx=true; } else if(key=="y"){ if(!parse_int(s,i,y)) return false; goty=true; } else if(key=="s"){ if(!parse_int(s,i,sv)) return false; gots=true; } else { return false; } size_t k=i; if(parse_comma(s,k)){ i=k; continue; } if(parse_object_end(s,i)) break; else return false; } if(!(gotx&&goty&&gots)) return false; out.push_back(game::CellUpdate{x,y,static_cast<game::CellState>(sv)}); size_t k=i; if(parse_comma(s,k)){ i=k; continue; } if(parse_array_end(s,i)) break; else return false; }
    return true;
}

static bool parse_full(const std::string& s, size_t& i, FullMsg& out){
    if(!parse_object_start(s,i)) return false; bool gotw=false,goth=false,gotcells=false; for(;;){ std::string key; if(!parse_string(s,i,key)) return false; if(!parse_colon(s,i)) return false; if(key=="w"){ if(!parse_int(s,i,out.w)) return false; gotw=true; } else if(key=="h"){ if(!parse_int(s,i,out.h)) return false; goth=true; } else if(key=="cells"){ if(!parse_cells_array(s,i,out.cells)) return false; gotcells=true; } else if(key=="cell_px"){ parse_int(s,i,out.cell_px); } else if(key=="ox"){ parse_int(s,i,out.origin_x); } else if(key=="oy"){ parse_int(s,i,out.origin_y); } else if(key=="mines_total"){ parse_int(s,i,out.mines_total); } else if(key=="rect_l"){ parse_double(s,i,out.rect_l); } else if(key=="rect_t"){ parse_double(s,i,out.rect_t); } else if(key=="rect_w"){ parse_double(s,i,out.rect_w); } else if(key=="rect_h"){ parse_double(s,i,out.rect_h); } else if(key=="vv_x"){ parse_double(s,i,out.vv_x); } else if(key=="vv_y"){ parse_double(s,i,out.vv_y); } else if(key=="vv_scale"){ parse_double(s,i,out.vv_scale); } else if(key=="dpr"){ parse_double(s,i,out.dpr); } else {
            int dummy; if(!parse_int(s,i,dummy)) { std::string tmp; if(!parse_string(s,i,tmp)) { double d; if(!parse_double(s,i,d)) return false; } }
        } size_t k=i; if(parse_comma(s,k)){ i=k; continue; } if(parse_object_end(s,i)) break; else return false; }
    return gotw&&goth&&gotcells;
}

static bool parse_delta(const std::string& s, size_t& i, DeltaMsg& out){
    if(!parse_object_start(s,i)) return false; bool gotu=false; for(;;){ std::string key; if(!parse_string(s,i,key)) return false; if(!parse_colon(s,i)) return false; if(key=="updates"){ if(!parse_updates_array(s,i,out.updates)) return false; gotu=true; } else if(key=="cell_px"){ parse_int(s,i,out.cell_px); } else if(key=="ox"){ parse_int(s,i,out.origin_x); } else if(key=="oy"){ parse_int(s,i,out.origin_y); } else if(key=="mines_total"){ parse_int(s,i,out.mines_total); } else if(key=="rect_l"){ parse_double(s,i,out.rect_l); } else if(key=="rect_t"){ parse_double(s,i,out.rect_t); } else if(key=="rect_w"){ parse_double(s,i,out.rect_w); } else if(key=="rect_h"){ parse_double(s,i,out.rect_h); } else if(key=="vv_x"){ parse_double(s,i,out.vv_x); } else if(key=="vv_y"){ parse_double(s,i,out.vv_y); } else if(key=="vv_scale"){ parse_double(s,i,out.vv_scale); } else if(key=="dpr"){ parse_double(s,i,out.dpr); } else { int dummy; if(!parse_int(s,i,dummy)) { std::string tmp; if(!parse_string(s,i,tmp)) { double d; if(!parse_double(s,i,d)) return false; } } } size_t k=i; if(parse_comma(s,k)){ i=k; continue; } if(parse_object_end(s,i)) break; else return false; }
    return gotu;
}

static bool parse_bind(const std::string& s, size_t& i, BindMsg& out){
    if(!parse_object_start(s,i)) return false; bool gotp=false;
    for(;;){ std::string key; if(!parse_string(s,i,key)) return false; if(!parse_colon(s,i)) return false; if(key=="pid"){ if(!parse_int(s,i,out.pid)) return false; gotp=true; }
        else { int dummy; if(!parse_int(s,i,dummy)) { std::string tmp; if(!parse_string(s,i,tmp)) { double d; if(!parse_double(s,i,d)) return false; } } }
        size_t k=i; if(parse_comma(s,k)){ i=k; continue; } if(parse_object_end(s,i)) break; else return false; }
    return gotp;
}

static bool parse_root(const std::string& s, ParsedMessage& out){
    if(s.find("\"type\":\"full\"") != std::string::npos){
        size_t r=0; bool ok = parse_full(s,r,out.full); if(ok){ out.type=MsgType::Full; } return ok;
    }
    if(s.find("\"type\":\"delta\"") != std::string::npos){
        size_t r=0; bool ok = parse_delta(s,r,out.delta); if(ok){ out.type=MsgType::Delta; } return ok;
    }
    if(s.find("\"type\":\"bind\"") != std::string::npos){
        size_t r=0; bool ok = parse_bind(s,r,out.bind); if(ok){ out.type=MsgType::Bind; } return ok;
    }
    size_t i=0; if(!parse_object_start(s,i)) return false; std::string type; bool gott=false;
    size_t j=i; for(;;){ std::string key; if(!parse_string(s,j,key)) return false; if(!parse_colon(s,j)) return false; if(key=="type"){ if(!parse_string(s,j,type)) return false; gott=true; } else {
            int dummy; if(!parse_int(s,j,dummy)) { std::string tmp; if(!parse_string(s,j,tmp)) return false; }
        }
        size_t k=j; if(parse_comma(s,k)){ j=k; continue; } if(parse_object_end(s,j)) break; else return false; }
    if(!gott) return false;
    size_t r=0; if(type=="full"){ bool ok = parse_full(s,r,out.full); if(ok){ out.type=MsgType::Full; } return ok; }
    if(type=="delta"){ bool ok = parse_delta(s,r,out.delta); if(ok){ out.type=MsgType::Delta; } return ok; }
    if(type=="bind"){ bool ok = parse_bind(s,r,out.bind); if(ok){ out.type=MsgType::Bind; } return ok; }
    return false;
}

bool parse_message(const std::string& json, ParsedMessage& out){ return parse_root(json, out); }

}

