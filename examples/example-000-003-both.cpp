
#include <xsom.hpp>
#include <array>
#include <deque>
#include <utility>
#include <functional>
#include <iostream>
#include <fstream>
#include <cstdlib>
#include <iostream>
#include <ctime>
#include <cmath>

#include <ccmpl.hpp>

using namespace std::placeholders;

#define MAP_SIZE          500
#define NB_UPDATES         50
#define T_MATCH_SIGMA      .1
#define C_MATCH_SIGMA      .1
#define ALPHA              .1
#define ALPHA_OUT          .1
#define BETA               .5
#define BUMP_RADIUS        .2
#define TRACE_SIZE         50
#define STEP_PERIOD        50
#define SAVE_PERIOD       100
#define MAX_STEP            4
#define INPUT_SEQUENCE    "ABCDEFEDCB"

#define T_MATCH_SIGMA_2 ((T_MATCH_SIGMA)*(T_MATCH_SIGMA))
#define C_MATCH_SIGMA_2 ((C_MATCH_SIGMA)*(C_MATCH_SIGMA))

typedef xsom::Point2D<double>  Pos;
typedef Pos                    CWeight;
typedef double                 Input;
typedef Input                  TWeight;
typedef double                 Act;

template<typename X>
double euclidian(const X& x1,const X& x2) {
  X d = x1-x2;
  return d*d;
}

Pos random_pos(void) {
  Pos p;

  do {
    p.x = -1 + std::rand()/(1.+RAND_MAX)*2;
    p.y = -1 + std::rand()/(1.+RAND_MAX)*2;
  } while(p*p > 1);

  return p;
}


// Let us define a Layer as a set of 2D positions, tossed randomly in a
// 1-radius disk. Maps are made of a stack of such layers.
template<typename CONTENT>
class Layer {
public:
  typedef std::pair<Pos,CONTENT> element_type;
  std::array<element_type,MAP_SIZE> data;

private:

  const element_type& locate(const Pos& pos) const {
    // This is unrealistically time consuming....
    auto it     = data.begin();
    auto _it    = data.begin();
    double dmin = euclidian(pos,it->first);
    double d;
    for(++it; it != data.end(); ++it)
      if((d = euclidian(pos,it->first)) < dmin) {
	dmin =  d;
	_it = it;
      }
    return *_it;
  }

public:

  Layer(void) {
    for(auto& e : data) e.first = random_pos();
  }

  CONTENT operator()(const Pos& pos) const {
    return locate(pos).second;
  }

  Pos bmu() const {
    auto it    = data.begin();
    Pos res    = it->first;
    double max = it->second;
    double m;
    for(++it; it != data.end(); ++it)
      if((m = it->second) > max) {
    	max =  m;
    	res  = it->first;
      }
    return res;
  }
};
  
typedef Layer<Act>             Acts;
typedef Layer<CWeight>         CWeights;
typedef Layer<TWeight>         TWeights;


std::ostream& operator<<(std::ostream& os, const CWeights& l) {
  for(auto& e : l.data)
    os << e.first.x << ' ' << e.first.y << ' ' << e.second.x << e.second.y << std::endl;
  return os;
}

std::istream& operator>>(std::istream& is, CWeights& l) {
  for(auto& e : l.data) 
    is >> e.first.x >> e.first.y >> e.second.x >> e.second.y;
  return is;
}


std::ostream& operator<<(std::ostream& os, const TWeights& l) {
  for(auto& e : l.data)
    os << e.first.x << ' ' << e.first.y << ' ' << e.second << std::endl;
  return os;
}

std::istream& operator>>(std::istream& is, TWeights& l) {
  for(auto& e : l.data)
    is >> e.first.x >> e.first.y >> e.second;
  return is;
}


// This handles the simulation state.
class State {
public:

  Acts        ta;   // thalamic matching
  TWeights    tw;   // thalamic weights
  Acts        ca;   // cortical matching
  CWeights    cw;   // cortical weights
  Acts        ma;   // map activities


  Acts        ta_;  // thalamic matching
  CWeights    tw_;  // thalamic weights (positions)
  Acts        ca_;  // cortical matching
  CWeights    cw_;  // cortical weights
  Acts        ma_;  // map activities

  Input       x;    // The current input
  Pos         x_;   // The current input (a position)
  Pos         win;  // The winning position
  Pos         win_; // The winning position

  void clear() {
    for(auto& p : tw.data)  p.second = std::rand()/(1.+RAND_MAX);
    for(auto& p : cw.data)  p.second = random_pos();
    for(auto& p : tw_.data) p.second = random_pos();
    for(auto& p : cw_.data) p.second = random_pos();
    win.x  = 0;
    win.y  = 0;
    win_.x = 0;
    win_.y = 0;
    x      = 0;
    x_.x   = 0;
    x_.y   = 0;
  }

  void next(const Input& o) {
    x  = o;
    x_ = win; // The input of map_ is the previous map position.
  }

  void load(std::string filename) {
    std::ifstream file;
    clear();
    file.open(filename.c_str());
    if(!file) {}
    else {
      file >> tw >> cw >> tw_ >> cw_;
      file.close();
    }
  }
  
  void save(std::string filename) const {
    std::ofstream file;
    file.open(filename.c_str());
    file << tw  << std::endl
	 << cw  << std::endl
	 << tw_ << std::endl
	 << cw_ << std::endl;
    file.close();
  }
  
};

// Let us define an object providing a sequence of inputs in
// [0,1]. The sequence is periodical, and a noise can be added. Values
// are represented bu letters, A for 0, F for 1.
class InputSampler {
private:
  std::string::iterator current;
public:
  std::string seq;
  InputSampler(const std::string& sequence) : seq(sequence) {
    current = seq.begin();
  }

  double operator()() {
    double res;
    switch(*current) {
    case 'A': res = 0.0; break;
    case 'B': res = 0.2; break;
    case 'C': res = 0.4; break;
    case 'D': res = 0.6; break;
    case 'E': res = 0.8; break;
    case 'F': res = 1.0; break;
    default : res = 0;
    }
    if(++current == seq.end())
      current = seq.begin();
    return res;
  }
};

double t_dist(const TWeight& w,const Input& xi) {
  return euclidian(w,xi);
}

double c_dist(const CWeight& w,const Pos& p) {
  return euclidian(w,p);
}

double m_dist(const Pos& w1,const Pos& w2) {
  return euclidian(w1,w2);
}

double t_match(const Pos& local, const TWeight& w, const Input& x) {
  // local is not used since the matching rule do not depend on the
  // position of the unit where it is applied.
  return xsom::gaussian(w,x,T_MATCH_SIGMA_2,t_dist);
}

double c_match(const Pos& local, const CWeight& w, const Pos& p) {
  // local is not used since the matching rule do not depend on the
  // position of the unit where it is applied.
  return xsom::gaussian(w,p,C_MATCH_SIGMA_2,c_dist);
}

double h(const Pos& winner,const Pos& p) {
  return xsom::linear(p,winner,BUMP_RADIUS,c_dist);
}

void set_act(Acts& act,
	     std::function<Act (const Pos&)> activity_at) {
  for(auto& e : act.data) e.second = activity_at(e.first);
}

void set_tweight(TWeights& wgt,
		 std::function<TWeight (const Pos&)> weight_at) {
  for(auto& e : wgt.data) e.second = weight_at(e.first);
}

TWeight get_tweight(const TWeights& wgt,const Pos& p) {
  return wgt(p);
}

void set_cweight(CWeights& wgt,
		 std::function<CWeight (const Pos&)> weight_at) {
  for(auto& e : wgt.data) e.second = weight_at(e.first);
}

CWeight get_cweight(const CWeights& wgt,const Pos& p) {
  return wgt(p);
}

double merge(const Acts& thal, const Acts& cort, const Pos& p) {
  double t = thal(p);
  double c = cort(p);
  return sqrt(t*(BETA*t+(1-BETA)*c));
}

void update_output(Pos& out,const Pos& request) {
  out += (request-out)*ALPHA_OUT;
}

// Data transfer to plotting structures.

void fill_srf(const Layer<double>& l, std::vector<ccmpl::ValueAt>& points) {
  points.clear();
  for(auto& e : l.data) points.push_back({e.first.x,e.first.y,e.second});
}

void fill_pal(const Layer<CWeight>& l, std::vector<ccmpl::ColorAt>& points) {
  points.clear();
  for(auto& e : l.data) 
    points.push_back({e.first.x,e.first.y,
	  ccmpl::color::from_point(e.second,
				  -1.0,1.0,
				  -1.0,1.0)});
}

void fill_dot(const Pos& p, ccmpl::Point& pt) {
  pt.x = p.x;
  pt.y = p.y;
}

void fill_trc(const std::deque<Pos>& trace, std::vector<ccmpl::Point>& points) {
  points.clear();
  for(auto& e : trace)
    points.push_back({e.x,e.y});
}

void fill_cod(std::vector<ccmpl::ColorAt>& points) {
  points.clear();
  for(int i=0; i<MAP_SIZE; ++i) {
    Pos s = random_pos();
    points.push_back({s.x,s.y,
	  ccmpl::color::from_point(s,
				  -1.0,1.0,
				  -1.0,1.0)});
  }
}

// The main function.
    
    
#define ACT_VIEW_FILE "viewer-act-000-003.py"
#define WGT_VIEW_FILE "viewer-wgt-000-003.py"
#define DATA_FILE "example-003.save"


int main(int argc, char* argv[]) {
  if(argc != 3) {
    std::cout << "Usage : " << std::endl
  	      << argv[0] << " act generate" << std::endl
  	      << argv[0] << " act run | ./" << ACT_VIEW_FILE << std::endl
  	      << argv[0] << " wgt generate" << std::endl 
  	      << argv[0] << " wgt run | ./" << WGT_VIEW_FILE << std::endl;
    return 0;
  }

  bool generate_mode = std::string(argv[2])=="generate";
  bool act_mode = std::string(argv[1])=="act";

  std::srand(std::time(0));


  // The state of the simulation

  State state;
  std::deque<Pos> trace; // history trace.

  // Plotting
  double xsize = 9.0;
  double ysize = 9.0;

  auto layout = {"##-","##-","--#"};
  if(act_mode) { 
    layout = {"###","###"}; 
    ysize = 6.0; 
  }

  auto display = ccmpl::layout(xsize, ysize, layout);
  if(act_mode) {
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder); 
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Input map";
    display().xtitle  = "thalamic activity";
    display()        += ccmpl::surface("cmap='jet'",  0, 1,   std::bind(fill_srf, std::ref(state.ta),   _1));
    display()        += ccmpl::dot    ("c='w',lw=1,s=20",     std::bind(fill_dot, std::ref(state.win),  _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder); 
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Input map";
    display().xtitle  = "cortical activity";
    display()        += ccmpl::surface("cmap='jet'",  0, 1,   std::bind(fill_srf, std::ref(state.ca),   _1));
    display()        += ccmpl::dot    ("c='k',lw=1,s=20",     std::bind(fill_dot, std::ref(state.win_), _1));
    display()        += ccmpl::dot    ("c='w',lw=1,s=20",     std::bind(fill_dot, std::ref(state.win),  _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder); 
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Input map";
    display().xtitle  = "map activity";
    display()        += ccmpl::surface("cmap='jet'",  0, 1,   std::bind(fill_srf, std::ref(state.ma),   _1));
    display()        += ccmpl::dot    ("c='w',lw=1,s=20",     std::bind(fill_dot, std::ref(state.win),  _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder);
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Recurrent map";
    display().xtitle  = "thalamic activity";
    display()        += ccmpl::surface("cmap='jet'",  0, 1,   std::bind(fill_srf, std::ref(state.ta_),  _1));
    display()        += ccmpl::dot    ("c='w',lw=1,s=20",     std::bind(fill_dot, std::ref(state.win_), _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder);
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Recurrent map";
    display().xtitle  = "cortical activity";
    display()        += ccmpl::surface("cmap='jet'",  0, 1,   std::bind(fill_srf, std::ref(state.ca_),  _1));
    display()        += ccmpl::dot    ("c='k',lw=1,s=20",     std::bind(fill_dot, std::ref(state.win),  _1));
    display()        += ccmpl::dot    ("c='w',lw=1,s=20",     std::bind(fill_dot, std::ref(state.win_), _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder);
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Recurrent map";
    display().xtitle  = "map activity";
    display()        += ccmpl::surface("cmap='jet'",  0, 1,   std::bind(fill_srf, std::ref(state.ma_),  _1));
    display()        += ccmpl::dot    ("c='w',lw=1,s=20",     std::bind(fill_dot, std::ref(state.win_), _1));
  } 
  else {
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder);
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Input map";
    display().xtitle  = "thalamic weight";
    display()        += ccmpl::surface("cmap='binary'", 0, 1, std::bind(fill_srf, std::ref(state.tw), _1));
    display()        += ccmpl::line   ("'ro-'",               std::bind(fill_trc, std::ref(trace),    _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder);
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Input map";
    display().xtitle  = "cortical weight";
    display()        += ccmpl::palette("",                    std::bind(fill_pal, std::ref(state.cw),   _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder);
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Cortical map";
    display().xtitle  = "thalamic weight";
    display()        += ccmpl::palette("",                    std::bind(fill_pal, std::ref(state.tw_),  _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder);
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Cortical map";
    display().xtitle  = "cortical weight";
    display()        += ccmpl::palette("",                    std::bind(fill_pal, std::ref(state.cw_),  _1));
    display++;
    display()         = ccmpl::view2d({-1.0, 1.0}, {-1.0, 1.0}, ccmpl::aspect::equal, ccmpl::span::placeholder);
    display()         = ccmpl::show_tics(false,false);
    display().title   = "Position to Color";  
    display()        += ccmpl::palette("",                    fill_cod                                     );
  }

  if(generate_mode) {
    if(act_mode)
      display.make_python(ACT_VIEW_FILE,true);
    else
      display.make_python(WGT_VIEW_FILE,true);
    return 0;
  }

  InputSampler sample(INPUT_SEQUENCE);
  
  // Building layers and maps.

  auto tl 
    = xsom::layer<TWeight,
  		  Pos>(ALPHA,
  		       state.x,
  		       t_match,
  		       std::bind(h,           std::ref(state.win), _1),
  		       std::bind(get_tweight, std::ref(state.tw),  _1),
  		       std::bind(set_act,     std::ref(state.ta),  _1),
  		       std::bind(set_tweight, std::ref(state.tw),  _1));

  auto cl 
    = xsom::layer<CWeight,
  		  Pos>(ALPHA,
  		       state.win_,
  		       c_match,
  		       std::bind(h,           std::ref(state.win), _1),
  		       std::bind(get_cweight, std::ref(state.cw),  _1),
  		       std::bind(set_act,     std::ref(state.ca),  _1),
  		       std::bind(set_cweight, std::ref(state.cw),  _1));

  auto map 
    = xsom::map(state.win,
  		std::bind(set_act,    std::ref(state.ma),                     _1),
  		std::bind(merge,      std::ref(state.ta), std::ref(state.ca), _1),
  		std::bind(&Acts::bmu, std::ref(state.ma)                        ),
  		update_output);

  auto tl_ 
    = xsom::layer<CWeight,
  		  Pos>(ALPHA,
  		       state.x_,
  		       c_match,
  		       std::bind(h,           std::ref(state.win_), _1),
  		       std::bind(get_cweight, std::ref(state.tw_),  _1),
  		       std::bind(set_act,     std::ref(state.ta_),  _1),
  		       std::bind(set_cweight, std::ref(state.tw_),  _1));

  auto cl_ 
    = xsom::layer<CWeight,
  		  Pos>(ALPHA,
  		       state.win,
  		       c_match,
  		       std::bind(h,           std::ref(state.win_),_1),
  		       std::bind(get_cweight, std::ref(state.cw_) ,_1),
  		       std::bind(set_act,     std::ref(state.ca_) ,_1),
  		       std::bind(set_cweight, std::ref(state.cw_) ,_1));

  auto map_ 
    = xsom::map(state.win_,
  		std::bind(set_act,    std::ref(state.ma_),                      _1),
  		std::bind(merge,      std::ref(state.ta_), std::ref(state.ca_), _1),
  		std::bind(&Acts::bmu, std::ref(state.ma_)                         ),
  		update_output);

  auto archi = xsom::setup::network();
  
  map    += tl;
  map    += cl;
  archi  += map;

  map_   += tl_;
  map_   += cl_;
  archi  += map_;


  // Running

  state.load(DATA_FILE);

  unsigned int step = 0;
  unsigned int ustep = 0;
  unsigned int wstep = 0;
  while(true) {
    
    state.next(sample());

    // Updates

    for(unsigned int u=0; u < NB_UPDATES; ++u, ++ustep) {

      // As x and x_ are constant for all the updates, thalamic layers
      // need only to be updated once.
      tl.update_locked  = (u != 0);
      tl_.update_locked = (u != 0);

      archi.update();

      if(act_mode) {
	if(u==0)
	  std::cout << display("-#""###""##""-#""###""##", ccmpl::nofile(), ccmpl::filename("update",ustep,"png"));
	else
	  std::cout << display("##""###""##""##""###""##", ccmpl::nofile(), ccmpl::filename("update",ustep,"png"));
      }
    }

    if(act_mode && step == MAX_STEP) {
      std::cout << ccmpl::stop;
      return 0;
    }

    // Update the trace

    trace.push_back(state.win);
    if(trace.size()>TRACE_SIZE)
      trace.pop_front();

    // Learn

    archi.learn();

    if(!act_mode && step % STEP_PERIOD == 0) {
      if(step == 0)
	std::cout << display("##""#""#""#""#", ccmpl::nofile(), ccmpl::filename("som",wstep,"png"));
      else
	std::cout << display("##""#""#""#""-", ccmpl::nofile(), ccmpl::filename("som",wstep,"png"));
      ++wstep;
    }

    if(step % SAVE_PERIOD == 0) {
      state.save(DATA_FILE);
      std::cerr << '\"' << DATA_FILE << "\" saved." << std::endl;
    }

    ++step;
  }

  return 0;
}
