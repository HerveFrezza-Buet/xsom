
#include <xsom.hpp>
#include <ccmpl.hpp>
#include <vector>
#include <random>
#include <iostream>
#include <cmath>
#include <algorithm>
#include <functional>
#include <array>
#include <deque>

using namespace std::placeholders;

#define MAP_SIZE 500
#define T_MATCH_SIGMA .4
#define C_MATCH_SIGMA .4
#define ALPHA .1
#define H_RADIUS .05
#define SEQUENCE "AAAFFFBCDE"
#define BETA 0.5


// For the upper layer convolution.
#define KERNEL_TYPE xsom::tab::fft::KernelType::Gaussian  // We use a Gaussian convolution kernel
#define SIGMA_CONV .0125                                  // This is its standard deviation un Pos units.
#define GRID_SIGMA (SIGMA_CONV * MAP_SIZE)                // This is its standard deviation in table index units (since MAP_SIZE is mapped into [0,1]).

// Other macros
#define T_MATCH_SIGMA_2 ((T_MATCH_SIGMA)*(T_MATCH_SIGMA))
#define C_MATCH_SIGMA_2 ((C_MATCH_SIGMA)*(C_MATCH_SIGMA))



using Pos      = double; // The map is 1D, so unit positions are scalar. Here, it is a double (in [0, 1]).
using Input    = double; // Inputs are sequences of 1D values.
using TWeight  = Input;  // "Thalamic" weights (i.e weights matching the input).
using CWeight  = Pos;    // "Cortical" weights (i.e. weights matching the previous bmu position).

using Act      = double; // The map activities are scalar.

// An activity distribution over the map is an array of activities. 
// Let us use tabular tools from xsom (1D here)

using Acts     = xsom::tab1d::Table<Act>;
using CWeights = xsom::tab1d::Table<CWeight>;
using TWeights = xsom::tab1d::Table<TWeight>;

// As we will use convolution at the final layer, wee need a table that performs convolution.
using Map      = xsom::tab1d::fft::Convolution;

// This compute the squared euclidian distance ofr any type X which
// supports - and *.
template<typename X>
double euclidian(const X& x1,const X& x2) {
  X d = x1-x2;
  return d*d;
}

// This is the "thalamic" distance, that compares the input to the
// prototype.
double t_dist(const TWeight& w, const Input& xi) {return euclidian(w,xi);}

// This the 'cortical' distance, that compares a position in a remote
// map to a cortical weight.
double c_dist(const CWeight& w, const Input& p)  {return euclidian(w,p);}


// This tells how much the input x matches against the weight w. This
// functions returns 1 when the input and the prototype match
// perfectly, and is 0 when they do not match. The matching rule may
// depend on the position of the considered weight in the map
// (local). This is not the case here.
double t_match(const Pos& local, const TWeight& w, const Input& x) {return xsom::gaussian(w,x,T_MATCH_SIGMA_2,t_dist);}


// This tells how much the remote position p matches against the
// weight w. This functions returns 1 when the input and the prototype
// match perfectly, and is 0 when they do not match. The matching rule
// may depend on the position of the considered weight in the map
// (local). This is not the case here.
double c_match(const Pos& local, const CWeight& w, const Pos& p) {return xsom::gaussian(w,p,C_MATCH_SIGMA_2,c_dist);}

// This is the "map" distance, that compares two map positions.
double m_dist(const Pos& w1,const Pos& w2) {return euclidian(w1,w2);}

// This computes the learning strenght (in [0,1]) to be applyed at
// position p in the map, knowing that BMU position is winner.
double h(const Pos& bmu,const Pos& p) {return xsom::linear(p, bmu, H_RADIUS, m_dist);}

// This merges the "thalamic" and "cortical" matchings into a single
// scalar.
double merge(const Acts& thal, const Acts& cort, const Pos& p) {
  double t = thal(p);
  double c = cort(p);
  return BETA * t + (1 - BETA) * c;
}

// This updates the output position of the map from a position request
// (which is the BMU).
void update_output(Pos& out, const Pos& request) {out = request;}


class State{
private :
  std::size_t seq_size;
  
public:
  Input       x     = 0;        // The current input
  Pos         bmu_  = 0;        // The previous bmu
  Pos         bmu   = 0;        // The bmu position

  xsom::tab1d::Mapping mapping; // converts positions to table index.
  
  Acts        ta;               // thalamic matching
  TWeights    tw;               // thalamic weights
  Acts        ca;               // cortical matching
  CWeights    cw;               // cortical weights
  Map         ma;               // map activities

  std::deque<Pos> bmus;
  
  template<typename RANDOM_DEVICE>
  State(RANDOM_DEVICE& rd)
    : seq_size(std::string(SEQUENCE).size()),
      mapping(0, 1, MAP_SIZE),
      ta(mapping),
      tw(mapping),
      ca(mapping),
      cw(mapping),
      ma(mapping, GRID_SIGMA, KERNEL_TYPE) {

    // Random initialization.
    auto uniform = std::uniform_real_distribution<double>(0.0, 1.0);
    for(auto& c : tw.content) c = uniform(rd);
    for(auto& c : cw.content) c = uniform(rd);
  }

  void next(const Input& o) {
    // This is for next step computation.
    x  = o;
    bmu_ = bmu;
    bmus.push_front(bmu);
    if(bmus.size() > seq_size)
      bmus.pop_back();
  }
};



// Let us define an object providing a sequence of inputs in
// [0,1]. The sequence is periodical, and a noise can be added. Values
// are represented bu letters, A for 0, F for 1.
class InputSampler {
private:
  std::string::iterator current;

  void compute_display() {
    display = seq;
    for(auto& c : display) c += 32;
    display[std::distance(seq.begin(), current)] -= 32;
  }
  
public:
  std::string seq;
  std::string display;
  
  InputSampler(const std::string& sequence) : seq(sequence) {
    current = seq.begin();
    compute_display();
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
    
    compute_display();
    
    return res;
  }
};


// Plot functions

void fill_bar(const double& value, double& bar_pos) {bar_pos = value;}

void print_label(const std::string& label, ccmpl::Point& pos, std::string& text) {pos = ccmpl::Point(0.02, 0.9); text = label;}

void fill_bmus(std::deque<Pos>& data, std::vector<ccmpl::Point>& dots) {
  dots.clear();
  auto out = std::back_inserter(dots);
  for(auto d : data)
    *(out++) = {d, 1.0};
}


#define VIEWER_PREFIX "recurrent_SOM"


int main(int argc, char* argv[]) {

  // random seed initialization
  std::random_device rd;
  std::mt19937 gen(rd());
  
  State state(gen);
  InputSampler sequence(SEQUENCE);
  
  ccmpl::Main m(argc, argv, VIEWER_PREFIX);
  
  auto display = ccmpl::layout(12.0, 8.0,
                               {">.",
				"##",
				"##",
			        ">."},
                               ccmpl::RGB(1., 1., 1.));
  display.set_ratios({10., 10.},{1., .5, .5, .5});

  std::string flags = "";
  
  // Display global merge (raw and convoluated), and bars for bmu(t) and bmu(t-1) (dashed).
  display().title = "Merged output";
  display()   = ccmpl::view2d({0., 1.}, {0., 1.1}, ccmpl::aspect::fit, ccmpl::span::placeholder);
  display()  += ccmpl::vbar("'k'",         std::bind(fill_bar, std::ref(state.bmu), _1));                   flags += "#";
  display()  += ccmpl::vbar("'k--'",       std::bind(fill_bar, std::ref(state.bmu_), _1));                  flags += "#";
  display()  += ccmpl::line("'b'",         std::bind(&Map::fill_line, std::ref(state.ma), _1));             flags += "#";
  display()  += ccmpl::line("'b--'",       std::bind(&Map::fill_line_noconv, std::ref(state.ma), _1));      flags += "#";
  display()  += ccmpl::text("",            std::bind(print_label, std::ref(sequence.display), _1, _2));     flags += "#";
  display++;

  // Display thalamic matching and a bar for bmu(t).
  display().title = "Thalamic activity";
  display()   = ccmpl::view2d({0., 1.}, {0., 1.1}, ccmpl::aspect::fit, ccmpl::span::placeholder);
  display()  += ccmpl::vbar("'k'", std::bind(fill_bar, std::ref(state.bmu), _1));                           flags += "#";
  display()  += ccmpl::line("'r'", std::bind(&Acts::fill_line, std::ref(state.ta), _1));                    flags += "#";
  display++;

  // Display cortical matching and a bar for bmu(t-1).
  display().title = "Cortical activity";
  display()   = ccmpl::view2d({0., 1.}, {0., 1.1}, ccmpl::aspect::fit, ccmpl::span::placeholder);
  display()  += ccmpl::vbar("'k--'", std::bind(fill_bar, std::ref(state.bmu_), _1));                        flags += "#";
  display()  += ccmpl::line("'g'",   std::bind(&Acts::fill_line, std::ref(state.ca), _1));                  flags += "#";
  display++;

  // Display thalamic weights, a horizontal bar for the current input and a vertical bar for bmu(t).
  display().title  = "Thalamic weights";
  display().ytitle = "observation";
  display()   = ccmpl::view2d({0., 1.}, {0., 1.1}, ccmpl::aspect::fit, ccmpl::span::placeholder);
  display()  += ccmpl::vbar("'k'",              std::bind(fill_bar, std::ref(state.bmu), _1));              flags += "#";
  display()  += ccmpl::line("'k'",              std::bind(&TWeights::fill_line, std::ref(state.tw), _1));   flags += "#";
  display()  += ccmpl::hbar("'r', linewidth=3", std::bind(fill_bar, std::ref(state.x), _1));                flags += "#";
  display++;

  // Display cortical weights, a horizontal bar for the current input (bmu(t-1)) and a vertical bar for bmu(t-1) as well.
  display().title = "Cortical weights";
  display()   = ccmpl::view2d({0., 1.}, {0., 1.1}, ccmpl::aspect::fit, ccmpl::span::placeholder);
  display()  += ccmpl::vbar("'k--'",            std::bind(fill_bar, std::ref(state.bmu_), _1));             flags += "#";
  display()  += ccmpl::line("'k'",              std::bind(&CWeights::fill_line, std::ref(state.cw), _1));   flags += "#";
  display()  += ccmpl::hbar("'g', linewidth=3", std::bind(fill_bar, std::ref(state.bmu_), _1));             flags += "#";
  display++;

  // Display both weights.
  display().title = "All weights";
  display()   = ccmpl::view2d({0., 1.}, {0., 1.1}, ccmpl::aspect::fit, ccmpl::span::placeholder);
  display()  += ccmpl::line("'r'",              std::bind(&TWeights::fill_line, std::ref(state.tw), _1));   flags += "#";
  display()  += ccmpl::line("'g'",              std::bind(&TWeights::fill_line, std::ref(state.cw), _1));   flags += "#";
  display()  += ccmpl::dots("c='b', s=30",      std::bind(fill_bmus, std::ref(state.bmus), _1));            flags += "#";
  
  // the ccmpl::Main object handles generation here.
  m.generate(display, true); // true means "use GUI".
  
  auto tl = xsom::layer<TWeight, Pos>(ALPHA,
				      state.x,
				      t_match,
				      std::bind(h,                std::ref(state.bmu), _1),
				      std::bind(&TWeights::get,   std::ref(state.tw),  _1),
				      std::bind(&Acts::learn,     std::ref(state.ta),  _1),
				      std::bind(&TWeights::learn, std::ref(state.tw),  _1));

  auto cl = xsom::layer<CWeight, Pos>(ALPHA,
				      state.bmu_,
				      c_match,
				      std::bind(h,                std::ref(state.bmu), _1),
				      std::bind(&CWeights::get,   std::ref(state.cw),  _1),
				      std::bind(&Acts::learn,     std::ref(state.ca),  _1),
				      std::bind(&CWeights::learn, std::ref(state.cw),  _1));

  auto map = xsom::map(state.bmu,
		       std::bind(&Map::learn,               std::ref(state.ma), _1),
		       std::bind(merge, std::ref(state.ta), std::ref(state.ca), _1),
		       [&ma = state.ma, &gen]() {return ma.random_bmu(gen);},
		       // std::bind(&Map::random_bmu<std::mt19937>, std::ref(state.ma), std::ref(gen)),  // This is another formulation of the previous line.
		       update_output);
  
  auto archi  = xsom::setup::network(); // This is the global arcgitecture, containing a single map here.
  map        += tl;                     // We add the "thalamic" layer to the map.
  map        += cl;                     // We add the "cortical" layer to the map.
  archi      += map;                    // We add the map to the architecture (one map here).

  // Let us build a sequencer for synchronizing the computation
  auto seq    = xsom::setup::sequencer(archi, display);

  // See the doxygen documentation of xsom::setup::Sequencer for an
  // exhaustive list of sequencer functionalities.


  /* */ seq.__def("step");
  /* */   seq.__step([&state, &sequence](){state.next(sequence());});
  /* */   seq.__update_and_learn();
  /* */ seq.__fed();
  
  /* */ seq.__def("step by step");
  /* */   seq.__for(30);
  /* */     seq.__call("step");
  /* */     seq.__plot([&flags](){return flags;});
  /* */   seq.__rof();
  /* */ seq.__fed();

  /* */ seq.__def("chunk");
  /* */   seq.__for(10);
  /* */     seq.__call("step");
  /* */   seq.__rof();
  /* */   seq.__plot([&flags](){return flags;});
  /* */ seq.__fed();

  /* */ seq.__loop();
  /* */   seq.__call("step by step");
  /* */   seq.__for(10);
  /* */     seq.__call("chunk");
  /* */   seq.__rof();
  /* */ seq.__pool();

  // Now we can run the simulation.
  seq.run();
  
  return 0;
}

