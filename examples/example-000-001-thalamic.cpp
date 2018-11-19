
#include <xsom.hpp>
#include <array>
#include <functional>
#include <iostream>
#include <cstdlib>
#include <random>
#include <ctime>

#include <ccmpl.hpp>

using namespace std::placeholders;

#define MAP_SIZE         50
#define T_MATCH_SIGMA    .1
#define ALPHA            .1
#define BUMP_SIZE         3

#define T_MATCH_SIGMA_2 ((T_MATCH_SIGMA)*(T_MATCH_SIGMA))

typedef double                       Pos;
typedef xsom::Point2D<double>        Input;
typedef Input                        TWeight;
typedef double                       Act;
typedef std::array<Act,MAP_SIZE>     Acts; 
typedef std::array<TWeight,MAP_SIZE> TWeights;

template<typename X>
double euclidian(const X& x1,const X& x2) {
  X d = x1-x2;
  return d*d;
}

Pos idx2pos(unsigned int i) {
  return (Pos)(i+.5);
}

unsigned int pos2idx(Pos p) {
  if(p<0)
    return 0;
  unsigned int i = (unsigned int)(p);
  if(i>=MAP_SIZE)
    return MAP_SIZE-1;
  return i;
}

template<typename RANDOM_DEVICE>
Input sample(RANDOM_DEVICE& rd) {
  auto d = std::uniform_real_distribution<double>(0,1);
  return {d(rd), d(rd)};
}

double t_dist(const TWeight& w,const Input& xi) {
  return euclidian(w,xi);
}

double m_dist(const Pos& w1,const Pos& w2) {
  return euclidian(w1,w2);
}

double t_match(const Pos& local, const TWeight& w, const Input& x) {
  // local is not used since the matching rule do not depend on the
  // position of the unit where it is applied.
  return xsom::gaussian(w,x,T_MATCH_SIGMA_2,t_dist);
}

double h(const Pos& winner,const Pos& p) {
  return xsom::linear(p,winner,BUMP_SIZE,m_dist);
}

void set_act(Acts& act,
	     std::function<Act (const Pos&)> activity_at) {
  auto a = act.begin();
  for(unsigned int i=0; i<MAP_SIZE; ++i)
    *(a++) = activity_at(idx2pos(i));
}

void set_tweight(TWeights& wgt,
		 std::function<TWeight (const Pos&)> weight_at) {
  auto a = wgt.begin();
  for(unsigned int i=0; i<MAP_SIZE; ++i)
    *(a++) = weight_at(idx2pos(i));
}

TWeight get_tweight(const TWeights& wgt, const Pos& p) {
  return wgt[pos2idx(p)];
}

double merge(const Acts& thal, const Pos& p) {
  return thal[pos2idx(p)];
}

Pos winner(const Acts& act) {
  Act max = -1;
  Pos res = 0;
  Pos pos = 0;
  for(auto a : act) {
    if(a > max) {
      max = a;
      res = pos;
    }
    ++pos;
  }
  return res;
}

void update_output(Pos& out, const Pos& request) {
  out = request;
}

void fill_bar(const Pos& p, double& bar) {
  bar = p;
}

void fill_in(const Input& xi, ccmpl::Point& pt) {
  pt.x = xi.x;
  pt.y = xi.y;
}

void fill_win(const Pos& win, const TWeights& tw, ccmpl::Point& pt) {
  TWeight wwin = tw[pos2idx(win)];
  pt.x = wwin.x;
  pt.y = wwin.y;
}

void fill_wgt(const TWeights& tab, std::vector<ccmpl::Point>& curve) {
  curve.clear();
  for(auto v : tab) 
    curve.push_back({v.x,v.y});
}

void fill_act(const Acts& tab, std::vector<ccmpl::Point>& curve) {
  curve.clear();
  unsigned int i=0;
  for(auto v : tab) 
    curve.push_back({idx2pos(i++),v});
}

#define VIEW_PREFIX "viewer-000-001"
int main(int argc, char* argv[]) {
  
  // random seed initialization
  std::random_device rd;
  std::mt19937 gen(rd());
  
  ccmpl::Main m(argc,argv,VIEW_PREFIX);
  
  // Simulation state.

  Acts     ta;    // thalamic activities
  Acts     ma;    // map activities
  TWeights tw;    // thalamic weights
  Input    xi;    // the current input
  Pos      win;   // the winning position.

  // Plot
  
  // The layout string arguments represent with a # the position of
  // the two graphs within a grid, line by line.
  auto display = ccmpl::layout(5.0, 5.0,
			       {"#-","-#"},
			       ccmpl::RGB(1., 1., 1.));

  // We set the ratios width_ratios and height_ratios of gridspec. 
  // This function call is optional
  display.set_ratios({1., 1.}, {1., 1.});

  // Let us detail the first graph. Default limits are x=[0,1], y=[0,1]
  display().title   = "Kohonen SOM";
  display()         = ccmpl::view2d({0., 1.}, {0., 1.}, ccmpl::aspect::fit, ccmpl::span::placeholder); 
  display()        += ccmpl::line("'b-'", std::bind(fill_wgt,            std::ref(tw),                _1));
  display()        += ccmpl::dot ("c='r',lw=1,s=20", std::bind(fill_win, std::ref(win), std::ref(tw), _1));
  display()        += ccmpl::dot ("c='r',lw=1,s=20", std::bind(fill_in,  std::ref(xi),                _1));

  display++; // Skip to next graph.
  display()         = ccmpl::view2d({0., MAP_SIZE}, {0., 1.}, ccmpl::aspect::fit, ccmpl::span::placeholder); 
  display().title   = "Map activity";
  display()        += ccmpl::line("'b-'", std::bind(fill_act, std::ref(ta),  _1));
  display()        += ccmpl::vbar("'r-'", std::bind(fill_bar, std::ref(win), _1));

  // the ccmpl::Main object handles generation here
  m.generate(display, true); // true means "use GUI"
  
  // Execution

  for(auto& w : tw) w = sample(gen);

  auto thal 
    = xsom::layer<TWeight,
		  Pos>(ALPHA,
		       xi,                                         // Input
		       t_match,                                    // match(pos,w,xi) 
		       std::bind(h,           std::ref(win), _1),  // h(pos)
		       std::bind(get_tweight, std::ref(tw),  _1),  // weight(pos)
		       std::bind(set_act,     std::ref(ta),  _1),  // set_act(f), with f(pos) the desired activity value at pos.
		       std::bind(set_tweight, std::ref(tw),  _1)); // set_weight(f), with f(pos) the desired weight at pos.

  auto map
    = xsom::map(win,
  		std::bind(set_act, std::ref(ma), _1), // set_act(f), with f(pos) the desired activity value at pos.
  		std::bind(merge,   std::ref(ta), _1), // merge(p) = f(act1(p),act2(p),...actn(p))
  		std::bind(winner,  std::ref(ma)    ), // best() = argmax_p act(p)
  		update_output);                 


  auto archi = xsom::setup::network();
  map   += thal;
  archi += map;
    

  while(true) {
    
    // Update the architecture.
    
    xi = sample(gen);
    archi.update(); 
    archi.learn();
      
    // Fill the plot data and send it to the plotting script.  There
    // are 5 stuff to plot, on all axes. To update only the first and
    // the third in the plot, provide "#-#--". Here, we plot all of
    // them. Other args, if non "" (or nofile), are respectively the pdf name and
    // the png name.

    std::cout << display("###" "##", ccmpl::nofile(), ccmpl::nofile());
  }


  return 0;
}
