

#include <xsom.hpp>
#include <array>
#include <functional>
#include <iostream>
#include <cstdlib>
#include <iostream>
#include <ctime>

#include <ccmpl.hpp>

using namespace std::placeholders;


#define MAP_SIZE      50
#define C_MATCH_SIGMA  5
#define ALPHA         .1
#define DURATION     100

#define C_MATCH_SIGMA_2 ((C_MATCH_SIGMA)*(C_MATCH_SIGMA))

typedef double                       Pos;
typedef Pos                          CWeight;
typedef CWeight                      Input;
typedef double                       Act;
typedef std::array<Act,MAP_SIZE>     Acts; 
typedef std::array<CWeight,MAP_SIZE> CWeights;

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

double c_dist(const CWeight& w,const Input& p) {
  return euclidian(w,p);
}

double m_dist(const Pos& w1,const Pos& w2) {
  return euclidian(w1,w2);
}

double c_match(const Pos& local, const CWeight& w, const Pos& p) {
  // local is not used since the matching rule do not depend on the
  // position of the unit where it is applied.
  return xsom::gaussian(w,p,C_MATCH_SIGMA_2,c_dist);
}

// Not relevant in this example since we do not learn.
double h(const Pos& p) {
  return 0;
}

double merge(const Acts& cort, const Pos& p) {
  return cort[pos2idx(p)];
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

void update_output(Pos& out,const Pos& request) {
  out += .1*(request-out);
}

void set_act(Acts& act,
	     std::function<Act (const Pos&)> activity_at) {
  auto a = act.begin();
  for(unsigned int i=0; i<MAP_SIZE; ++i)
    *(a++) = activity_at(idx2pos(i));
}

void set_cweight(CWeights& wgt,
		 std::function<CWeight (const Pos&)> weight_at) {
  auto a = wgt.begin();
  for(unsigned int i=0; i<MAP_SIZE; ++i)
    *(a++) = weight_at(idx2pos(i));
}

CWeight get_cweight(const CWeights& wgt, const Pos& p) {
  return wgt[pos2idx(p)];
}

void fill_bar(const Pos& p, double& bar) {
  bar = p;
}

void fill_x(const Pos& p, ccmpl::Point& pt) {
  pt.x = p;
  pt.y = 0;
}

void fill_wgt(const CWeights& tab, std::vector<ccmpl::Point>& curve) {
  curve.clear();
  unsigned int i=0;
  for(auto v : tab) 
    curve.push_back({idx2pos(i++),v/MAP_SIZE});
}

void fill_act(const Acts& tab, std::vector<ccmpl::Point>& curve) {
  curve.clear();
  unsigned int i=0;
  for(auto v : tab) 
    curve.push_back({idx2pos(i++),v});
}

#define VIEW_FILE "viewer-000-002.py"
int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cout << "Usage : " << std::endl
  	      << argv[0] << " generate" << std::endl
  	      << argv[0] << " run | ./" << VIEW_FILE << std::endl;
    return 0;
  }

  bool generate_mode = std::string(argv[1])=="generate";

  std::srand(std::time(0));
  
  // Simulation state.

  Acts     ca_1,ca_2;    // cortical activities
  Acts     ma_1,ma_2;    // map activities
  CWeights cw_1,cw_2;    // cortical weights
  Pos      win_1,win_2;  // map outputs
  
  // Plot
  
  auto display     = ccmpl::layout(10.0, 3.0, {"##"});
  display()        = ccmpl::view2d({0., MAP_SIZE}, {0., 1.}, ccmpl::aspect::equal , ccmpl::span::placeholder);
  display().title  = "map 1";
  display()       += ccmpl::line("'k--'",            std::bind(fill_wgt, std::ref(cw_1),  _1));
  display()       += ccmpl::line("'b-'",             std::bind(fill_act, std::ref(ca_1),  _1));
  display()       += ccmpl::vbar("'r-'",             std::bind(fill_bar, std::ref(win_1), _1));
  display()       += ccmpl::dot ("c='r',lw=1,s=20",  std::bind(fill_x,   std::ref(win_2), _1));
  display++;
  display()        = ccmpl::view2d({0., MAP_SIZE}, {0., 1.}, ccmpl::aspect::equal , ccmpl::span::placeholder);
  display().title  = "map 2";
  display()       += ccmpl::line("'k--'",            std::bind(fill_wgt, std::ref(cw_2),  _1));
  display()       += ccmpl::line("'b-'",             std::bind(fill_act, std::ref(ca_2),  _1));
  display()       += ccmpl::vbar("'r-'",             std::bind(fill_bar, std::ref(win_2), _1));
  display()       += ccmpl::dot ("c='r',lw=1,s=20",  std::bind(fill_x,   std::ref(win_1), _1));

  if(generate_mode) {
    display.make_python(VIEW_FILE,true);
    return 0;
  }


  // Execution

  for(auto& a : ca_1) a=0;
  for(auto& a : ca_2) a=0;
  for(auto& a : ma_1) a=0;
  for(auto& a : ma_2) a=0;

  double i;
  i= 0;         for(auto& w : cw_1) w=idx2pos(i++);
  i=MAP_SIZE-1; for(auto& w : cw_2) w=idx2pos(i--);

  auto cl_1 
    = xsom::layer<CWeight,
  		  Pos>(ALPHA, 
  		       win_2,                                       // remote w*
  		       c_match,                                     // match(w,pos) 
  		       h,                                           // h(pos)... unused in this example
  		       std::bind(get_cweight, std::ref(cw_1), _1),  // weight(pos)
  		       std::bind(set_act,     std::ref(ca_1), _1),  // set_act(f), with f(pos) the desired activity value at pos.
  		       std::bind(set_cweight, std::ref(cw_1), _1)); // set_weight(f), with f(pos) the desired weight at pos.

  auto cl_2
    = xsom::layer<CWeight,
  		  Pos>(ALPHA,
  		       win_1,
  		       c_match,
  		       h,
  		       std::bind(get_cweight, std::ref(cw_2), _1),
  		       std::bind(set_act,     std::ref(ca_2), _1),
  		       std::bind(set_cweight, std::ref(cw_2), _1));
		
  auto map_1
    = xsom::map(win_1,
  		std::bind(set_act, std::ref(ma_1), _1), // set_act(f), with f(pos) the desired activity value at pos.
  		std::bind(merge,   std::ref(ca_1), _1), // merge(p) = f(act1(p),act2(p),...actn(p))
  		std::bind(winner,  std::ref(ma_1)    ), // best() = argmax_p act(p).
  		update_output);

  auto map_2
    = xsom::map(win_2,
  		std::bind(set_act, std::ref(ma_2), _1),
  		std::bind(merge,   std::ref(ca_2), _1),
  		std::bind(winner,  std::ref(ma_2)    ),
  		update_output);


  auto archi = xsom::setup::network();

  map_1 += cl_1;
  map_2 += cl_2;
  archi += map_1;
  archi += map_2;

  win_1 = idx2pos(0);
  win_2 = idx2pos(MAP_SIZE-1);
    
  int e=0;
  int f=0;
  while(true) {

    // Fill the plot data and send it to the plotting script.

    if(e<50)
      std::cout << display("####""####", ccmpl::nofile(), ccmpl::filename("frame",e,"png"));
    else 
      std::cout << display("####""####", ccmpl::nofile(), ccmpl::nofile());

    // Update the architecture.
    
    archi.update();
    if(f >= 50) {
      f=0;
      win_1 = std::rand()/(1.+RAND_MAX)*MAP_SIZE;
      win_2 = std::rand()/(1.+RAND_MAX)*MAP_SIZE;
    }

    if(e>50) ++f;
    ++e;
  }
  
  return 0;
}
