
#include <xsom.hpp>
#include <array>
#include <functional>
#include <iostream>
#include <random>
#include <ctime>
#include <algorithm>

#include <ccmpl.hpp>

using namespace std::placeholders;

#define PERIOD           50

#define MAP_SIZE         50
#define T_MATCH_SIGMA    .1
#define ALPHA            .1
#define BUMP_SIZE         3

#define T_MATCH_SIGMA_2 ((T_MATCH_SIGMA)*(T_MATCH_SIGMA))


// The map is 1D, so unit positions are scalar. Here, it is a double (in [0, 50]).
typedef double                        Pos;

// The inputs will be points in [0,1]^2
typedef xsom::Point2D<double>         Input;

// The prototypes (thalamic weights) are identical to inputs. So the
// type of prototypes is the same as the type of inputs.
typedef Input                         TWeight;

// The map activities are scalars. Here, an activity is a double.
typedef double                        Act;

// An activity distribution over the map is an array of activities.
typedef std::array<Act, MAP_SIZE>     Acts;

// The prototypes of the map is an array of weights.
typedef std::array<TWeight, MAP_SIZE> TWeights;

// This compute the squared euclidian distance ofr any type X which
// supports - and *.
template<typename X>
double euclidian(const X& x1,const X& x2) {
  X d = x1-x2;
  return d*d;
}

// Positions in the map are indexed with a double (i.e. POS), but the
// map is handled by an array, so only MAP_SIZE values are stored in
// arrays. For any continuous position, there is an index value in the
// array corresponding to it. The two following functions compute the
// conversion from scalar positions to array indices.

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

// This tosses a random input sample in [0,1]^2
template<typename RANDOM_DEVICE>
Input sample(RANDOM_DEVICE& rd) {
  auto d = std::uniform_real_distribution<double>(0,1);
  return {d(rd), d(rd)};
}

// This is the "thalamic" distance, that compares the input to the prototype.
double t_dist(const TWeight& w,const Input& xi) {
  return euclidian(w,xi);
}

// This is the "map" distance, that compares two map positions.
double m_dist(const Pos& w1,const Pos& w2) {
  return euclidian(w1,w2);
}

// This tells how much the input x matches against the weight w. This
// functions returns 1 when the input and the prototype match
// perfectly, and is 0 when they do not match. The matching rule may
// depend on the position of the considered weight in the map
// (local). This is not the case here.
double t_match(const Pos& local, const TWeight& w, const Input& x) {
  return xsom::gaussian(w,x,T_MATCH_SIGMA_2,t_dist);
}

// This computes the learning strenght (in [0,1]) to be applyed at
// position p in the map, knowing that BMU position is winner.
double h(const Pos& winner,const Pos& p) {
  return xsom::linear(p,winner,BUMP_SIZE,m_dist);
}

// This sets act to be equal to a function, i.e. each value stored in
// act is the result of the function for the corresponding position.
void set_act(Acts& act,
	     std::function<Act (const Pos&)> activity_at) {
  auto a = act.begin();
  for(unsigned int i=0; i<MAP_SIZE; ++i)
    *(a++) = activity_at(idx2pos(i));
}

// This is similar to set_act for weights.
void set_tweight(TWeights& wgt,
		 std::function<TWeight (const Pos&)> weight_at) {
  auto a = wgt.begin();
  for(unsigned int i=0; i<MAP_SIZE; ++i)
    *(a++) = weight_at(idx2pos(i));
}

// This returns the weight value stored in wgt at a given position.
TWeight get_tweight(const TWeights& wgt, const Pos& p) {
  return wgt[pos2idx(p)];
}

// This function takes many activity arrays as argument and tell how
// to compute a global activity at that position. Here, we consider a
// single activity array, so this merging only consists in returning
// the value of that single array at this position.
double merge(const Acts& thal, const Pos& p) {
  return thal[pos2idx(p)];
}

// This compute the best matching unit position from an array of value.
Pos winner(const Acts& act) {
  auto max_iter = std::max_element(act.begin(), act.end());
  return idx2pos(std::distance(act.begin(), max_iter));
}

// The ouput of the map is a 'selected position', determined from the
// activities. More precisely, the activites suggest a position
// (i.e. the request), that is usually the position where the merge is
// maximal, and the actual map output has to be updated from this
// questested position. Here, this update is an affectation, so the
// output of the map is directly the BMU position.
void update_output(Pos& out, const Pos& request) {
  out = request;
}

// The following are ccmpl graphics filling functions.

// This shows the output (BMU) position in the activity distribution visualization.
void fill_bar(const Pos& p, double& bar) {
  bar = p;
}

// This shows where the input is.
void fill_in(const Input& xi, ccmpl::Point& pt) {
  pt.x = xi.x;
  pt.y = xi.y;
}

// This shows the BMU weight.
void fill_win(const Pos& win, const TWeights& tw, ccmpl::Point& pt) {
  TWeight wwin = tw[pos2idx(win)];
  pt.x = wwin.x;
  pt.y = wwin.y;
}

// This show the weights (this is a curve in [0, 1]^2).
void fill_wgt(const TWeights& tab, std::vector<ccmpl::Point>& curve) {
  curve.clear();
  for(auto v : tab) 
    curve.push_back({v.x,v.y});
}

// This displays the activity distribution.
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
  
  // The simulation consists in updating some values.... these values
  // are the simulation state. Here, the values are the following.

  Acts     ta;    // thalamic activities
  Acts     ma;    // map activities
  TWeights tw;    // thalamic weights
  Input    xi;    // the current input
  Pos      win;   // the winning position.
  
  for(auto& w : tw) w = sample(gen); // We initialize the weights randomly.

  // Next is the plotting.
  
  auto display = ccmpl::layout(5.0, 5.0,
			       {"#-","-#"},
			       ccmpl::RGB(1., 1., 1.));
  display.set_ratios({1., 1.}, {1., 1.});

  display().title   = "Kohonen SOM";
  display()         = ccmpl::view2d({0., 1.}, {0., 1.}, ccmpl::aspect::fit, ccmpl::span::placeholder); 
  display()        += ccmpl::line("'b-'", std::bind(fill_wgt,            std::ref(tw),                _1));
  display()        += ccmpl::dot ("c='r',lw=1,s=20", std::bind(fill_win, std::ref(win), std::ref(tw), _1));
  display()        += ccmpl::dot ("c='r',lw=1,s=20", std::bind(fill_in,  std::ref(xi),                _1));

  display++;
  
  display()         = ccmpl::view2d({0., MAP_SIZE}, {0., 1.}, ccmpl::aspect::fit, ccmpl::span::placeholder); 
  display().title   = "Map activity";
  display()        += ccmpl::line("'b-'", std::bind(fill_act, std::ref(ta),  _1));
  display()        += ccmpl::vbar("'r-'", std::bind(fill_bar, std::ref(win), _1));

  // the ccmpl::Main object handles generation here
  m.generate(display, true); // true means "use GUI"
  
  
  // Let us define an architecure that handles all the computation.

  // This is an activity layer. Maps are made of layers.
  auto thal 
    = xsom::layer<TWeight,                                         // Internal computation
		  Pos>(ALPHA,                                      // --------------------
		       xi,                                         // Input
		       t_match,                                    // match(pos,w,xi) 
		       std::bind(h,           std::ref(win), _1),  // h(pos)
		       std::bind(get_tweight, std::ref(tw),  _1),  // weight(pos)
		       std::bind(set_act,     std::ref(ta),  _1),  // set_act(f), with f(pos) the desired activity value at pos.
		       std::bind(set_tweight, std::ref(tw),  _1)); // set_weight(f), with f(pos) the desired weight at pos.

  // This is our map. It gathers here only one layer, merges the
  // result, and updates its ouputs. The layer will be added to the
  // map just next.
  auto map                                            // Internal computation
    = xsom::map(win,                                  // --------------------
  		std::bind(set_act, std::ref(ma), _1), // set_act(f), with f(pos) the desired activity value at pos.
  		std::bind(merge,   std::ref(ta), _1), // merge(p) = f(act1(p),act2(p),...actn(p))
  		std::bind(winner,  std::ref(ma)    ), // best() = argmax_p act(p)
  		update_output);                       // update(output, best) changes the output.


  auto archi = xsom::setup::network(); // This is the global arcgitecture, containing a single map here.
  map   += thal;                       // We add the layers to the maps (one layer into one map here).
  archi += map;                        // We add the maps to the architecture (one map here).
    

  // Now, we can run a simulation.
  
  while(true) {
    
    // Update the architecture.

    for(unsigned int i = 0; i < PERIOD; ++i) {
      xi = sample(gen);
      archi.update(); // This computes all the matchings.
      archi.learn();  // Then, weights are updated.
    }
    
    // Ths sends the display information to the ccmpl viewer (every PERIOD updates).
    std::cout << display("###" "##", ccmpl::nofile(), ccmpl::nofile());
  }


  return 0;
}
