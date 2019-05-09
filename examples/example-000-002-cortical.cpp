

#include <xsom.hpp>
#include <array>
#include <functional>
#include <iostream>
#include <iostream>
#include <ctime>
#include <random>
#include <algorithm>

#include <ccmpl.hpp>

using namespace std::placeholders;


#define MAP_SIZE      50
#define C_MATCH_SIGMA  5
#define ALPHA         .1
#define DURATION     100

#define C_MATCH_SIGMA_2 ((C_MATCH_SIGMA)*(C_MATCH_SIGMA))

// The map is 1D, so unit positions are scalar. Here, it is a double (in [0, 50]).
typedef double                        Pos;

// The 'cortical' weight is a position (in another map).
typedef Pos                           CWeight;

// The inputs for learning cortical weights are a position as well.
typedef Pos                           Input;

// The map activities are scalars. Here, an activity is a double.
typedef double                        Act;

// An activity distribution over the map is an array of activities.
typedef std::array<Act, MAP_SIZE>     Acts;

// The prototypes of a map is an array of cortical weights.
typedef std::array<CWeight, MAP_SIZE> CWeights;


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

// This the 'cortical' distance, that compares a position in a remote
// map to a cortical weight.
double c_dist(const CWeight& w,const Input& p) {
  return euclidian(w,p);
}

// This is the "map" distance, that compares two map positions.
double m_dist(const Pos& w1,const Pos& w2) {
  return euclidian(w1,w2);
}

// This tells how much the remote position p matches against the
// weight w. This functions returns 1 when the input and the prototype
// match perfectly, and is 0 when they do not match. The matching rule
// may depend on the position of the considered weight in the map
// (local). This is not the case here.
double c_match(const Pos& local, const CWeight& w, const Pos& p) {
  return xsom::gaussian(w,p,C_MATCH_SIGMA_2,c_dist);
}

// This computes the learning strenght (in [0,1]) to be applyed at
// position p in the map, knowing that BMU position is winner.
double h(const Pos& p) {
  return 0;
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
void set_cweight(CWeights& wgt,
		 std::function<CWeight (const Pos&)> weight_at) {
  auto a = wgt.begin();
  for(unsigned int i=0; i<MAP_SIZE; ++i)
    *(a++) = weight_at(idx2pos(i));
}

// This returns the weight value stored in wgt at a given position.
CWeight get_cweight(const CWeights& wgt, const Pos& p) {
  return wgt[pos2idx(p)];
}

// This function takes many activity arrays as argument and tell how
// to compute a global activity at that position. Here, we consider a
// single activity array, so this merging only consists in returning
// the value of that single array at this position.
double merge(const Acts& cort, const Pos& p) {
  return cort[pos2idx(p)];
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
// questested position. Here, the update consists of making the input
// slightly closer to the request.
void update_output(Pos& out,const Pos& request) {
  out += .1*(request-out);
}

// The following are ccmpl graphics filling functions.

// This shows the output position (which is not the BMU) in the
// activity distribution visualization.
void fill_bar(const Pos& p, double& bar) {
  bar = p;
}

// Shows the remote output as a dot.
void fill_x(const Pos& p, ccmpl::Point& pt) {
  pt.x = p;
  pt.y = 0;
}

// This displays the weights.
void fill_wgt(const CWeights& tab, std::vector<ccmpl::Point>& curve) {
  curve.clear();
  unsigned int i=0;
  for(auto v : tab) 
    curve.push_back({idx2pos(i++),v/MAP_SIZE});
}

// This displays the activity distribution.
void fill_act(const Acts& tab, std::vector<ccmpl::Point>& curve) {
  curve.clear();
  unsigned int i=0;
  for(auto v : tab) 
    curve.push_back({idx2pos(i++),v});
}

#define VIEW_PREFIX "viewer-000-002"
int main(int argc, char* argv[]) {
  
  // random seed initialization
  std::random_device rd;
  std::mt19937 gen(rd());
  
  ccmpl::Main m(argc,argv,VIEW_PREFIX);
  
  
  // The simulation consists in updating some values.... these values
  // are the simulation state. We consider two maps (_1 and _2
  // suffixes). Here, the values are the following.

  Acts     ca_1,ca_2;    // cortical activities
  Acts     ma_1,ma_2;    // map activities
  CWeights cw_1,cw_2;    // cortical weights
  Pos      win_1,win_2;  // map outputs

  for(auto& a : ca_1) a=0;
  for(auto& a : ca_2) a=0;
  for(auto& a : ma_1) a=0;
  for(auto& a : ma_2) a=0;

  // Cortical weights are set beforehand, thay won't be changed
  // (i.e. no learning here).
  unsigned int i;
  i= 0;         for(auto& w : cw_1) w=idx2pos(i++);
  i=MAP_SIZE-1; for(auto& w : cw_2) w=idx2pos(i--);

  // Let us initialize winner positions in each map.
  win_1 = idx2pos(0);
  win_2 = idx2pos(MAP_SIZE-1);
  
  // Next is the plotting.
  
  auto display     = ccmpl::layout(m.hostname, m.port,
				   10.0, 3.0, {"##"},
				   ccmpl::RGB(1., 1., 1.));
  
  display()        = ccmpl::view2d({0., MAP_SIZE}, {0., 1.}, ccmpl::aspect::fit , ccmpl::span::placeholder);
  display().title  = "map 1";
  display()       += ccmpl::line("'k--'",            std::bind(fill_wgt, std::ref(cw_1),  _1));
  display()       += ccmpl::line("'b-'",             std::bind(fill_act, std::ref(ca_1),  _1));
  display()       += ccmpl::vbar("'r-'",             std::bind(fill_bar, std::ref(win_1), _1));
  display()       += ccmpl::dot ("c='r',lw=1,s=20",  std::bind(fill_x,   std::ref(win_2), _1));
  
  display++;
  
  display()        = ccmpl::view2d({0., MAP_SIZE}, {0., 1.}, ccmpl::aspect::fit , ccmpl::span::placeholder);
  display().title  = "map 2";
  display()       += ccmpl::line("'k--'",            std::bind(fill_wgt, std::ref(cw_2),  _1));
  display()       += ccmpl::line("'b-'",             std::bind(fill_act, std::ref(ca_2),  _1));
  display()       += ccmpl::vbar("'r-'",             std::bind(fill_bar, std::ref(win_2), _1));
  display()       += ccmpl::dot ("c='r',lw=1,s=20",  std::bind(fill_x,   std::ref(win_1), _1));


  // the ccmpl::Main object handles generation here
  m.generate(display, true); // true means "use GUI"
  
  // Let us define an architecure that handles all the computation.

  // This is the activity layer for map #1
  auto cl_1 
    = xsom::layer<CWeight,                                          // Internal computation
  		  Pos>(ALPHA,                                       // --------------------
  		       win_2,                                       // The input is the remote w*
  		       c_match,                                     // match(pos, w, win_2) 
  		       h,                                           // h(pos)
  		       std::bind(get_cweight, std::ref(cw_1), _1),  // weight(pos)
  		       std::bind(set_act,     std::ref(ca_1), _1),  // set_act(f), with f(pos) the desired activity value at pos.
  		       std::bind(set_cweight, std::ref(cw_1), _1)); // set_weight(f), with f(pos) the desired weight at pos.

  // This is the activity layer for map #2 
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
  		update_output);                         // update(output, best) changes the output.

  auto map_2
    = xsom::map(win_2,
  		std::bind(set_act, std::ref(ma_2), _1),
  		std::bind(merge,   std::ref(ca_2), _1),
  		std::bind(winner,  std::ref(ma_2)    ),
  		update_output);


  auto archi = xsom::setup::network();

  map_1 += cl_1;
  archi += map_1;
  
  map_2 += cl_2;
  archi += map_2;

    
  int e=0;
  int f=0;
  auto u = std::uniform_real_distribution<double>(idx2pos(0), idx2pos(MAP_SIZE-1));
  while(true) {

    // Fill the plot data and send it to the plotting script.

    if(e<50)
      display("####""####", ccmpl::nofile(), ccmpl::filename("frame",e,"png"));
    else 
      display("####""####", ccmpl::nofile(), ccmpl::nofile());

    // Update the architecture.
    
    archi.update();
    if(f >= 50) {
      f=0;
      win_1 = u(gen);
      win_2 = u(gen);
    }

    if(e>50) ++f;
    ++e;
  }
  
  return 0;
}
