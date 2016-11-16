#include <xsom.hpp>
#include <utility>
#include <vector>
#include <functional>
#include <algorithm>

#include <ccmpl.hpp>

using namespace std::placeholders;

typedef xsom::Point2D<double>     Pos;
typedef double                    Weight;
typedef double                    Input;
typedef double                    Act;
typedef xsom::tab2d::Table<Act>     Acts;
typedef xsom::tab2d::Table<Weight>  Weights;

Pos random_pos(void) {
  Pos p;

  do {
    p.x = -1 + std::rand()/(1.+RAND_MAX)*2;
    p.y = -1 + std::rand()/(1.+RAND_MAX)*2;
  } while(p*p > 1);

  return p;
}

bool inside_the_map(const Pos& p) {
  return p*p <= 1;
}

// This is a fake non-smooth activity function.

double activity(const Pos& p) {
  if(p*p < .25) {
    if(p.x*p.y > 0)
      return 1;
    else
      return 0;
  }

  if(p.x*p.y > 0)
    return 0;
  else
    return 1;
}

double merge(const Acts& acts,const Pos& p) {
  return acts(p);
}


void update_output(Pos& out,const Pos& request) {
  out = request;
}

void set_weight(std::function<Weight (const Pos&)> weight_at) {
  // We handle no weights here.
}

Weight get_weight(const Pos& p) {
  // We handle no weights here.
  return 0;
}

double match(const Pos& p, Weight w, double x) {
  // We handle no weights here.
  return 0;
}

double h(const Pos& p) {
  // This is dummy.
  return 0;
}


// Data transfer to plotting structures.

void fill_activity(std::vector<ccmpl::ValueAt>& points) {
  points.clear();
  for(unsigned int i=0; i < 5000; ++i) {
    auto p = random_pos();
    points.push_back({p.x,p.y,activity(p)});
  }
}

void fill_dot(const Pos& p, ccmpl::Point& pt) {
  pt.x = p.x;
  pt.y = p.y;
}

void fill_plot(const Acts& l, std::vector<ccmpl::ValueAt>& points) {
  points.clear();

  xsom::Index2D idx;
  unsigned int width  = l.mapping.size.w;
  unsigned int height = l.mapping.size.h;
  auto c   = l.content.begin();
  for(idx.h = 0; idx.h < height; ++idx.h)
    for(idx.w = 0; idx.w < width; ++idx.w, ++c) {
      xsom::Point2D<double> pos = l.mapping.index2pos(idx);
      if(l.pos_is_valid(pos))
	points.push_back({pos.x,pos.y,*c});
    }
}


// The space [-1,1]x[-1,1] is mapped to a [0,501[x[0,501[ grid, where
// the convolution is performed. The mapping provides conversion
// functions.

#define CONVOLUTION_SIDE 501
#define GRID_SIGMA       (CONVOLUTION_SIDE/20.0)

#define VIEW_FILE "viewer-000-004.py"
int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cout << "Usage : " << std::endl
  	      << argv[0] << " generate" << std::endl
  	      << argv[0] << " run | ./" << VIEW_FILE << std::endl;
    return 0;
  }

  bool generate_mode = std::string(argv[1])=="generate";

  ///////////
  //       //
  // State //
  //       //
  ///////////
  
  // Conversion from positions to row-column numbers.

  auto big_mapping 
    = xsom::tab2d::mapping<CONVOLUTION_SIDE,
			 CONVOLUTION_SIDE>({-1,-1},{1,1}); // min_pos, max_pos

  auto small_mapping 
    = xsom::tab2d::mapping<10,20>({-1,-1},{1,1}); 

  // This stores double values within a grid.
  auto grid = xsom::tab2d::table<double>(small_mapping,inside_the_map); 
  // grid initialization.
  for(auto& e : grid.content) e=0; 

  // This stores the convolution of some bi-dimentional signal. The
  // convolution is based on a tabular representation, since a
  // discrete Fourrier Transform is used.
  auto fft1 = xsom::tab2d::fft::convolution(GRID_SIGMA,xsom::tab::fft::Convolution::KernelType::Gaussian, big_mapping,inside_the_map);
  auto fft2 = xsom::tab2d::fft::convolution(GRID_SIGMA,xsom::tab::fft::Convolution::KernelType::Gaussian, big_mapping,inside_the_map);
  Pos win1;
  Pos win2;
  double input;

  //////////
  //      //
  // Plot //
  //      //
  //////////
  

  auto display      = ccmpl::layout(10,10, {"##","##"});
  display()         = {-1.,1.,-1.,1};
  display()         = ccmpl::show_tics(false,false);
  display().title   = "10x20 Grid";
  display()        += ccmpl::surface("cmap='gray'", 0, 1, std::bind(&xsom::tab2d::fill_plot_values, std::ref(grid), _1));
  display()        += ccmpl::dot("c='r',lw=1,s=20,zorder=2", std::bind(fill_dot, std::ref(win1), _1));
  display++;
  display()         = {-1.,1.,-1.,1};
  display()         = ccmpl::show_tics(false,false);
  display().title   = "501x501 convolution of the 10x10 grid";
  display()        += ccmpl::surface("cmap='jet'", 0, 1, std::bind(&xsom::tab2d::fft::Convolution::fill_plot, std::ref(fft1), _1));
  display++;
  display()         = {-1.,1.,-1.,1};
  display()         = ccmpl::show_tics(false,false);
  display().title   = "Raw values";
  display()        += ccmpl::surface("cmap='gray'", 0, 1, fill_activity);
  display()        += ccmpl::dot("c='r',lw=1,s=20,zorder=2", std::bind(fill_dot, std::ref(win2), _1));
  display++;
  display()         = {-1.,1.,-1.,1};
  display()         = ccmpl::show_tics(false,false);
  display().title   = "501x501 convolution raw values";
  display()        += ccmpl::surface("cmap='jet'", 0, 1, std::bind(&xsom::tab2d::fft::Convolution::fill_plot, std::ref(fft2), _1));
  
  if(generate_mode) {
    display.make_python(VIEW_FILE,true);
    return 0;
  }

  /////////
  //     //
  // Run //
  //     //
  /////////

  // Building layers and maps.

  auto layer
    = xsom::layer<Weight,
  		  Pos>(0,input,match,h,get_weight,
  		       std::bind(&xsom::tab2d::Table<double>::learn, std::ref(grid), _1),
  		       set_weight);

  auto map 
    = xsom::map(win1,
  		std::bind(&xsom::tab2d::fft::Convolution::learn, std::ref(fft1), _1),
  		std::bind(merge, std::ref(grid), _1),
  		std::bind(&xsom::tab2d::fft::Convolution::bmu, std::ref(fft1)),
  		update_output);


  auto archi = xsom::setup::network();
  // map   += layer;  I do not add the layer since the dummy match produces a 0 value.
  archi += map;
				   
  // Run the computation. 

  grid.learn(activity); // Let us set explicitly the grid values, since no grid evaluation is performed.
  archi.update();       // updating the maps calls fft1.bmu, which computes the convolution.

  // This is a convolution without any map/layer context.
  fft2.learn(activity);
  win2 = fft2.bmu();
	   	     
  std::cout << display("##""#""##""#","example-000-004.pdf", ccmpl::nofile())
  	    << ccmpl::stop;

  return 0;
}
