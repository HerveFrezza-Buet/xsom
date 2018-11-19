#include <xsom.hpp>
#include <ccmpl.hpp>
#include <cmath>

using namespace std::placeholders;

// This shows the use of gridplot functionalities. It is suitable for plotting 2D Kohonen grids.

// let us define a mapping f from R^2 to R^2
xsom::Point2D<double> twirl(const xsom::Point2D<double>& p) {
  double r2    = p*p;
  double theta = std::atan2(p.y, p.x) + 3.141592654*std::exp(-2*r2);
  double r     = sqrt(r2);
  return {r*std::cos(theta), r*std::sin(theta)};
}


// Even is the map is MAP_SIZExMAP_SIZE, we will display a NB_WxNB_H grid. For each line in the grid, the values in the map are taken at each STEP interval.
#define MAP_SIZE 100
#define NB_W 15
#define NB_H 15
#define STEP 2


#define VIEW_PREFIX "viewer-002-003"
int main(int argc, char* argv[]) {
  
  ccmpl::Main m(argc,argv,VIEW_PREFIX);
  
  
  auto mapping = xsom::tab2d::mapping({-1, -1}, {1, 1}, {MAP_SIZE, MAP_SIZE});
  auto table   = xsom::tab2d::table<xsom::Point2D<double>>(mapping);
  table.learn(twirl);
  
  auto display      = ccmpl::layout(10,10, {"#"});

  display()         = ccmpl::view2d({-1.1, 1.1}, {-1.1, 1.1}, ccmpl::aspect::equal, ccmpl::span::placeholder);
  display()         = ccmpl::show_tics(true, true);
  display().title   = "Grid plot";
  display()         = "equal";
  display()        += ccmpl::lines("'b-', linewidth=1.0, zorder=0",
				   std::bind(xsom::tab2d::fill_lines<NB_W, NB_H, STEP>, std::ref(table), _1));
  display()        += ccmpl::line("'r-', linewidth=7.0, zorder=1",
				  std::bind(xsom::tab2d::fill_x, std::ref(table), 0.1, STEP, _1));
  display()        += ccmpl::line("'g-', linewidth=7.0, zorder=1",
				  std::bind(xsom::tab2d::fill_y, std::ref(table), 0.1, STEP, _1));
  

  // the ccmpl::Main object handles generation here
  m.generate(display, false); // true means "use GUI"
  
  std::cout << display("###", "example-002-003.pdf", ccmpl::nofile())
  	    << ccmpl::stop;
  return 0;
}
