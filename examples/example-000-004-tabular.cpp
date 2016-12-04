#include <xsom.hpp>
#include <utility>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>

#include <ccmpl.hpp>

using namespace std::placeholders;

#define PI 3.1415926544


// This shows how to handle tabular data. xsom::tab2d::* classes and
// functions are used, since 2D arrays are considered. Similar code
// can also be obtained for 1D array thanks to xsom::tab1d::*.

// Let us define a parametrized surface, torus : (u,v) -> (x,y,z)

struct Point3D {
  double x,y,z;
};
// This comparison allows for the further use of bmu methods for a
// table containing Point3D values.
bool operator<(const Point3D& a, const Point3D& b) {
  return a.x < b.x; // Comparison based on x only...
}

#define torus_r  .5
#define torus_R  (1-torus_r)
Point3D torus(const xsom::Point2D<double>& uv) { // u,v in [0,2pi]
  double r = torus_R + torus_r*std::cos(uv.y);
  return {r*std::cos(uv.x), r*std::sin(uv.x), torus_r*std::sin(uv.y)};
}

ccmpl::RGB color_of_3d(const Point3D& p) {return {(p.x+1)/2, (p.y+1)/2, (p.z+torus_r)/(2*torus_r)};}
double     value_of_3d(const Point3D& p) {return  (p.x+1)/2                                       ;}

#define NB_U 200
#define NB_V  50

#define NB_TRIANGULATION_POINTS 1000

// let us also define a function, f : (u,v) -> value

double gabor(const xsom::Point2D<double>& uv) {
  return std::sin(20*uv.x)*std::exp(-5*(uv*uv));
}

// Let us play with tabular functions.

#define VIEW_FILE "viewer-000-004.py"
int main(int argc, char* argv[]) {
  if(argc < 2) {
    std::cout << "Usage : " << std::endl
  	      << argv[0] << " generate" << std::endl
  	      << argv[0] << " run | ./" << VIEW_FILE << std::endl;
    return 0;
  }

  bool generate_mode = std::string(argv[1])=="generate";

  auto display      = ccmpl::layout(20,10, {"####", "###."});
  std::string flags = "";

  // Let us start with the more general case. We will provide a
  // tabular version for the torus function. We define the function
  // for some (u,v) pairs only (the lambda function returns true for
  // them). Just omit this argument if you want to allow for values at
  // each position.
  
  auto torus_mapping = xsom::tab2d::mapping({0, 0}, {2*PI, 2*PI}, {NB_U, NB_V});
  auto tabular_torus = xsom::tab2d::table<Point3D>(torus_mapping,
						   [](const xsom::Point2D<double>& uv) {auto d = uv - xsom::Point2D<double>(PI,PI); return d*d <= PI*PI;});
  
  // We set the values in the grid from the torus function. The grid
  // contains Poin3D values. Those which are actually set are the ones
  // for which the above lambda-function returns true.  For the other
  // to be actually initialized, let us set first (and once) the value
  // of all the points in the grid.
  tabular_torus.clear({0,0,0}); // initialization to the (0,0,0) Point3D value.
  tabular_torus.learn(torus); 

  // This displays a color image from the table. We have to provide
  // fill_image_rgb with a function that converts the content of the
  // map (Point3D) into a color.

  display()         = {0, 2*PI, 0, 2*PI};
  display()         = ccmpl::show_tics(true, true);
  display().title   = "torus (clipped, fill_image_rgb)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        += ccmpl::image("interpolation='bilinear'",
  				   [&tabular_torus](std::vector<double>& x,
  						    std::vector<double>& y,
  						    std::vector<double>& z,
  						    unsigned int& width,
  						    unsigned int& depth) {tabular_torus.fill_image_rgb(color_of_3d, x, y, z, width, depth);});
  flags += '#';
  
  // The following does the same with a 1-depth image.
  display++;
  display()         = {0, 2*PI, 0, 2*PI};
  display()         = ccmpl::show_tics(true, true);
  display().title   = "torus (clipped, fill_image_gray)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        += ccmpl::image("cmap='jet', interpolation='bilinear', clim=(0,1)",
  				   [&tabular_torus](std::vector<double>& x,
  						    std::vector<double>& y,
  						    std::vector<double>& z,
  						    unsigned int& width,
  						    unsigned int& depth) {tabular_torus.fill_image_gray(value_of_3d, x, y, z, width, depth);});
  flags += '#';
  // If the content supports the < operator, the highest position in the map (called bmu) can be obtained.
  display()        += ccmpl::dot("c='w',lw=1,s=50,zorder=1", [&tabular_torus](ccmpl::Point& dot) {dot = tabular_torus.bmu();});
  flags += '#';

  // The two following does the same as rgb and 1-depth previous
  // images, but surface and palette are rather used. This can be
  // slower since a Delaunay Triangulation is computed by
  // matplotlib. Nevertheless, the clipping is rendered in a better
  // way.

  display++;
  display()         = {0, 2*PI, 0, 2*PI};
  display()         = ccmpl::show_tics(true, true);
  display().title   = "torus (clipped, fill_palette)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        += ccmpl::palette("", [&tabular_torus](std::vector<ccmpl::ColorAt>& points) {tabular_torus.fill_palette(color_of_3d, NB_TRIANGULATION_POINTS, points);});
  flags += '#';

  display++;
  display()         = {0, 2*PI, 0, 2*PI};
  display()         = ccmpl::show_tics(true, true);
  display().title   = "torus (clipped, fill_surface)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        += ccmpl::surface("cmap='jet'", 0, 1, [&tabular_torus](std::vector<ccmpl::ValueAt>& points) {tabular_torus.fill_surface(value_of_3d, NB_TRIANGULATION_POINTS, points);});
  flags += '#';

  // Let us now consider the more specific but more usual usual case
  // of tabular functions returning a float. Converter to floats, such
  // as value_of_3d, are now useless, which allows us to use method
  // pointers directly instead of lambda-functions.

  
  auto gabor_mapping = xsom::tab2d::mapping({-1, -1}, {1, 1}, {NB_U, NB_V});
  auto tabular_gabor = xsom::tab2d::table<double>(gabor_mapping);
  tabular_gabor.learn(gabor);

  display++;
  display()         = {-1, 1, -1, 1};
  display()         = ccmpl::show_tics(false, false);
  display().title   = "gabor (fill_image_gray)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        += ccmpl::image("cmap='jet', interpolation='bilinear', clim=(-1,1), zorder=0",
  				    std::bind(&xsom::tab2d::Table<double>::fill_image_gray, std::ref(tabular_gabor), _1, _2, _3, _4, _5));
  flags += "#";
  display()        += ccmpl::dot("c='w',lw=1,s=50,zorder=1", [&tabular_gabor](ccmpl::Point& dot) {dot = tabular_gabor.bmu();});
  flags += "#";

  display++;
  display()         = {-1, 1, -1, 1};
  display()         = ccmpl::show_tics(false, false);
  display().title   = "gabor (fill_surface)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        +=  ccmpl::surface("cmap='jet'", -1, 1, std::bind(&xsom::tab2d::Table<double>::fill_surface, std::ref(tabular_gabor), NB_TRIANGULATION_POINTS, _1));
  flags += '#';


  // Last, let us mention that tabular functions are indeed functions
  // f : pos -> content. They can thus be used as an argument of
  // learn.
  
  auto gabor_lowres_mapping = xsom::tab2d::mapping({-1, -1}, {1, 1}, {10, 10});
  auto tabular_gabor_lowres = xsom::tab2d::table<double>(gabor_lowres_mapping);
  tabular_gabor_lowres.learn(tabular_gabor);

  display++;
  display()         = {-1, 1, -1, 1};
  display()         = ccmpl::show_tics(false, false);
  display().title   = "gabor (lower resolution, fill_image_gray)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        += ccmpl::image("cmap='jet', interpolation='bilinear', clim=(-1,1), zorder=0",
  				    std::bind(&xsom::tab2d::Table<double>::fill_image_gray, std::ref(tabular_gabor_lowres), _1, _2, _3, _4, _5));
  flags += '#';
  
  if(generate_mode) {
    display.make_python(VIEW_FILE,false);
    return 0;
  }
  
  std::cout << display(flags,"example-000-004.pdf", ccmpl::nofile())
  	    << ccmpl::stop;
  return 0;
}



