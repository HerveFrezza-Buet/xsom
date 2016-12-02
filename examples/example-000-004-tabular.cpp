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

#define torus_r  .5
#define torus_R  (1-torus_r)
Point3D torus(double u,double v) { // u,v in [0,2pi]
  double r = torus_R + torus_r*std::cos(v);
  return {r*std::cos(u), r*std::sin(u), torus_r*std::sin(v)};
}

ccmpl::RGB color_of_3d(const Point3D& p) {return {(p.x+1)/2, (p.y+1)/2, (p.z+torus_r)/(2*torus_r)};}
double     value_of_3d(const Point3D& p) {return  (p.x+1)/2                                       ;}

#define NB_U_BIG 100
#define NB_V_BIG  30

// let us also define a function, f : (u,v) -> value

double v(const xsom::Point2D<double>& uv) {
  if(uv*uv < .25) {
    if(uv.x*uv.y > 0)
      return 1;
    else
      return 0;
  }
  if(uv.x*uv.y > 0)
    return 0;
  else
    return 1;
}

// Let us play with tabular functions.

#define VIEW_FILE "viewer-000-004.py"
int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cout << "Usage : " << std::endl
  	      << argv[0] << " generate" << std::endl
  	      << argv[0] << " run | ./" << VIEW_FILE << std::endl;
    return 0;
  }

  bool generate_mode = std::string(argv[1])=="generate";

  auto display      = ccmpl::layout(10,5, {"##"});
  std::string flags = "";

  // Let us start with the more general case. We will provide a
  // tabular version for the torus function. We define the function
  // for all (u,v) pairs (the lambda function always return true).
  
  auto torus_mapping = xsom::tab2d::mapping({0, 0}, {2*PI, 2*PI}, {NB_U_BIG, NB_V_BIG});
  auto tabular_torus = xsom::tab2d::table<Point3D>(torus_mapping,
						   [](const xsom::Point2D<double>& uv) {return true;});
  
  // We set the values in the grid from the torus function. The grid contains Poin3D values.
  tabular_torus.learn([](const xsom::Point2D<double>& p) {return torus(p.x, p.y);}); 

  // This displays a color image from the table. We have to provide
  // fill_image_rgb with a function that converts the content of the
  // map (Point3D) into a color.
 
  display()         = {0, 2*PI, 0, 2*PI};
  display()         = ccmpl::show_tics(true, true);
  display().title   = "torus";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        +=  ccmpl::image("interpolation='bilinear'",
  				    [&tabular_torus](std::vector<double>& x,
  						     std::vector<double>& y,
  						     std::vector<double>& z,
  						     unsigned int& width,
  						     unsigned int& depth) {tabular_torus.fill_image_rgb(color_of_3d, x, y, z, width, depth);});  flags += '#';
  
  // The following does the same with a 1-depth image.
  
  display++;
  display()         = {0, 2*PI, 0, 2*PI};
  display()         = ccmpl::show_tics(true, true);
  display().title   = "torus";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        +=  ccmpl::image("cmap='jet', interpolation='bilinear', clim=(0,1)",
				    [&tabular_torus](std::vector<double>& x,
						     std::vector<double>& y,
						     std::vector<double>& z,
						     unsigned int& width,
						     unsigned int& depth) {tabular_torus.fill_image_gray(value_of_3d, x, y, z, width, depth);});  flags += '#';

  
  if(generate_mode) {
    display.make_python(VIEW_FILE,false);
    return 0;
  }
  
  std::cout << display(flags,"example-000-004.pdf", ccmpl::nofile())
  	    << ccmpl::stop;
  return 0;
}



