#include <utility>
#include <vector>
#include <functional>
#include <algorithm>
#include <cmath>
#include <string>

#include <xsom.hpp>
#include <ccmpl.hpp>

using namespace std::placeholders;

// This shows how to use convolution in xsom. Read first
// example-000-004-tabular.cc, since the convolution data in xsom is
// derived from xsom::tab1d::Table<double> and
// xsom::tab2d::Table<double>

// Let us define a function

double f(const xsom::Point2D<double>& uv) {
  if(uv*uv < .01)
    return 1;
    
  if(uv*uv < .25) {
    if(uv.x*uv.y > 0)
      return .8;
    else
      return 0;
  }
  
  if(uv.x*uv.y > 0)
    return 0;
  
  return .8;
}

// Let us define a convolution layer. It consists of a tabular object
// (xsom::tab*d::Table<double>), that deals with two fields of scalar
// values. The second is the convolution of the first one.

#define NB_U 200
#define NB_V  50

#define SIGMA (NB_U/30)

#define VIEW_FILE "viewer-002-002.py"
int main(int argc, char* argv[]) {
  if(argc < 2) {
    std::cout << "Usage : " << std::endl
  	      << argv[0] << " generate" << std::endl
  	      << argv[0] << " run | ./" << VIEW_FILE << std::endl;
    return 0;
  }

  bool generate_mode = std::string(argv[1])=="generate";

  auto display      = ccmpl::layout(10,10, {"##", "##"});
  std::string flags = "";

  // Let us compute a convolution function.
  
  auto f_mapping = xsom::tab2d::mapping({-1, -1}, {1, 1}, {NB_U, NB_V});
  auto tabular_f = xsom::tab2d::fft::convolution(f_mapping, SIGMA, xsom::tab::fft::KernelType::Gaussian);

  tabular_f.learn(f);
  tabular_f.convolve(); // This does the convolution.

  auto conv_bmu       = tabular_f.bmu();     // This call invokes  tabular_f.convolve() internally.
  auto conv_bmu_val   = tabular_f(conv_bmu);
  auto noconv_bmu     = tabular_f.bmu_noconv();   
  auto noconv_bmu_val = tabular_f.get_noconv(noconv_bmu);
  
  
  // The convolution offers the same functions as usual tabulars.
  
  display()         = {-1, 1, -1, 1};
  display()         = ccmpl::show_tics(false, false);
  display().title   = "convolution";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        += ccmpl::image(std::string("cmap='jet', interpolation='bilinear', clim=(0,") + std::to_string(conv_bmu_val) + "), zorder=0",
				   std::bind(&xsom::tab2d::fft::Convolution::fill_image_gray, std::ref(tabular_f), _1, _2, _3, _4, _5));
  flags += "#";
  display()        += ccmpl::dot("c='g',lw=1,s=50,zorder=1", [conv_bmu](ccmpl::Point& dot) {dot = conv_bmu;}); 
  flags += "#";
  display()        += ccmpl::dot("c='w',lw=1,s=50,zorder=1", [noconv_bmu](ccmpl::Point& dot) {dot = noconv_bmu;}); 
  flags += "#";

  display++;
  display()         = {-1, 1, -1, 1};
  display()         = ccmpl::show_tics(false, false);
  display().title   = "no convolution";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        +=  ccmpl::image("cmap='binary', interpolation='bilinear', clim=(0,1)",
				    std::bind(&xsom::tab2d::fft::Convolution::fill_image_gray_noconv, std::ref(tabular_f), _1, _2, _3, _4, _5));
  flags += '#';

  // Let us now do the same with a position mask.

  auto g_mapping = xsom::tab2d::mapping({-1, -1}, {1, 1}, {NB_U, NB_U});  // NB_U twice here...
  auto tabular_g = xsom::tab2d::fft::convolution(g_mapping, SIGMA, xsom::tab::fft::KernelType::Gaussian,
						 [](const xsom::Point2D<double> uv) {return uv*uv <= .6;});
  tabular_g.clear(0);                             // Since we use a mask, we clear for initialization.
  tabular_g.learn(f);
  auto g_conv_bmu       = tabular_g.bmu();        // This call invokes  tabular_g.convolve() internally.
  auto g_conv_bmu_val   = tabular_g(g_conv_bmu);
  auto g_noconv_bmu     = tabular_g.bmu_noconv();   
  
  
  display++;
  display()         = {-1, 1, -1, 1};
  display()         = ccmpl::show_tics(false, false);
  display().title   = "convolution (masked)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        += ccmpl::image(std::string("cmap='jet', interpolation='bilinear', clim=(0,") + std::to_string(g_conv_bmu_val) + "), zorder=0",
				   std::bind(&xsom::tab2d::fft::Convolution::fill_image_gray, std::ref(tabular_g), _1, _2, _3, _4, _5));
  flags += "#";
  display()        += ccmpl::dot("c='g',lw=1,s=50,zorder=1", [g_conv_bmu](ccmpl::Point& dot) {dot = g_conv_bmu;}); 
  flags += "#";
  display()        += ccmpl::dot("c='w',lw=1,s=50,zorder=1", [g_noconv_bmu](ccmpl::Point& dot) {dot = g_noconv_bmu;}); 
  flags += "#";

  display++;
  display()         = {-1, 1, -1, 1};
  display()         = ccmpl::show_tics(false, false);
  display().title   = "no convolution (masked)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        +=  ccmpl::image("cmap='binary', interpolation='bilinear', clim=(0,1)",
				    std::bind(&xsom::tab2d::fft::Convolution::fill_image_gray_noconv, std::ref(tabular_g), _1, _2, _3, _4, _5));
  flags += '#';


  
  if(generate_mode) {
    display.make_python(VIEW_FILE,false);
    return 0;
  }
  
  std::cout << display(flags,"example-002-002.pdf", ccmpl::nofile())
  	    << ccmpl::stop;
  return 0;
}



