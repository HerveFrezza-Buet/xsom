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

  if(uv.y <= -1.5 - uv.x)
	  return 0.8;
  if((uv.x-1)*(uv.x-1) + (uv.y-1)*(uv.y-1) <= 0.5*0.5)
	  return 0.8;
  
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

#define VIEW_PREFIX "viewer-002-002"
int main(int argc, char* argv[]) {
  
  ccmpl::Main m(argc,argv,VIEW_PREFIX);

  auto display      = ccmpl::layout(m.hostname, m.port,
				    10,10, {"##", "##"});
  std::string flags = "";

  // Let us compute a convolution function.
  
  auto f_mapping = xsom::tab2d::mapping({-1, -1}, {1, 1}, {NB_U, NB_V});
  auto tabular_f = xsom::tab2d::fft::convolution(f_mapping, SIGMA, xsom::tab::fft::KernelType::Gaussian, xsom::tab::fft::PaddingType::Constant);

  tabular_f.learn(f);
  tabular_f.convolve(); // This does the convolution.

  auto conv_bmu       = tabular_f.bmu();     // This call invokes  tabular_f.convolve() internally.
  auto conv_bmu_val   = tabular_f(conv_bmu);
  auto noconv_bmu     = tabular_f.bmu_noconv();   
  auto noconv_bmu_val = tabular_f.get_noconv(noconv_bmu);
  
  
  // The convolution offers the same functions as usual tabulars.
  
  display()         = ccmpl::view2d({-1, 1}, {-1, 1}, ccmpl::aspect::equal, ccmpl::span::placeholder);
  display()         = ccmpl::show_tics(false, false);
  display().title   = "convolution (constant pad)";
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
  display()         = ccmpl::view2d({-1, 1}, {-1, 1}, ccmpl::aspect::equal, ccmpl::span::placeholder);
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
  auto tabular_g = xsom::tab2d::fft::convolution(g_mapping, SIGMA, xsom::tab::fft::KernelType::Gaussian, xsom::tab::fft::PaddingType::Zero,
						 [](const xsom::Point2D<double> uv) {return uv*uv <= .6;});
  tabular_g.clear(0);                             // Since we use a mask, we clear for initialization.
  tabular_g.learn(f);
  auto g_conv_bmu       = tabular_g.bmu();        // This call invokes  tabular_g.convolve() internally.
  auto g_conv_bmu_val   = tabular_g(g_conv_bmu);
  auto g_noconv_bmu     = tabular_g.bmu_noconv();   
  
  
  display++;
  display()         = ccmpl::view2d({-1, 1}, {-1, 1}, ccmpl::aspect::equal, ccmpl::span::placeholder);
  display()         = ccmpl::show_tics(false, false);
  display().title   = "convolution (masked, zero pad)";
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
  display()         = ccmpl::view2d({-1, 1}, {-1, 1}, ccmpl::aspect::equal, ccmpl::span::placeholder);
  display()         = ccmpl::show_tics(false, false);
  display().title   = "no convolution (masked)";
  display().xtitle  = "u";
  display().ytitle  = "v";
  display()         = "equal";
  display()        +=  ccmpl::image("cmap='binary', interpolation='bilinear', clim=(0,1)",
				    std::bind(&xsom::tab2d::fft::Convolution::fill_image_gray_noconv, std::ref(tabular_g), _1, _2, _3, _4, _5));
  flags += '#';



  // the ccmpl::Main object handles generation here
  m.generate(display, false); // true means "use GUI"
  
  display(flags,"example-002-002.pdf", ccmpl::nofile());
  !display;
  
  return 0;
}



