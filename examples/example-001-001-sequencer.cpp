#include <iostream>
#include <functional>

#include <xsom.hpp>
#include <ccmpl.hpp>
using namespace std::placeholders;

/* This example shows how to use a sequencer in order to control executions. */


// This draws a gabor function depending on time.
void fill_data(std::vector<ccmpl::Point>& curve, double& time) {
  curve.clear();
  for(double x = -3; x <= 3 ; x += .05) 
    curve.push_back({x,std::cos(10*(x+.01*time))*std::exp(-x*x)});
}


std::string flags() {return "#";}


#define VIEW_PREFIX "viewer-001-001"
int main(int argc, char* argv[]) {
  
  ccmpl::Main m(argc,argv,VIEW_PREFIX);
  
  double current_time = 0;
  
  // Plot
  
  auto display  = ccmpl::layout(5.0, 5.0, {"#"});
  display()     = {-5, 5, -1, 1};   
  display()    += ccmpl::line("'b-'", std::bind(fill_data, _1, std::ref(current_time))); 

  m.generate(display,true);

  // Computation : let us build a fake architecture

  auto archi  = xsom::setup::network();
  auto map    = xsom::setup::debug("fake map");
  archi      += map;

  // Let us build a sequencer for synchronizing the computation
  auto seq    = xsom::setup::sequencer(archi,display);

  // See the doxygen documentation of xsom::setup::Sequencer for an
  // exhaustive list of sequencer functionalities.

  /* */ seq.__def("tick");
  /* */   seq.__step([&current_time](){++current_time;});
  /* */ seq.__fed();
  
  /* */ seq.__def("update step");
  /* */   seq.__update();
  /* */   seq.__call("tick");
  /* */   seq.__plot_png(flags, "update-frame");
  /* */ seq.__fed();
  
  /* */ seq.__def("learn step");
  /* */   seq.__update_and_learn();
  /* */   seq.__call("tick");
  /* */   seq.__plot(flags, "learn-frame", "update-frame");
  /* */ seq.__fed();
  
  /* */ seq.__for(10);
  /* */   seq.__for(9);
  /* */     seq.__call("update step");
  /* */   seq.__rof();
  /* */   seq.__call("learn step");
  /* */ seq.__rof();
  
  // Now we can run the simulation.
  seq.run();
  
}
