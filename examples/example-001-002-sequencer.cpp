#include <iostream>
#include <functional>

#include <xsom.hpp>
#include <ccmpl.hpp>
using namespace std::placeholders;

/* This example shows how to use an interactive sequencer in order to control executions. */


// This draws a gabor function depending on time.
void fill_data(std::vector<ccmpl::Point>& curve, double& time) {
  curve.clear();
  for(double x = -3; x <= 3 ; x += .05) 
    curve.push_back({x,std::cos(10*(x+.01*time))*std::exp(-x*x)});
}

std::string flags() {return "#";}


#define VIEW_PREFIX "viewer-001-002"
#define PIPE_NAME   "/tmp/ccmpl"
int main(int argc, char* argv[]) {

  if(argc == 1) {
    std::cout << std::endl
	      << std::endl
	      << "Usage :" << std::endl
	      << "  mkfifo " << PIPE_NAME << std::endl
	      << "  " << argv[0] << " display" << std::endl
	      << "  python3 ./" << VIEW_PREFIX << ".py " << PIPE_NAME << std::endl
	      << "  " << argv[0] << " run    (on another terminal)" << std::endl
	      << std::endl
	      << std::endl
	      << std::endl;
    ::exit(0);
  }
  
  ccmpl::Main m(argc,argv,VIEW_PREFIX);
  
  double current_time = 0;
  
  // Plot
  
  auto display  = ccmpl::layout(5.0, 5.0, {"#"},
				ccmpl::RGB(1., 1., 1.));
  display()     = ccmpl::view2d({-5, 5}, {-1, 1}, ccmpl::aspect::fit, ccmpl::span::placeholder); 
  display()    += ccmpl::line("'b-'", std::bind(fill_data, _1, std::ref(current_time))); 

  m.generate(display,true);


  // Let us build a sequencer for synchronizing the computation
  auto seq    = xsom::setup::sequencer(display);

  // Let us turn the sequence into a keyboard interaction mode. We add a custom menue item.
  seq.add_menu_item('x',                    // The key to be pressed.
		    "the 'x' key",          // The key description.
		    "Prints 'Hello world'", // The key description.
		    [&seq]() {              // What that key does.
		      seq.msg_info("Hello World");
		    });
  seq.interactive(true, "/tmp/ccmpl"); // call this after having added menus.
  
  // See the doxygen documentation of xsom::setup::Sequencer for an
  // exhaustive list of sequencer functionalities.

  /* */ seq.__def("tick");
  /* */   seq.__();
  /* */   seq.__step([&current_time](){++current_time;});
  /* */ seq.__fed();
  
  /* */ seq.__def("update step");
  /* */   seq.__update();
  /* */   seq.__plot_png(flags, "update-frame");
  /* */   seq.__call("tick");
  /* */ seq.__fed();
  
  /* */ seq.__def("learn step");
  /* */   seq.__update_and_learn();
  /* */   seq.__plot(flags, "learn-frame", "update-frame");
  /* */   seq.__call("tick");
  /* */ seq.__fed();
  
  /* */ seq.__loop();
  /* */   seq.__for(9);
  /* */     seq.__call("update step");
  /* */   seq.__rof();
  /* */   seq.__call("learn step");
  /* */ seq.__pool();
  
  // Now we can run the simulation.
  seq.run();
  
}



