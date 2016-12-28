#include <iostream>

#include <xsom.hpp>
#include <ccmpl.hpp>

/* This example shows how to use a sequencer in order to control executions. */



#define VIEW_FILE "viewer-001-001.py"
int main(int argc, char* argv[]) {
  if(argc != 2) {
    std::cout << "Usage : " << std::endl
	      << argv[0] << " generate" << std::endl
	      << argv[0] << " run | ./" << VIEW_FILE << std::endl;
    return 0;
  }

  bool generate_mode = std::string(argv[1])=="generate";
  
  // Plot
  
  auto display = ccmpl::layout(5.0, 5.0, {"#"});
  
  if(generate_mode) {
    display.make_python(VIEW_FILE,true);
    return 0;
  }

  // Computation : let us build a fake architecture

  auto archi  = xsom::setup::network();
  auto map    = xsom::setup::debug("fake map");
  archi      += map;

  // Let us build a sequencer for synchronizing the computation
  auto seq    = xsom::setup::sequencer(archi,display);


  unsigned int i = 0;
  seq.__for(20);
  seq.__if([&i](){return (i/2)%2 == 0;});
  /* */ seq.__if([&i](){return i%2 == 0;});
  /* */ seq.__else();
  /*   */ seq.__print("BBBB");
  /* */ seq.__fi();
  seq.__else();
  /* */ seq.__if([&i](){return i%2 == 0;});
  /*   */ seq.__print("CCCC");
  /* */ seq.__else();
  /* */ seq.__fi();
  seq.__fi();
  seq.__step([&i](){++i;});
  seq.__rof();

  
  

  // Now we can run the simulation.
  seq.run();
  
}
