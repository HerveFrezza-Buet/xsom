#include <iostream>
#include <xsomInputs.hpp>

int main(int argc, char* argv[]) {

  {
    std::cout << "Sampling the MackeyGlass serie" << std::endl;
    xsom::input::continuous::MackeyGlass input;
    while(input.time() < 200) {
      std::cout << input.input() << "(" << input.time() << ") ";
      input.shift();
    }
    std::cout << std::endl;
  }

  {
    std::cout << "Sampling the binary automata" << std::endl;
    xsom::input::discrete::BinaryAutomata input(0.2, 0.5);
    int n0 = 0;
    int n1 = 0;
    while(input.time() < 10000) {
      std::cout << input.input() << " ";
      if(input.input())
	n1 += 1;
      else
	n0 += 1;
      input.shift();
    }
    std::cout << std::endl;
    std::cout << "0 (" << double(n0)/double(n0+n1) << ")" << std::endl
	      << "1 (" << double(n1)/double(n0+n1) << ")" << std::endl;
    std::cout << std::endl;
  }
}
