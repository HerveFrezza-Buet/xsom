#include <iostream>
#include <xsomInputs.hpp>
#include <fstream>

int main(int argc, char* argv[]) {

  {
    std::cout << "Sampling the MackeyGlass serie" << std::endl;
    xsom::input::continuous::MackeyGlass input;
    std::ofstream outfile("mackey.data");
    while(input.time() < 50) {
      std::cout << input.input() << "(" << input.time() << ") ";
      outfile << input.time() << "\t" << input.input() << std::endl;
      input.shift();
      
    }
    outfile.close();
    std::cout << std::endl;
    std::cout << "Samples saved in mackey.data" << std::endl;
  }

  {
    // As in (Voegtlin, 2002)
    std::cout << "Sampling the binary automata" << std::endl;
    xsom::input::discrete::BinaryAutomata input(0.7, 0.6);
    while(input.time() < 10) {
      std::cout << input.input() << " ";
      input.shift();
    }
    std::cout << std::endl;
  }
}
