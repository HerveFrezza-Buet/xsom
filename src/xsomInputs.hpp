#pragma once

#include <deque>
#include <cmath>
#include <random>

namespace xsom {
  namespace input {
    namespace discrete {
      class BinaryAutomata {
	bool state;
	double t;
	std::random_device rd;
	std::mt19937 gen;
	std::bernoulli_distribution d_0;
	std::bernoulli_distribution d_1;
	
      public:
	BinaryAutomata(double p00, double p11):
	  state(false), t(0.0),
	  gen(rd()), d_0(p00), d_1(p11) {
	  
	}

	double time() const {
	  return t;
	}
	bool input() const {
	  return state;
	}
	void shift() {
	  t+=1;
	  
	  if(state) {
	    bool stay = d_1(gen);
	    if(! stay)
	      state = false;
	  }
	  else {
	    bool stay = d_0(gen);
	    if(! stay)
	      state = true;
	  }
	}	
      };
    }
    
    namespace continuous {
    class MackeyGlass {
    private:
      double a,b,d;
      double t, dt;
      double x, x_d;
      std::deque<double> history;
    public:
      MackeyGlass(): MackeyGlass(0.2, -0.1, 17.0, 0.1, 0.5) {
      }
      MackeyGlass(double a, double b, double d, double dt, double x_d):
	a(a), b(b), d(d), t(0), dt(dt), x(x_d), x_d(x_d){
	history.resize(int(d/dt));
	std::fill(history.begin(), history.end(), x_d);
      }
      double time() const {
	return t;
      }
      
      double input() const {
	return x;
      }
      void shift() {
	x_d = history.front();
	history.pop_front();
	x = x + dt * (b * x + a * x_d / (1.0 + pow(x_d, 10.0)));
	history.push_back(x);
	t += dt;
	
      }
    };
    }
  }
}
