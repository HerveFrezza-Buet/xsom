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
	std::bernoulli_distribution d_0; // probability P(x_{t+1} = 0 / x_t = 0)
	std::bernoulli_distribution d_1; // probability P(x_{t+1} = 1 / x_t = 1)
	
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
	  
	  if(state) { // x_t = 1
	    bool stay = d_1(gen);
	    if(! stay)
	      state = false;
	  }
	  else {     // x_t = 0
	    bool stay = d_0(gen);
	    if(! stay)
	      state = true;
	  }
	}	
      };
    }
    
    namespace continuous {
      class MackeyGlass {
	// http://www.scholarpedia.org/article/Mackey-Glass_equation
      private:
	double gamma, beta, tau, n;
	double t, dt;
	double x, x_tau;
	std::deque<double> history;
      public:
	MackeyGlass(): MackeyGlass(1., 2., 2., 9.65, 0.1, 0.5) {
	}
	MackeyGlass(double gamma, double beta, double tau, double n, double dt, double x_tau):
	  gamma(gamma), beta(beta), tau(tau), n(n), t(0), dt(dt), x(x_tau), x_tau(x_tau){
	  history.resize(int(tau/dt));
	  std::fill(history.begin(), history.end(), x_tau);
	}
	double time() const {
	  return t;
	}
      
	double input() const {
	  return x;
	}
	void shift() {
	  x_tau = history.front();
	  history.pop_front();
	  x = x + dt * (beta * x_tau / (1.0 + pow(x_tau, 10.0)) - gamma * x);
	  history.push_back(x);
	  t += dt;
	
	}
      };
    }
  }
}
