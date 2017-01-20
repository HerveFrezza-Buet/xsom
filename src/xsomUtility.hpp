#pragma once

#include <cmath>
#include <string>
#include <sstream>
#include <iomanip>

namespace xsom {
  template<typename VALUE,
	   typename fctSQUARE_DIST>
  double gaussian(const VALUE& val1,
	     const VALUE& val2,
	     double sigma2,
	     const fctSQUARE_DIST& dist2) {
    return ::exp(-dist2(val1,val2)/(2*sigma2));
  }

  template<typename VALUE,
	   typename fctSQUARE_DIST>
  double linear(const VALUE& val1,
		const VALUE& val2,
		double r2,
		const fctSQUARE_DIST& dist2) {
    double d2 = dist2(val1,val2);
    if(d2 > r2)
      return 0;
    return 1-std::sqrt(d2/r2);
  }

  template<typename VALUE,
	   typename fctSQUARE_DIST>
  double parabolic(const VALUE& val1,
		   const VALUE& val2,
		   double r2,
		   const fctSQUARE_DIST& dist2) {
    double d2 = dist2(val1,val2);
    if(d2 > r2)
      return 0;
    return 1-d2/r2;
  }


  namespace random {

    /**
     * @return a random value in [min,max[
     */
    inline double uniform(double min,double max) {
      return min + (max-min)*(std::rand()/(RAND_MAX+1.0));
    }

    /**
     * @return a random integer in [0,max[
     */
    template<typename VALUE>
    inline VALUE uniform(VALUE max) {
      return (VALUE)(max*(std::rand()/(RAND_MAX+1.0)));
    }


    /**
     * @param p in [0,1]
     * @return true with the probability proba.
     */
    inline bool proba(double p) {
      return uniform(0,1)<p;
    }
  }
}
