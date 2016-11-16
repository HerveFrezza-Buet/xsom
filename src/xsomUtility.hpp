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
}
