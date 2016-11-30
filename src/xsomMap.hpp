#pragma once

#include <xsomEntity.hpp>
#include <xsomPoint.hpp>

#include <ccmpl.hpp>

namespace xsom {

  template<typename POS>
  class Map : public xsom::Container {
  private:

    POS* out;
    std::function<void (std::function<double (const POS&)>)> learn_activity;
    std::function<double (const POS&)>                       merge_activity;
    std::function<POS ()>                                    best_position;
    std::function<void (POS&, const POS&)>                   learn_output;


  public:

    Map() : xsom::Container() {}

    template<typename fctLEARN_ACTIVITY,
	     typename fctMERGE_ACTIVITY,
	     typename fctBEST_POSITION,
	     typename fctLEARN_OUTPUT>
    Map(POS& output,
	const fctLEARN_ACTIVITY& fct_learn_activity,
	const fctMERGE_ACTIVITY& fct_merge_activity,
	const fctBEST_POSITION&  fct_best_position,
	const fctLEARN_OUTPUT&   fct_learn_output) 
      : xsom::Container(false),
	out(&output),
	learn_activity(fct_learn_activity),
	merge_activity(fct_merge_activity),
	best_position(fct_best_position),
	learn_output(fct_learn_output) {}
    
    Map(const Map& cp) 
      : xsom::Container(cp),
	out(cp.out),
	learn_activity(cp.learn_activity),
	merge_activity(cp.merge_activity),
	best_position(cp.best_position),
	learn_output(cp.learn_output) {}

    Map& operator=(const Map& cp) {
      if(&cp != this) {
	this->xsom::Container::operator=(cp);
	out = cp.out;
	learn_activity = cp.learn_activity;
	merge_activity = cp.merge_activity;
	best_position = cp.best_position;
	learn_output = cp.learn_output;
      }

      return *this;
    }
    
    virtual ~Map() {}


    virtual void on_update() {
      this->xsom::Container::on_update();
      learn_activity(merge_activity);
      learn_output(*out,best_position());
    }
  };

  template<typename POS,
	   typename fctLEARN_ACTIVITY,
	   typename fctMERGE_ACTIVITY,
	   typename fctBEST_POSITION,
	   typename fctLEARN_OUTPUT>
  Map<POS> map(POS& output,
	       const fctLEARN_ACTIVITY& fct_learn_activity,
	       const fctMERGE_ACTIVITY& fct_merge_activity,
	       const fctBEST_POSITION&  fct_best_position,
	       const fctLEARN_OUTPUT&   fct_learn_output) {
    return Map<POS>(output,
		    fct_learn_activity,
		    fct_merge_activity,
		    fct_best_position,
		    fct_learn_output);
  }
}

