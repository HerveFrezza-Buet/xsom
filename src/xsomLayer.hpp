#pragma once

#include <functional>

#include <xsomEntity.hpp>


namespace xsom {


  template<typename WEIGHT,
	   typename POSITION,
	   typename INPUT>
  class Layer : public xsom::Entity {
  public:

    /**
     * This is the layer's learning rate.
     */
    double alpha;

  private:

    INPUT* input;
      
    std::function<double (const POSITION&,const WEIGHT&, const INPUT&)>    match;
    std::function<double (const POSITION&)>                                neighbour;
    std::function<WEIGHT (const POSITION&)>                                get_weight;
    std::function<void (std::function<double (const POSITION&)>)>          learn_activity;
    std::function<void (std::function<WEIGHT (const POSITION&)>)>          learn_weight;

  public:

    template<typename fctMATCH,
	     typename fctNEIGHBOUR,
	     typename fctGET_WEIGHT,
	     typename fctLEARN_ACTIVITY,
	     typename fctLEARN_WEIGHT>
    Layer(double learning_rate,
	  INPUT& external_input,
	  const fctMATCH&          fct_match,
	  const fctNEIGHBOUR&      fct_neighbour,
	  const fctGET_WEIGHT&     fct_get_weight,
	  const fctLEARN_ACTIVITY& fct_learn_activity,
	  const fctLEARN_WEIGHT&   fct_learn_weight) 
      : xsom::Entity(), alpha(learning_rate),
	input(&external_input),
	match(fct_match),
	neighbour(fct_neighbour),
	get_weight(fct_get_weight),
	learn_activity(fct_learn_activity),
	learn_weight(fct_learn_weight) {}
    
    Layer() : xsom::Entity() {}


    Layer(const Layer& cp)
      : xsom::Entity(cp), 
	alpha(cp.alpha),
	input(cp.input),
	match(cp.match),
	neighbour(cp.neighbour),
	get_weight(cp.get_weight),
	learn_activity(cp.learn_activity),
	learn_weight(cp.learn_weight) {}

    Layer& operator=(const Layer& cp) {
      if(this != &cp) {
	this->xsom::Entity::operator=(cp);
	alpha = cp.alpha;
	input = cp.input;
	match = cp.match;
	neighbour = cp.neighbour;
	get_weight = cp.get_weight;
	learn_activity = cp.learn_activity;
	learn_weight = cp.learn_weight;
      }
      return *this;
    }
    
    
    virtual ~Layer() {}

    virtual void on_update() {
      auto& _input      = *input;
      auto& _match      = match;
      auto& _get_weight = get_weight;
      learn_activity([_input,_match,_get_weight](const POSITION& p) -> double {return _match(p,_get_weight(p),_input);});
    }

    virtual void on_learn() {
      double _alpha      = this->alpha;
      auto&  _input      = *input;
      auto&  _neighbour  = neighbour;
      auto&  _get_weight = get_weight;
      learn_weight([_alpha,_input,_neighbour,_get_weight](const POSITION& p) 
		   -> WEIGHT {
		     WEIGHT w = _get_weight(p);
		     double coef = _neighbour(p);
		     if(coef == double(0))
		       return w;
		     else
		       return w+(_input-w)*(_alpha*coef);
		   });
    }
  };

  
  template<typename WEIGHT,
	   typename POSITION,
	   typename INPUT,
	   typename fctMATCH,
	   typename fctNEIGHBOUR,
	   typename fctGET_WEIGHT,
	   typename fctLEARN_ACTIVITY,
	   typename fctLEARN_WEIGHT>
  Layer<WEIGHT,
	POSITION,
	INPUT> layer(double learning_rate,
		     INPUT& external_input,
		     const fctMATCH&          fct_match,
		     const fctNEIGHBOUR&      fct_neighbour,
		     const fctGET_WEIGHT&     fct_get_weight,
		     const fctLEARN_ACTIVITY& fct_learn_activity,
		     const fctLEARN_WEIGHT&   fct_learn_weight) {
    return Layer<WEIGHT,
		 POSITION,
		 INPUT>(learning_rate,
			external_input,
			fct_match,
			fct_neighbour,
			fct_get_weight,
			fct_learn_activity,
			fct_learn_weight);
  }
}
