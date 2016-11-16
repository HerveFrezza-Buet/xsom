#pragma once 

#include <vector>
#include <algorithm>

namespace xsom {


  class Entity {
    
  public:
    
    /**
     * Freezes all updates whatever the lock values.
     */
    bool freeze;
    bool update_locked;
    bool learn_locked;
    Entity() : freeze(false), update_locked(false), learn_locked(false) {}
    Entity(const Entity& cp) : freeze(cp.freeze), update_locked(cp.update_locked), learn_locked(cp.learn_locked) {}
    Entity& operator=(const Entity& cp) {
      if(this != &cp) {
	freeze        = cp.freeze;
	update_locked = cp.update_locked;
	learn_locked  = cp.learn_locked;
      }
      return *this;
    }
    virtual ~Entity() {}
    virtual void on_learn()=0;
    virtual void on_update()=0;
    void learn() {
      if(!freeze && !learn_locked) on_learn();
    }

    void update() {
      if(!freeze && !update_locked) on_update();
    }
    
  };


  class Container : public xsom::Entity { 
  private:

    std::vector<xsom::Entity*> elems;
    bool shuffup;

  public:
    
    Container() : elems(), shuffup(false) {}
    Container(bool shuffled_update) : xsom::Entity(), elems(), shuffup(shuffled_update)  {}
    Container(const Container& cp) : xsom::Entity(cp), elems(cp.elems), shuffup(cp.shuffup) {}
    Container& operator=(const Container& cp) {
      if(&cp != this) {
	this->xsom::Entity::operator=(cp);
	elems   = cp.elems;
	shuffup = cp.shuffup;
      }
      return *this;
    }
    virtual ~Container() {}

    Container& operator+=(xsom::Entity& entity) {
      elems.push_back(&entity);
      return *this;
    }

    virtual void on_learn() {
      for(auto elem_ptr : elems) elem_ptr->on_learn();
    }

    virtual void on_update() {
      if(shuffup)
	std::random_shuffle(elems.begin(),elems.end());
      for(auto elem_ptr : elems) elem_ptr->update();
    }
  };
}
