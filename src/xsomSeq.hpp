#pragma once

#include <xsomNetwork.hpp>
#include <ccmpl.hpp>
#include <map>
#include <iterator>
#include <string>
#include <iostream>
#include <sstream>
#include <iomanip>
#include <list>
#include <stack>
#include <memory>
#include <functional>

namespace xsom {
  namespace setup {
    class Sequencer;
  }
  
  namespace instr {

    class Instruction;
   
    
    using Instr = std::shared_ptr<Instruction>;

    struct Done {
      Done() {}
    };

    class Instruction {
    protected:
      bool has_next;
      Instruction() : has_next(true) {}
      virtual void execute() = 0;

    public:

      virtual Instr deep_copy() = 0;
      void next() {
	if(has_next)
	  execute();
	else
	  throw Done();
      }
    };

    class Step : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      std::function<void ()> f;
      
      template<typename Fun>
      Step(const Fun& f) : Instruction(), f(f) {}

    protected:
      
      virtual void execute() {
	f();
	has_next = false;
      }

    public:
      
      virtual Instr deep_copy() {
	return Instr(new Step(f));
      }
      
    };

    class Call : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      Instr instr;
      Call(Instr instr) : Instruction(), instr(instr->deep_copy()) {}

    protected:
      
      virtual void execute() {
	try {
	  instr->next();
	}
	catch(Done) {
	  has_next = false;
	}
      }

    public:
      
      virtual Instr deep_copy() {
	return Instr(new Call(instr)); // deep copy is performed by Call(...).
      }
    };
    
    class Seq : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      std::list<Instr> seq;
      std::list<Instr>::iterator current;
      
      Seq(const std::list<Instr>& seq) : Instruction(), seq(seq) {
	current = this->seq.begin();
      }

    protected:

      virtual void execute() {
	while(current != seq.end()) {
	  try {
	    (*current)->next();
	    break;
	  }
	  catch(Done d) {
	    ++current;
	  }
	}
	has_next = current != seq.end();
      }
      
    public:
      
      virtual Instr deep_copy() {
	std::list<Instr> _seq;
	auto out = std::back_inserter(_seq);
	for(auto i : seq) *(out++) = i->deep_copy();
	return Instr(new Seq(_seq));
      }
    };

    class Loop : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      Instr body;
      
      Loop(Instr body) : Instruction(), body(body) {}

    protected:
      
      virtual void execute() {
	try {
	  body->next();
	}
	catch(Done e) {
	  body = body->deep_copy();
	  try {
	    body->next();
	  }
	  catch(Done e) {
	    std::cerr << "Sequencer error : loop has an invalid instruction" << std::endl;
	  }
	}
      }

    public:
      
      virtual Instr deep_copy() {
	return Instr(new Loop(body->deep_copy()));
      }
      
    };


    class For : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      Instr body;
      unsigned int nb,i;
      
      For(Instr body, unsigned int nb_times) : Instruction(), body(body), nb(nb_times), i(0) {}

    protected:
      
      virtual void execute() {
	if(i == nb)
	  has_next = false;
	else {
	  try {
	    body->next();
	  }
	  catch(Done e) {
	    ++i;
	    body = body->deep_copy();
	    execute();
	  }
	}
      }

    public:
      
      virtual Instr deep_copy() {
	return Instr(new For(body->deep_copy(), nb));
      }
      
    };


    


    class While : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      Instr body;
      std::function<bool ()> test;
      bool first;

      template<typename TEST>
      While(Instr body, const TEST& test) : Instruction(), body(body), test(test), first(true) {}

    protected:
      
      virtual void execute() {
	if(first) {
	  has_next = test();
	  first = false;
	}

	if(has_next) {
	  try {
	    body->next();
	  }
	  catch(Done e) {
	    has_next = test();
	    if(has_next) {
	      body = body->deep_copy();
	      execute();
	    }
	  }
	}
      }

    public:
      
      virtual Instr deep_copy() {
	return Instr(new While(body->deep_copy(), test));
      }
      
    };

    
  }

  namespace setup {

    class Sequencer {
    private:
      xsom::Container* archi;
      ccmpl::chart::Layout* display;

      std::map<std::string,unsigned int>         frame_pdf;
      std::map<std::string,unsigned int>         frame_png;
      std::map<std::string,xsom::instr::Instr>   macros;
      std::stack<std::list<xsom::instr::Instr>>  context;
      std::stack<std::string>                    macro_names;
      std::stack<unsigned int>                   for_times;
      std::stack<std::function<bool ()>>         tests;
      

      std::string next(std::map<std::string,unsigned int>& frame,
		       const std::string& name,
		       const std::string& suffix) {
	auto kv = frame.find(name);
	unsigned int id = 0;
	if(kv == frame.end())
	  frame[name] = 1;
	else
	  id = (kv->second)++;
	std::ostringstream ostr;
	ostr << name << '-' << std::setw(6) << std::setfill('0') << id
	     << '.' << suffix;
	return ostr.str();
      }

      std::string next_pdf(const std::string& name) {
	return next(frame_pdf, name, "pdf");
      }

      std::string next_png(const std::string& name) {
	return next(frame_png, name, "png");
      }

      
      Sequencer(xsom::Container* archi, ccmpl::chart::Layout* display)
	: archi(archi), display(display) {
	context.push(std::list<xsom::instr::Instr>());
      }

    public:

      Sequencer(const Sequencer&) = default;
      Sequencer& operator=(const Sequencer&) = default;

      Sequencer(xsom::Container& archi, ccmpl::chart::Layout& display)
	: Sequencer(&archi, &display) {}

      Sequencer(xsom::Container& archi)
	: Sequencer(&archi, nullptr) {}
    
      Sequencer(ccmpl::chart::Layout& display)
	: Sequencer(nullptr, &display) {}

      Sequencer()
	: Sequencer(nullptr,nullptr) {}

      /**
       * Add a step that calls a macro.
       */
      void __call(const std::string& macro_name) {
	auto kv = macros.find(macro_name);
	if(kv == macros.end())
	  std::cerr << "Sequencer error : macro \"" << macro_name
		    << "\" undefined." << std::endl;
	else
	  context.top().push_back(xsom::instr::Instr(new xsom::instr::Call(kv->second)));
      }

      /**
       * defines a macro.
       */
      void __def(const std::string& macro_name) {
	macro_names.push(macro_name);
	context.push(std::list<xsom::instr::Instr>());
      }
      
      /**
       * ends a macro definition
       */
      void __fed() {
	macros[macro_names.top()] = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	context.pop();
	macro_names.pop();
      }

      /**
       * Add a step that calls f. f is "void f()".
       */
      template<typename Fun>
      void __step(const Fun& f) {
	context.top().push_back(xsom::instr::Instr(new xsom::instr::Step(f)));
      }

      /**
       * Add a step which is a sequence... until end is called.
       */
      void __begin() {
	context.push(std::list<xsom::instr::Instr>());
      }

      
      /**
       * Close the sequence opened by begin.
       */
      void __end() {
	auto seq = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	context.pop();
	context.top().push_back(seq);
      }

      /**
       * Add a step that prints a message on cerr
       */
      void __print(const std::string& message) {
	__step([message](){std::cerr << "Sequencer info : " << message << std::endl;});
      }

      /**
       * Define an infinite loop. end by pool
       */
      void __loop() {
	context.push(std::list<xsom::instr::Instr>());
      }
      
      /**
       * close the _loop call.
       */
      void __pool() {
	auto seq = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	context.pop();
	context.top().push_back(xsom::instr::Instr(new xsom::instr::Loop(seq)));
      }
      
      /**
       * Define an nb_times repetition. end by rof
       */
      void __for(unsigned int nb_times) {
	context.push(std::list<xsom::instr::Instr>());
	for_times.push(nb_times);
      }
      
      /**
       * close the _for call.
       */
      void __rof() {
	auto seq = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	context.pop();
	context.top().push_back(xsom::instr::Instr(new xsom::instr::For(seq,for_times.top())));
	for_times.pop();
      }
      
      /**
       * Define an while. test is "bool test()". end by elihw
       */
      template<typename TEST>
      void __while(const TEST& test) {
	context.push(std::list<xsom::instr::Instr>());
	tests.push(test);
      }
      
      /**
       * close the _while call.
       */
      void __elihw() {
	auto seq = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	context.pop();
	context.top().push_back(xsom::instr::Instr(new xsom::instr::While(seq,tests.top())));
	tests.pop();
      }

      /**
       * Update the architecture
       */
      void __update() {
	if(archi != null_ptr)
	  __step([archi](){archi->update();});
      }

      /**
       * Learn the architecture
       */
      void __learn() {
	if(archi != null_ptr)
	  __step([archi](){archi->learn();});
      }

      /**
       * Update and learn the architecture
       */
      void __update_and_learn() {
	if(archi != null_ptr)
	  __step([archi](){archi->update(); archi->learn();});
      }
      
      /**
       * Runs the whole program of the sequencer.
       */
      void run() {
	if(context.size() != 1)
	  std::cerr << "Sequencer error : there are badly closed instructions" << std::endl;
	else {
	  auto main = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	  try {
	    while(true)
	      main->next();
	  }
	  catch(xsom::instr::Done& e) {}
	  if(display != nullptr)
	    std::cout << ccmpl::stop;
	}
      }
    };

    Sequencer sequencer(xsom::Container& archi,
			ccmpl::chart::Layout& display) {return Sequencer(archi, display);}
    Sequencer sequencer(xsom::Container& archi)        {return Sequencer(archi);}
    Sequencer sequencer(ccmpl::chart::Layout& display) {return Sequencer(display);}
    Sequencer sequencer()                              {return Sequencer();}
    
  }
}
