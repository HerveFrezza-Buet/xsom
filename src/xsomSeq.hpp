#pragma once

#include <xsomNetwork.hpp>
#include <ccmpl.hpp>
#include <map>
#include <iterator>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>
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


    class RepUntil : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      Instr body;
      std::function<bool ()> test;

      template<typename TEST>
      RepUntil(Instr body, const TEST& test) : Instruction(), body(body), test(test) {}

    protected:
      
      virtual void execute() {
	try {
	  body->next();
	}
	catch(Done e) {
	  has_next = !(test());
	  if(has_next) {
	    body = body->deep_copy();
	    execute();
	  }
	}
      }

    public:
      
      virtual Instr deep_copy() {
	return Instr(new RepUntil(body->deep_copy(), test));
      }
      
    };

    

    class If : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      Instr body_then;
      Instr body_else;
      Instr body;
      std::function<bool ()> test;
      bool first;

      template<typename TEST>
      If(const TEST& test, Instr bthen, Instr belse)
	: Instruction(), body_then(bthen), body_else(belse), body(nullptr), test(test), first(true) {}

    protected:
      
      virtual void execute() {
	if(first) {
	  if(test())
	    body = body_then;
	  else
	    body = body_else;
	}

	has_next = body != nullptr;

	if(has_next) {
	  try {
	    body->next();
	  }
	  catch(Done e) {
	    has_next = false;
	  }
	}
      }

    public:
      
      virtual Instr deep_copy() {
	Instr bt = body_then;
	Instr be = body_else;
	if(bt != nullptr) bt = body_then->deep_copy();
	if(be != nullptr) be = body_else->deep_copy();
	return Instr(new If(test, bt, be));
      }
      
    };
    
  }

  namespace setup {

    /**
     * Sequencer offers a mini-programming language in order to schedule the computation.
     * @code{.cpp}
     * // auto seq = xsom::setup::sequencer();
     * // auto seq = xsom::setup::sequencer(archi);
     * // auto seq = xsom::setup::sequencer(display);
     * auto seq = xsom::setup::sequencer(archi,display);
     * 
     * // Instruction
     * seq.__step(f); // void f()
     *
     * // Predefined print on cerr
     * seq.__print("Hello world !");
     *
     * // Print and increment the counter "cntr". It displays something 
     * // like "value is 000054/100.", starting from "value is 000015/100."
     * seq.__counter("cntr", "value is ", 15, 100);
     *
     * // Sequence (not really useful)
     * seq.__begin();
     * ...
     * seq.__end();
     *
     * // Macro definition and call
     * seq.__def("foo")
     * ...
     * seq.__fed()
     * ...
     * seq.__call("foo")
     *
     * // If
     * seq.__if(check); // bool check()
     * ... // then code here, can be empty.
     * seq.__else();
     * ... // else code here, can be empty.
     * seq.__fi();
     *
     * // Infinite loop
     * seq.__loop();
     * ...
     * seq.__pool();
     *
     * // Finite loop (repeate n times)
     * seq.__for(n);
     * ...
     * seq.__rof();
     *
     * // While loop 
     * seq.__while(check); // bool check()
     * ...
     * seq.__elihw();
     *
     * // Repeat until loop
     * seq.__repeat();
     * ...
     * seq.__until(check); // bool check()
     *
     * // Run the architecture
     * seq.__update();
     * seq.__learn();
     * seq.__update_and_learn();
     *
     * // Plot the data.  std::string flags()
     * seq.__plot(flags);  
     * seq.__plot(flags, "pdf-frame" "png-frame");
     * seq.__plot_pdf(flags, "pdf-frame");
     * seq.__plot_png(flags, "png-frame");
     *
     * // save and load
     * seq.__on_load(load_fct); // void load_fct(std::istream&)
     * seq.__on_save(save_fct); // void save_fct(std::ostream&)
     * ...
     * seq.__load("status.data"); 
     * seq.__save("status.data"); 
     * seq.__autosave("status", "data"); // saves "status-000000.data"
     * seq.__autosave("status", "data"); // saves "status-000001.data"
     * seq.__autosave("status", "data"); // saves "status-000002.data"
     * ...
     *
     *
     * @endcode
     */
    class Sequencer {
    private:
      xsom::Container* archi;
      ccmpl::chart::Layout* display;

      std::map<std::string,unsigned int>         frame_pdf;
      std::map<std::string,unsigned int>         frame_png;
      std::map<std::string,unsigned int>         frame_file;
      std::map<std::string,unsigned int>         frame_counter;
      std::map<std::string,xsom::instr::Instr>   macros;
      std::stack<std::list<xsom::instr::Instr>>  context;
      std::stack<std::string>                    macro_names;
      std::stack<unsigned int>                   for_times;
      std::stack<std::function<bool ()>>         tests;
      std::stack<xsom::instr::Instr>             thens;
      std::function<void (std::ostream&)>        save_fct;
      std::function<void (std::istream&)>        load_fct;
      
      static std::string next(std::map<std::string,unsigned int>& frame,
			      unsigned int offset,
			      const std::string& name,
			      const std::string& prefix,
			      const std::string& suffix) {
	auto kv = frame.find(name);
	unsigned int id = 0;
	if(kv == frame.end())
	  frame[name] = 1;
	else
	  id = (kv->second)++;
	std::ostringstream ostr;
	ostr << prefix << std::setw(6) << std::setfill('0') << id+offset << suffix;
	return ostr.str();
      }
      
      static std::string next(std::map<std::string,unsigned int>& frame,
			      const std::string& name,
			      const std::string& prefix,
			      const std::string& suffix) {
	return next(frame, 0, name, prefix, suffix);
      }

      std::string next_counter(const std::string& name, unsigned int start, unsigned int max) {
	return next(frame_counter, start, name, "", std::string("/")+std::to_string(max)+".");
      }
      
      std::string next_filename(const std::string& prefix,
				const std::string& suffix) {
	return next(frame_file, prefix, prefix+"-", std::string(".")+suffix);
      }


      std::string next_pdf(const std::string& name) {
	return next(frame_pdf, name, name+"-", ".pdf");
      }

      std::string next_png(const std::string& name) {
	return next(frame_png, name, name+"-", ".png");
      }

      void print_save_info(const std::string& filename) {
	std::cerr << "Sequencer info : \"" << filename << "\" saved." << std::endl;
      }
      
      void print_load_info(const std::string& filename) {
	std::cerr << "Sequencer info : \"" << filename << "\" loaded." << std::endl;
      }


      void load(const std::string& filename) {
	if(this->load_fct) {
	  std::ifstream file;
	  file.open(filename);
	  if(file) {
	    this->load_fct(file);
	    file.close();
	    print_load_info(filename);
	  }
	  else
	    std::cerr << "Sequencer error : cannot open \"" << filename << "\" for reading." << std::endl;
	}
	else
	  std::cerr << "Sequencer error : no load function defined (have you called __on_load ?)."<< std::endl;
	  
      }


      void save(const std::string& filename) {
	if(this->save_fct) {
	  std::ofstream file;
	  file.open(filename);
	  if(file) {
	    this->save_fct(file);
	    file.close();
	    print_save_info(filename);
	  }
	  else
	    std::cerr << "Sequencer error : cannot open \"" << filename << "\" for writing." << std::endl;
	}
	else
	  std::cerr << "Sequencer error : no save function defined (have you called __on_save ?)."<< std::endl;
	  
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
       * Ingrements a counter and prints it.
       */
      void __counter(const std::string& name, const std::string& prefix, unsigned int start, unsigned int max) {
	__step([name, prefix, start, max, this](){std::cerr << prefix << this->next_counter(name, start, max) << std::endl;});
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
       * close the __while call.
       */
      void __elihw() {
	auto seq = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	context.pop();
	context.top().push_back(xsom::instr::Instr(new xsom::instr::While(seq,tests.top())));
	tests.pop();
      }

      /**
       * Define an Repeat. 
       */
      void __repeat() {
	context.push(std::list<xsom::instr::Instr>());
      }
      
      /**
       * Define an until statement, closing __repeat. test is "bool test(). 
       */
      template<typename TEST>
      void __until(const TEST& test) {
	auto seq = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	context.pop();
	context.top().push_back(xsom::instr::Instr(new xsom::instr::RepUntil(seq,test)));
      }
      
      /**
       * if
       */
      template<typename TEST>
      void __if(const TEST& test) {
	context.push(std::list<xsom::instr::Instr>());
	tests.push(test);
      }
      
      /**
       * else
       */
      void __else() {
	xsom::instr::Instr bthen = nullptr;
	if(context.top().size() > 0)
	  bthen = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	thens.push(bthen);
	context.pop();
	context.push(std::list<xsom::instr::Instr>());
      }
      
      /**
       * close the __if call.
       */
      void __fi() {
	xsom::instr::Instr belse = nullptr;
	if(context.top().size() > 0)
	  belse = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	context.pop();
	context.top().push_back(xsom::instr::Instr(new xsom::instr::If(tests.top(), thens.top(), belse)));
	tests.pop();
	thens.pop();
      }

      /**
       * Update the architecture
       */
      void __update() {
	if(archi != nullptr)
	  __step([archi = this->archi](){archi->update();});
      }

      /**
       * Learn the architecture
       */
      void __learn() {
	if(archi != nullptr)
	  __step([archi = this->archi](){archi->learn();});
      }

      /**
       * Update and learn the architecture
       */
      void __update_and_learn() {
	if(archi != nullptr)
	  __step([archi = this->archi](){archi->update(); archi->learn();});
      }

      /**
       * Plot. flags is "bool flags()"
       */
      template<typename FLAGS>
      void __plot(const FLAGS& flags) {
	if(display != nullptr) {
	  std::function<std::string ()> f(flags);
	  __step([display = this->display, f]() {std::cout << (*display)(f(), ccmpl::nofile(), ccmpl::nofile());});
	}
      }
      
      /**
       * Plot. flags is "bool flags()"
       */
      template<typename FLAGS>
      void __plot_png(const FLAGS& flags, const std::string& png_prefix) {
	if(display != nullptr) {
	  std::function<std::string ()> f(flags);
	  __step([this, f, png_prefix]() {
	      auto png = this->next_png(png_prefix);
	      std::cout << (*(this->display))(f(), ccmpl::nofile(), png);
	      this->print_save_info(png);
	    });
	}
      }
      
      /**
       * Plot. flags is "bool flags()"
       */
      template<typename FLAGS>
      void __plot_pdf(const FLAGS& flags, const std::string& pdf_prefix) {
	if(display != nullptr) {
	  std::function<std::string ()> f(flags);
	  __step([this, f, pdf_prefix]() {
	      auto pdf = this->next_pdf(pdf_prefix);
	      std::cout << (*(this->display))(f(), pdf, ccmpl::nofile());
	      this->print_save_info(pdf);
	    });
	}
      }
      
      /**
       * Plot. flags is "bool flags()"
       */
      template<typename FLAGS>
      void __plot(const FLAGS& flags, const std::string& pdf_prefix, const std::string& png_prefix) {
	if(display != nullptr) {
	  std::function<std::string ()> f(flags);
	  __step([this, f, pdf_prefix, png_prefix]() {
	  auto pdf = this->next_pdf(pdf_prefix);
	  auto png = this->next_png(png_prefix);
	      std::cout << (*(this->display))(f(), pdf, png);
	      this->print_save_info(pdf);
	      this->print_save_info(png);
	    });
	}
      }

      /**
       * Set load function : void f(std::istream&)
       */
      template<typename F>
      void __on_load(const F& f) {
	std::function<void (std::istream&)> l(f);
	__step([this,l](){this->load_fct = l;});
      }

      /**
       * Set save function : void f(std::ostream&)
       */
      template<typename F>
      void __on_save(const F& f) {
	std::function<void (std::ostream&)> s(f);
	__step([this,s](){this->save_fct = s;});
      }

      /**
       * Load file
       */
      void __load(const std::string& filename) {
	__step([this, filename](){this->load(filename);});
      }

      /**
       * Save file
       */
      void __save(const std::string& filename) {
	__step([this, filename](){this->save(filename);});
      }
      
      /**
       * Save file with prefix-xxxxxx.suffix format.
       */
      void __autosave(const std::string& prefix, const std::string& suffix) {
	__step([this, prefix, suffix](){this->save(this->next_filename(prefix, suffix));});
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
