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
#include <stdexcept>

#include <unistd.h>
#include <fcntl.h>

#include <ncurses.h>

namespace xsom {
  namespace setup {
    class Sequencer;
  }

  namespace msg {
      struct ColorTag {
	int tag;
	ColorTag(int tag) : tag(tag) {}
	friend std::ostream& operator<<(std::ostream& os, const ColorTag& ct) {
	  return os << "\033[" << ct.tag << "m";
	}
      };


      inline std::ostream& deflt(std::ostream& os) {
	os << ColorTag(39);
	return os;
      }
      
      inline std::ostream& red(std::ostream& os) {
	os << ColorTag(31);
	return os;
      }
      
      inline std::ostream& green(std::ostream& os) {
	os << ColorTag(32);
	return os;
      }
      
      inline std::ostream& yellow(std::ostream& os) {
	os << ColorTag(33);
	return os;
      }
      
      inline std::ostream& blue(std::ostream& os) {
	os << ColorTag(34);
	return os;
      }
      
      inline std::ostream& magenta(std::ostream& os) {
	os << ColorTag(35);
	return os;
      }
      
      inline std::ostream& cyan(std::ostream& os) {
	os << ColorTag(36);
	return os;
      }
      
      inline std::ostream& endl(std::ostream& os) {
	os << deflt << std::endl;
	return os;
      }

      inline std::ostream& seq_error(std::ostream& os) {
	os << red     << "Sequencer error : ";
	return os;
      }

      inline std::ostream& seq_file_info(std::ostream& os) {
	os << magenta << "Sequencer file  : ";
	return os;
      }
      
      inline std::ostream& seq_cntr_info(std::ostream& os) {
	os << blue    << "Sequencer count : ";
	return os;
      }
      
      inline std::ostream& seq_value_info(std::ostream& os) {
	os << cyan    << "Sequencer value : ";
	return os;
      }
      
      inline std::ostream& seq_msg_info(std::ostream& os) {
	os << green   << "Sequencer info  : ";
	return os;
      }

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

    class KeyboardInteraction : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      xsom::setup::Sequencer* owner;
      KeyboardInteraction(xsom::setup::Sequencer* owner) : Instruction(), owner(owner) {}

    protected:
      
      virtual void execute(); // see definition at the end of the this file.

    public:
      
      virtual Instr deep_copy() {
	return Instr(new KeyboardInteraction(owner));
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
      Instr instruction;
      Call(Instr instruction) : Instruction(), instruction(instruction->deep_copy()) {}

    protected:
      
      virtual void execute() {
	try {
	  instruction->next();
	}
	catch(Done) {
	  has_next = false;
	}
      }

    public:
      
      virtual Instr deep_copy() {
	return Instr(new Call(instruction)); // deep copy is performed by Call(...).
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
      xsom::setup::Sequencer* owner;
      Instr body;
      
      Loop(xsom::setup::Sequencer* owner, Instr body) : Instruction(),  owner(owner), body(body) {}

    protected:
      
      virtual void execute();

    public:
      
      virtual Instr deep_copy() {
	return Instr(new Loop(owner, body->deep_copy()));
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
     * // Clear/Reset the counter.
     * seq.__counter_clear("cntr");
     *
     * // prints a numerical value. double f().
     * seq.__value("value is", f);
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
     * // Keyboard interation (if seq.interactive() has been called initially)
     * seq.__();
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

      friend class xsom::instr::KeyboardInteraction;

      bool inter = false;
      bool cont_mode = false;
      std::ofstream pipe;
      std::ostream* pipe_ptr;
      
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
      
      static void clear(std::map<std::string,unsigned int>& frame,
			       const std::string& name) {
	auto kv = frame.find(name);
	if(kv != frame.end())
	  frame.erase(kv);
      }
      
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
      
      void clear_counter(const std::string& name) {
	return clear(frame_counter, name);
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
	if(inter) 
	  inter_msg("Save", std::string("\"") + filename + std::string("\""));
	else
	  std::cerr << msg::seq_file_info << "\"" << filename << "\" saved." << msg::endl;
      }
      
      void print_load_info(const std::string& filename) {
	if(inter) 
	  inter_msg("Load", std::string("\"") + filename + std::string("\""));
	else
	  std::cerr << msg::seq_file_info << "\"" << filename << "\" loaded." << msg::endl;
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
	    std::cerr << msg::seq_error << "cannot open \"" << filename << "\" for reading." << msg::endl;
	}
	else
	  std::cerr << msg::seq_error << "no load function defined (have you called __on_load ?)."<< msg::endl;
	  
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
	    std::cerr << msg::seq_error << "cannot open \"" << filename << "\" for writing." << msg::endl;
	}
	else
	  std::cerr << msg::seq_error << "no save function defined (have you called __on_save ?)."<< msg::endl;
	  
      }

      
      Sequencer(xsom::Container* archi, ccmpl::chart::Layout* display)
	: archi(archi), display(display) {
	context.push(std::list<xsom::instr::Instr>());
	pipe_ptr = &(std::cout);
      }

      void inter_menu() {
	::clear();
	mvprintw(0, 0, "Key bindings");
	mvprintw(1, 0, "  <space> : next/pause");
	mvprintw(2, 0, "  c       : cont");
	mvprintw(3, 0, "  ESC     : quit");
	move    (4, 0);
      }
		       

      void inter_msg(const std::string& tag, const std::string& message) {
	inter_menu();
	mvprintw(5, 0, (std::string(tag) + " : " + message).c_str());
	move    (7, 0);
	refresh();
      }

    public:

      Sequencer(const Sequencer&) = delete;
      Sequencer& operator=(const Sequencer&) = delete;

      Sequencer(xsom::Container& archi, ccmpl::chart::Layout& display)
	: Sequencer(&archi, &display) {}

      Sequencer(xsom::Container& archi)
	: Sequencer(&archi, nullptr) {}
    
      Sequencer(ccmpl::chart::Layout& display)
	: Sequencer(nullptr, &display) {}

      Sequencer()
	: Sequencer(nullptr,nullptr) {}

      /**
       * This displays an error message
       */
      void msg_error(const std::string& message) {
	if(inter)
	  inter_msg("Error", message);
	else
	  std::cerr << msg::seq_error << message << msg::endl;
      }
      
      /**
       * This displays an information message
       */
      void msg_info(const std::string& message) {
	if(inter)
	  inter_msg("Info", message);
	else
	  std::cerr << msg::seq_msg_info << message << msg::endl;
      }

      /**
       * This sets the sequencer into a keybord interaction mode.
       * @step_mode Tells wether the interactive execution is started in setp-by-step mode or not.
       * @pipename The name of the system named pipe used for data exchange.
       */
      void interactive(bool step_mode, std::string pipename) {
	cont_mode = !step_mode;
	pipe.open(pipename.c_str(), std::fstream::app);
	if(!pipe)
	  throw std::runtime_error(std::string("Cannot open pipe \"") + pipename + "\".");
	pipe_ptr = &pipe;

	// Let us remove pipe buffer for more reactivity.
	auto f = open(pipename.c_str(), O_WRONLY);
	fcntl(f, F_SETPIPE_SZ, 0L);
	close(f);
      
	
	inter = true;

	initscr();
	noecho();
	raw();
	nodelay(stdscr, cont_mode);
	keypad(stdscr, TRUE);

	inter_menu();
	refresh();
      }

      /**
       * Add a step that calls a macro.
       */
      void __call(const std::string& macro_name) {
	auto kv = macros.find(macro_name);
	if(kv == macros.end())
	  std::cerr << msg::seq_error << "macro \"" << macro_name
		    << "\" undefined." << msg::endl;
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
       * Add a keyboard interaction point
       */
      void __() {
	context.top().push_back(xsom::instr::Instr(new xsom::instr::KeyboardInteraction(this)));
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
	__step([message](){std::cerr << msg::seq_msg_info << message << msg::endl;});
      }

      /**
       * Increments a counter and prints it.
       */
      void __counter(const std::string& name, const std::string& prefix, unsigned int start, unsigned int max) {
	__step([name, prefix, start, max, this](){std::cerr << msg::seq_cntr_info << prefix << this->next_counter(name, start, max) << msg::endl;});
      }
      
      /**
       * Clears a counter.
       */
      void __counter_clear(const std::string& name) {
	__step([name, this](){this->clear_counter(name);});
      }

      /**
       * Compute and display a value
       */
      template<typename Fun>
      void __value(const std::string& message, const Fun& f) {
	auto ff = std::bind(f);
	__step([message, ff]() {std::cerr << msg::seq_value_info << message << ff() << msg::endl;});
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
	context.top().push_back(xsom::instr::Instr(new xsom::instr::Loop(this, seq)));
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
	  __step([display = this->display, f, &os = *pipe_ptr]() {os << (*display)(f(), ccmpl::nofile(), ccmpl::nofile()) << std::flush;});
	}
      }
      
      /**
       * Plot. flags is "bool flags()"
       */
      template<typename FLAGS>
      void __plot_png(const FLAGS& flags, const std::string& png_prefix) {
	if(display != nullptr) {
	  std::function<std::string ()> f(flags);
	  __step([this, f, png_prefix, &os = *pipe_ptr]() {
	      auto png = this->next_png(png_prefix);
	      os << (*(this->display))(f(), ccmpl::nofile(), png) << std::flush;
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
	  __step([this, f, pdf_prefix, &os = *pipe_ptr]() {
	      auto pdf = this->next_pdf(pdf_prefix);
	      os << (*(this->display))(f(), pdf, ccmpl::nofile()) << std::flush;
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
	  __step([this, f, pdf_prefix, png_prefix, &os = *pipe_ptr]() {
	      auto pdf = this->next_pdf(pdf_prefix);
	      auto png = this->next_png(png_prefix);
	      os << (*(this->display))(f(), pdf, png) << std::flush;
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
	  std::cerr << msg::seq_error << "there are badly closed instructions" << msg::endl;
	else {
	  auto main = xsom::instr::Instr(new xsom::instr::Seq(context.top()));
	  try {
	    while(true)
	      main->next();
	  }
	  catch(xsom::instr::Done& e) {}
	  if(display != nullptr)
	    (*pipe_ptr) << ccmpl::stop;
	  if(inter)
	    endwin(); // closing ncurse context.
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



inline void xsom::instr::KeyboardInteraction::execute() {
  if(owner->inter) {
    bool inloop = true;
    while(inloop) {
      inloop = false;
      int key = getch();
      switch(key) {
      case 'c' :
	owner->cont_mode = true;
	nodelay(stdscr, owner->cont_mode);
	break;
      case ' ' :
	owner->cont_mode = false;
	nodelay(stdscr, owner->cont_mode);
	break;
      case 'q' : // ESC
       	throw Done();
      default:
	break;
      }
    }
  }
  
  has_next = false;
}

inline void xsom::instr::Loop::execute() {
  try {
    body->next();
  }
  catch(Done e) {
    body = body->deep_copy();
    try {
      body->next();
    }
    catch(Done e) {
      owner->msg_error("loop has an invalid instruction");
    }
  }
}
