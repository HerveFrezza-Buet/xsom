#pragma once

#include <xsomNetwork.hpp>
#include <ccmpl.hpp>
#include <map>
#include <string>
#include <sstring>
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
   
    
    using Instr = std::shared_ptr<const Instruction>;

    struct Done {
      Done() {}
    };

    class Instruction {
    protected:
      bool has_next;
      Instruction() : has_next(true) {}
      virtual void execute() = 0;

    public:
      
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
    };

    class Call : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      Instr instr;
      Call(Instr instr) : Instruction(), instr(instr) {}

    protected:
      
      virtual void execute() {
	try {
	  instr->next();
	}
	catch(Done) {
	  has_next = false;
	}
      }
    };
    
    class Seq : public Instruction {
    private:
      friend class xsom::setup::Sequencer;
      std::list<Instr> seq;
      std::list<Instr>::iterator current;
      
      Step(const std::list<Instr>& seq) : Instruction(), seq(seq), current(seq.begin()) {}

    protected:

      virtual void execute() {
	while(current != seq.end()) {
	  try {
	    current->next();
	    break;
	  }
	  catch(Done d) {
	    ++current;
	  }
	}
	has_next = current != seq.end();
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
	return osrt.str();
      }

      std::sring next_pdf(const std::string& name) {
	return next(frame_pdf, name, "pdf");
      }

      std::sring next_png(const std::string& name) {
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
	  context.top().push_back(new xsom::instr::Call(kv->second));
      }

      /**
       * Add a step that calls f. f is "void f()".
       */
      template<typename Fun>
      void __step(const Fun& f) {
	context.top().push_back(new xsom::instr::Step(f));
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
	auto seq = xsom::instr::Instr(new xsom::Instr::Seq(context.top()));
	context.pop();
	context.top().push_back(seq);
      }

      /**
       * Add a step that prints a message on cerr
       */
      void __print(const std::string& message) {
	__step([message](){std::cerr << "Sequencer info : " << message << std::endl;})
      }
      
      /**
       * Runs the whole program of the sequencer.
       */
      void run() {
	if(context.size() != 1)
	  std::cout << "Sequencer error : there are badly closed instructions" << std::endl;
	else {
	  auto main = xsom::instr::Instr(new xsom::Instr::Seq(context.top()));
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
