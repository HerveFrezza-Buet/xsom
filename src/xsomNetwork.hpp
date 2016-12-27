#pragma once

#include <xsomEntity.hpp>

namespace xsom {
  namespace setup {

    inline xsom::Container network() {return Container(false);}
    inline xsom::Debug debug(const string& name) {return Debug(name);}

  }
}
