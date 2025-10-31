#ifndef _f90stream_hpp_
#define _f90stream_hpp_

#include <ostream>

namespace f90stream
{
  extern std::ostream cout;
  int  get_cout_unit();
  void set_cout_unit( int const iunit );
}

#endif
