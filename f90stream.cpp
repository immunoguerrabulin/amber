extern "C"
{
  void f90stream_cout( int const n, char const * string );
  int  f90stream_get_cout_unit();
  void f90stream_set_cout_unit( int const iunit );
  void f90stream_flush_cout();
}

#include <sstream>
#include "f90stream.hpp"

class F90streamBuffer : public std::stringbuf
{
public:
  virtual int sync() 
  {
    int const n = (int)this->str().size();
    if ( n > 0 )
      {
	f90stream_cout( n, this->str().c_str() );
	this->str("");
      };
    return 0;
  }
  virtual ~F90streamBuffer() { sync(); }
};

namespace f90stream
{
  F90streamBuffer buffer_;
}
std::ostream f90stream::cout( &f90stream::buffer_ );

int f90stream::get_cout_unit()
{
  return f90stream_get_cout_unit();
}

void f90stream::set_cout_unit( int const iunit )
{
  f90stream_set_cout_unit( iunit );
}

void f90stream_flush_cout()
{
  f90stream::cout.flush();
}
