#ifndef CUT_UTILS_H
#define CUT_UTILS_H

#include <iostream>
#include <iterator>
#include <vector>
#include <cstring>

namespace GEO
{
  namespace CUT
  {

    inline void DumpDoubles( std::ostream & stream, const double * data, int length )
    {
      std::vector<int> c( length*sizeof( double ) / sizeof( int ) );
      std::memcpy( &c[0], data, c.size()*sizeof( int ) );
      stream << "  {";
      std::copy( c.begin(), c.end(), std::ostream_iterator<int>( stream, "," ) );
      stream << "};";
    }

  }
}

#endif
