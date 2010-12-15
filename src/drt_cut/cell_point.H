#ifndef CELL_POINT_H
#define CELL_POINT_H

#include <algorithm>
#include <set>

namespace GEO
{
namespace CELL
{
  class Line;

class Point
{
public:

  Point( double * x )
  {
    std::copy( x, x+3, x_ );
  }

  Line * FindLine( Point * other );

private:
  double x_[3];

  std::set<Line*> lines_;

};

}
}

#endif
