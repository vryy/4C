#ifndef CELL_LINE_H
#define CELL_LINE_H

#include "cell_point.H"

namespace GEO
{
namespace CELL
{
  class Facet;

class Line
{
public:

  Line( Point * p1, Point * p2 ) : p1_( p1 ), p2_( p2 ) {}

  bool MyPoint( Point * p )
  {
    return p1_==p or p2_==p;
  }

private:
  Point * p1_;
  Point * p2_;

  std::set<Facet*> facets_;
};

}
}

#endif
