
#include "cell_point.H"
#include "cell_line.H"

GEO::CELL::Line * GEO::CELL::Point::FindLine( Point * other )
{
  for ( std::set<Line*>::iterator i=lines_.begin(); i!=lines_.end(); ++i )
  {
    Line * l = *i;
    if ( l->MyPoint( other ) )
    {
      return l;
    }
  }
  return NULL;
}
