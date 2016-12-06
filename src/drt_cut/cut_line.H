/*---------------------------------------------------------------------*/
/*!
\file cut_line.H

\brief cut line

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#ifndef CUT_LINE_H
#define CUT_LINE_H

#include "cut_point.H"

namespace GEO
{
namespace CUT
{

/*!
\begin Line between two points. These lines result from cuts and there are no cut points on a line.
 */
class Line
{
public:

  Line( Point * p1, Point * p2, Side * cut_side1, Side * cut_side2, Element * cut_element );

  void AddSide( Side * cut_side );

  void AddElement( Element * cut_element );

  bool IsCut( Side* s1, Side* s2 )
  {
    return cut_sides_.count( s1 )>0 and cut_sides_.count( s2 )>0;
  }

  bool IsCut( Element * element )
  {
    return cut_elements_.count( element ) > 0;
  }

  bool IsCut( Side * side )
  {
    return cut_sides_.count( side ) > 0;
//     return ( cut_sides_.count( side ) > 0 and
//              BeginPoint()->IsCut( side ) and
//              EndPoint()->IsCut( side ) );
  }

  bool IsInternalCut( Side * side );

  bool OnEdge( Edge * edge )
  {
    return p1_->IsCut( edge ) and p2_->IsCut( edge );
  }

  Point* OtherPoint( Point* point )
  {
    if ( p1_==point )
      return p2_;
    else if ( p2_==point )
      return p1_;
    else
      throw std::runtime_error( "foreign point provided" );
  }

  Point* BeginPoint()
  {
    return p1_;
  }

  Point* EndPoint()
  {
    return p2_;
  }

  bool Between( Point * p1, Point * p2 )
  {
    return ( ( p1==p1_ and p2==p2_ ) or
             ( p1==p2_ and p2==p1_ ) );
  }

  void Plot( std::ofstream & f )
  {
    f << "# line\n";
    p1_->Plot( f );
    p2_->Plot( f );
    f << "\n\n";
  }

  void Intersection( plain_side_set & sides )
  {
    plain_side_set intersection;
    std::set_intersection( cut_sides_.begin(), cut_sides_.end(),
                           sides.begin(), sides.end(),
                           std::inserter( intersection, intersection.begin() ) );
    std::swap( sides, intersection );
  }

  const plain_side_set & CutSides() const { return cut_sides_; }

private:

  Point* p1_;
  Point* p2_;

  plain_side_set cut_sides_;

  plain_element_set cut_elements_;
};

}
}

#endif
