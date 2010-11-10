
#include "cut_line.H"
#include "cut_side.H"

GEO::CUT::Line::Line( Point * p1, Point * p2, Side * cut_side1, Side * cut_side2, Element * cut_element )
  : p1_( p1 ),
    p2_( p2 )
{
  if ( cut_side1 )
    AddSide( cut_side1 );
  if ( cut_side2 )
    AddSide( cut_side2 );
  if ( cut_element )
    AddElement( cut_element );

  p1_->Register( this );
  p2_->Register( this );

  // self-register at all sides

  const std::set<Side*> & cut_sides1 = p1->CutSides();
  const std::set<Side*> & cut_sides2 = p2->CutSides();

  std::set<Side*> sides;

  std::set_intersection( cut_sides1.begin(), cut_sides1.end(),
                         cut_sides2.begin(), cut_sides2.end(),
                         std::inserter( sides, sides.begin() ) );

  for ( std::set<Side*>::iterator i=sides.begin(); i!=sides.end(); ++i )
  {
    Side * s = *i;
    AddSide( s );
  }

#if 0
  std::vector<GEO::CUT::Edge*> edges = p1->CutEdges( p2 );
  for ( std::vector<GEO::CUT::Edge*>::iterator i=edges.begin(); i!=edges.end(); ++i )
  {
    Edge * e = *i;
    const std::set<Side*> & edge_sides = e->Sides();
    for ( std::set<Side*>::const_iterator i=edge_sides.begin(); i!=edge_sides.end(); ++i )
    {
      Side * s = *i;
      AddSide( s );
    }
  }
#endif
}

bool GEO::CUT::Line::IsInternalCut( Side * side )
{
  return cut_sides_.count( side )>0 and not side->OnEdge( this );
}
