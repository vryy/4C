
#include "cut_creator.H"
#include "cut_mesh.H"

void GEO::CUT::Creator::NewLine( Point * p1, Point * p2, Side * side, Side * other, Element * element )
{
  lines_.push_back( LineCmd() );
  LineCmd & line = lines_.back();
  line.Assign( p1, p2, side, other, element );
}

void GEO::CUT::Creator::Execute( Mesh & mesh )
{
  for ( std::vector<LineCmd>::iterator i=lines_.begin(); i!=lines_.end(); ++i )
  {
    LineCmd & line = *i;
    line.Execute( mesh );
  }
}

void GEO::CUT::Creator::LineCmd::Assign( Point * p1, Point * p2, Side * side, Side * other, Element * element )
{
  p1_ = p1;
  p2_ = p2;
  side_ = side;
  other_ = other;
  element_ = element;
}

void GEO::CUT::Creator::LineCmd::Execute( Mesh & mesh )
{
  mesh.NewLine( p1_, p2_, side_, other_, element_ );
}
