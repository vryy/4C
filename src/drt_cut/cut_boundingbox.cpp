
#include "cut_boundingbox.H"
#include "cut_tolerance.H"
#include "cut_node.H"
#include "cut_edge.H"
#include "cut_side.H"
#include "cut_element.H"

GEO::CUT::BoundingBox::BoundingBox( Edge & edge )
  : empty_( true )
{
  const std::vector<Node*> & nodes = edge.Nodes();
  AddPoints( nodes );
}

GEO::CUT::BoundingBox::BoundingBox( Side & side )
  : empty_( true )
{
  const std::vector<Node*> & nodes = side.Nodes();
  AddPoints( nodes );
}

GEO::CUT::BoundingBox::BoundingBox( Element & element )
  : empty_( true )
{
  const std::vector<Node*> & nodes = element.Nodes();
  AddPoints( nodes );
}

void GEO::CUT::BoundingBox::Assign( Side & side )
{
  empty_ = true;
  const std::vector<Node*> & nodes = side.Nodes();
  AddPoints( nodes );
}

void GEO::CUT::BoundingBox::Assign( Element & element )
{
  empty_ = true;
  const std::vector<Node*> & nodes = element.Nodes();
  AddPoints( nodes );
}

void GEO::CUT::BoundingBox::AddPoints( const std::vector<Node*> & nodes )
{
  for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
  {
    Node * n = *i;
    double x[3];
    n->Coordinates( x );
    AddPoint( x );
  }
}

void GEO::CUT::BoundingBox::AddPoint( const double * x )
{
  if ( empty_ )
  {
    empty_ = false;
    box_( 0, 0 ) = box_( 0, 1 ) = x[0];
    box_( 1, 0 ) = box_( 1, 1 ) = x[1];
    box_( 2, 0 ) = box_( 2, 1 ) = x[2];
  }
  else
  {
    box_( 0, 0 ) = std::min( box_( 0, 0 ), x[0] );
    box_( 1, 0 ) = std::min( box_( 1, 0 ), x[1] );
    box_( 2, 0 ) = std::min( box_( 2, 0 ), x[2] );
    box_( 0, 1 ) = std::max( box_( 0, 1 ), x[0] );
    box_( 1, 1 ) = std::max( box_( 1, 1 ), x[1] );
    box_( 2, 1 ) = std::max( box_( 2, 1 ), x[2] );
  }
}

bool GEO::CUT::BoundingBox::Within( const BoundingBox & b ) const
{
  return  ( ( b( 0, 1 )-box_( 0, 0 ) ) > -1e-5 and ( b( 0, 0 )-box_( 0, 1 ) ) < 1e-5 and
            ( b( 1, 1 )-box_( 1, 0 ) ) > -1e-5 and ( b( 1, 0 )-box_( 1, 1 ) ) < 1e-5 and
            ( b( 2, 1 )-box_( 2, 0 ) ) > -1e-5 and ( b( 2, 0 )-box_( 2, 1 ) ) < 1e-5 );
}

bool GEO::CUT::BoundingBox::Within( const double * x ) const
{
  return  ( ( x[0]-box_( 0, 0 ) ) > -1e-5 and ( x[0]-box_( 0, 1 ) ) < 1e-5 and
            ( x[1]-box_( 1, 0 ) ) > -1e-5 and ( x[1]-box_( 1, 1 ) ) < 1e-5 and
            ( x[2]-box_( 2, 0 ) ) > -1e-5 and ( x[2]-box_( 2, 1 ) ) < 1e-5 );
}

bool GEO::CUT::BoundingBox::Within( const Epetra_SerialDenseMatrix & xyz ) const
{
  BoundingBox bb;
  int numnode = xyz.N();
  for ( int i=0; i<numnode; ++i )
  {
    bb.AddPoint( &xyz( 0, i ) );
  }
  return Within( bb );
}

bool GEO::CUT::BoundingBox::Within( Element & element ) const
{
  BoundingBox bb( element );
  return Within( bb );
}
