/*---------------------------------------------------------------------*/
/*!
\file cut_boundingbox.cpp

\brief bounding box for cut

\level 2

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include "cut_boundingbox.H"
#include "cut_element.H"
#include "cut_volumecell.H"

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

GEO::CUT::BoundingBox::BoundingBox( VolumeCell & volcell )
  : empty_( true )
{
  const plain_facet_set & facete = volcell.Facets();
  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
      Facet *fac = *i;
      const std::vector<Point*> & corners = fac->CornerPoints();
      for(std::vector<Point*>::const_iterator k=corners.begin();k!=corners.end();k++)
      {
          const Point* po = *k;
          const double * coords = po->X();
          AddPoint(coords);
      }
  }
}

/*--------------------------------------------------------------------------------------------*
    construct bounding box over the volumecell in local coordinates with respect to elem1
*---------------------------------------------------------------------------------------------*/
GEO::CUT::BoundingBox::BoundingBox( VolumeCell & volcell, Element *elem1 )
  : empty_( true )
{

  const plain_facet_set & facete = volcell.Facets();
  double x[3];
  for(plain_facet_set::const_iterator i=facete.begin();i!=facete.end();i++)
  {
      Facet* fac = *i;
      std::vector<std::vector<double> > corLocal;
      fac->CornerPointsLocal(elem1,corLocal);
      for(std::vector<std::vector<double> >::const_iterator m=corLocal.begin();m!=corLocal.end();m++)
      {
          std::vector<double> loc = *m;
          for(int j=0;j<3;j++)
                  x[j] = loc[j];
          AddPoint(x);
      }
  }
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

void GEO::CUT::BoundingBox::Assign( Edge & edge )
{
  empty_ = true;
  const std::vector<Node*> & nodes = edge.Nodes();
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

bool GEO::CUT::BoundingBox::Within( double norm, const BoundingBox & b ) const
{
  if ( empty_ )
    return true;
  return ( InBetween( norm, minx(), maxx(), b.minx(), b.maxx() ) and
           InBetween( norm, miny(), maxy(), b.miny(), b.maxy() ) and
           InBetween( norm, minz(), maxz(), b.minz(), b.maxz() ) );
}

bool GEO::CUT::BoundingBox::Within( double norm, const double * x ) const
{
  if ( empty_ )
    return true;
  return ( InBetween( norm, minx(), maxx(), x[0], x[0] ) and
           InBetween( norm, miny(), maxy(), x[1], x[1] ) and
           InBetween( norm, minz(), maxz(), x[2], x[2] ) );
}

bool GEO::CUT::BoundingBox::Within( double norm, const Epetra_SerialDenseMatrix & xyz ) const
{
  BoundingBox bb;
  int numnode = xyz.N();
  for ( int i=0; i<numnode; ++i )
  {
    bb.AddPoint( &xyz( 0, i ) );
  }
  return Within( norm, bb );
}

bool GEO::CUT::BoundingBox::Within( double norm, Element & element ) const
{
  BoundingBox bb( element );
  return Within( norm, bb );
}

void GEO::CUT::BoundingBox::Print()
{
  if ( empty_ )
  {
    std::cout << "  BB: {}\n";
  }
  else
  {
    std::cout << "  BB: {("
              << box_( 0, 0 ) << ","
              << box_( 1, 0 ) << ","
              << box_( 2, 0 ) << ")-("
              << box_( 0, 1 ) << ","
              << box_( 1, 1 ) << ","
              << box_( 2, 1 )
              << ")}\n";
  }
}

void GEO::CUT::BoundingBox::CornerPoint( int i, double * x )
{
  x[0] = ( ( i & 1 )==1 ) ? maxx() : minx();
  x[1] = ( ( i & 2 )==2 ) ? maxy() : miny();
  x[2] = ( ( i & 4 )==4 ) ? maxz() : minz();
}
