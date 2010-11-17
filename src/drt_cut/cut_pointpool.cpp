
#include "cut_pointpool.H"
#include "cut_tolerance.H"

GEO::CUT::Point* GEO::CUT::OctTreeNode::NewPoint( const double * x, Edge * cut_edge, Side * cut_side )
{
  Teuchos::RCP<Point> p = GetPoint( x, cut_edge, cut_side );
  if ( p==Teuchos::null )
  {
    p = CreatePoint( points_.size(), x, cut_edge, cut_side );
#if 1
    if ( points_.size()%1000 == 0 )
    {
      Split( 0 );
    }
#endif
  }
  return &*p;
}

Teuchos::RCP<GEO::CUT::Point> GEO::CUT::OctTreeNode::GetPoint( const double * x, Edge * cut_edge, Side * cut_side )
{
  if ( not IsLeaf() )
  {
    for ( int i=0; i<8; ++i )
    {
      Teuchos::RCP<Point> p = nodes_[i]->GetPoint( x, cut_edge, cut_side );
      if ( p!=Teuchos::null )
      {
        return p;
      }
    }
  }
  else
  {
    LINALG::Matrix<3,1> px( x );
    LINALG::Matrix<3,1> nx;

    for ( std::set<Teuchos::RCP<Point>, PointPidLess>::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;

      p->Coordinates( nx.A() );
      nx.Update( -1, px, 1 );
      if ( nx.Norm2() < MINIMALTOL )
      {
        if ( cut_edge!=NULL )
        {
          p->AddEdge( cut_edge );
        }
        if ( cut_side!=NULL )
        {
          p->AddSide( cut_side );
        }
        return *i;
      }
    }
  }
  return Teuchos::null;
}

Teuchos::RCP<GEO::CUT::Point> GEO::CUT::OctTreeNode::CreatePoint( unsigned newid, const double * x, Edge * cut_edge, Side * cut_side )
{
  if ( not IsLeaf() )
  {
    Teuchos::RCP<Point> p = Leaf( x )->CreatePoint( newid, x, cut_edge, cut_side );
    AddPoint( x, p );
    return p;
  }
  else
  {
    Teuchos::RCP<Point> p = Teuchos::rcp( new Point( newid, x, cut_edge, cut_side ) );
    AddPoint( x, p );
    return p;
  }
}

GEO::CUT::OctTreeNode* GEO::CUT::OctTreeNode::Leaf( const double * x )
{
  int idx = 0;
  if ( x[0] > splitpoint_( 0 ) )
  {
    idx += 1;
  }
  if ( x[1] > splitpoint_( 1 ) )
  {
    idx += 2;
  }
  if ( x[2] > splitpoint_( 2 ) )
  {
    idx += 4;
  }
  return &*nodes_[idx];
}

void GEO::CUT::OctTreeNode::Split( int level )
{
  // We must not end up with a OctTreeNode that holds just nodes from the
  // cutter mesh. However, there is no real way to test this right now.

  if ( points_.size()>125 )
  {
    LINALG::Matrix<3,1> x;
    bool first = true;

    for ( std::set<Teuchos::RCP<Point>, PointPidLess>::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;
      if ( first )
      {
        first = false;
        p->Coordinates( splitpoint_.A() );
      }
      else
      {
        p->Coordinates( x.A() );
        splitpoint_.Update( 1, x, 1 );
      }
    }

    splitpoint_.Scale( 1./points_.size() );

    for ( int i=0; i<8; ++i )
    {
      nodes_[i] = Teuchos::rcp( new OctTreeNode() );
    }

    // avoid empty room (room not covered by boundary boxes)
    for ( int i=0; i<8; ++i )
    {
      // always have the split point in all boxes
      nodes_[i]->bb_.AddPoint( splitpoint_ );

      // always have the outmost point in each box
      double x[3];
      bb_.CornerPoint( i, x );
      Leaf( x )->bb_.AddPoint( x );
    }

    for ( std::set<Teuchos::RCP<Point>, PointPidLess>::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Teuchos::RCP<Point> p = *i;
      double x[3];
      p->Coordinates( x );
      Leaf( x )->AddPoint( x, p );
    }

    for ( int i=0; i<8; ++i )
    {
      nodes_[i]->Split( level+1 );
    }
  }
}

void GEO::CUT::OctTreeNode::AddPoint( const double * x, Teuchos::RCP<Point> p )
{
  points_.insert( p );
  bb_.AddPoint( x );
}

void GEO::CUT::OctTreeNode::CollectSides( const BoundingBox & sidebox, std::set<Side*> & sides )
{
  if ( not IsLeaf() )
  {
    if ( sidebox.Within( bb_ ) )
    {
      for ( int i=0; i<8; ++i )
      {
        nodes_[i]->CollectSides( sidebox, sides );
      }
    }
  }
  else
  {
    BoundingBox sbox;
    for ( std::set<Teuchos::RCP<Point>, PointPidLess>::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;
      const std::set<Side*> & sds = p->CutSides();
      for ( std::set<Side*>::iterator i=sds.begin(); i!=sds.end(); ++i )
      {
        Side * s = *i;
        if ( sides.count( s )==0 )
        {
          sbox.Assign( *s );
          if ( sbox.Within( sidebox ) )
          {
            sides.insert( s );
          }
        }
      }
    }
  }
}

void GEO::CUT::OctTreeNode::CollectElements( const BoundingBox & sidebox, std::set<Element*> & elements )
{
  if ( not IsLeaf() )
  {
    if ( sidebox.Within( bb_ ) )
    {
      for ( int i=0; i<8; ++i )
      {
        nodes_[i]->CollectElements( sidebox, elements );
      }
    }
  }
  else
  {
    BoundingBox elementbox;
    for ( std::set<Teuchos::RCP<Point>, PointPidLess>::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;
      const std::set<Element*> & els = p->Elements();
      for ( std::set<Element*>::iterator i=els.begin(); i!=els.end(); ++i )
      {
        Element * e = *i;
        if ( elements.count( e )==0 )
        {
          elementbox.Assign( *e );
          if ( elementbox.Within( sidebox ) )
          {
            elements.insert( e );
          }
        }
      }
    }
  }
}

void GEO::CUT::OctTreeNode::ResetOutsidePoints()
{
  for ( std::set<Teuchos::RCP<Point>, PointPidLess>::iterator i=points_.begin();
        i!=points_.end();
        ++i )
  {
    Point * p = &**i;
    if ( p->Position()==Point::outside )
    {
      p->Position( Point::undecided );
    }
  }
}
