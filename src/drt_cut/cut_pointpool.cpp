/*---------------------------------------------------------------------*/
/*!
\file cut_pointpool.cpp

\brief PointPool, stores a points in the cut and decides if points are merged or new points are created

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15249
</pre>

*----------------------------------------------------------------------*/

#include "cut_pointpool.H"
#include "../drt_lib/drt_globalproblem.H"


/*-----------------------------------------------------------------------------------------*
 * If a point with the coordinates "x" does not exists, it creates a new point correspondingly
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::Point* GEO::CUT::OctTreeNode::NewPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance )
{
  // check if the point already exists
  Point * p = GetPoint( x, cut_edge, cut_side, tolerance );

  if ( p==NULL )
  {
    p = &*CreatePoint( points_.size(), x, cut_edge, cut_side, tolerance); // create the point
#if 1
    if ( points_.size()%1000 == 0 ) // split the node starting from level 0
    {
      Split( 0 );
    }
#endif
  }
  return p;
}


/*-----------------------------------------------------------------------------------------*
 * Get the point with the specified coordinates "x" from the pointpool
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::Point* GEO::CUT::OctTreeNode::GetPoint( const double * x, Edge * cut_edge, Side * cut_side, double tolerance )
{
  // try to find the point in all of the 8 children nodes
  if ( not IsLeaf() )
  {
    // stop finding the point when point is not included in the current bounding box
    if(!bb_->Within(1.0, x)) return NULL;

    for ( int i=0; i<8; ++i )
    {
      Point * p = nodes_[i]->GetPoint( x, cut_edge, cut_side, tolerance );
      if ( p!=NULL )
      {
        return p;
      }
    }
  }
  else
  {
    LINALG::Matrix<3,1> px( x );
    LINALG::Matrix<3,1> nx;

    double tol = tolerance*norm_;

    // linear search for the node in the current leaf
    for ( RCPPointSet::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;

      p->Coordinates( nx.A() );
      nx.Update( -1, px, 1 );
      if ( nx.Norm2() <= (tol + p->Tolerance()) )
      {
        if ( cut_edge!=NULL )
        {
          p->AddEdge( cut_edge );
        }
        if ( cut_side!=NULL )
        {
          p->AddSide( cut_side );
        }
        p->EnlargeTolerance(tol+nx.Norm2());
        return p;
      }
    }
  }
  return NULL;
}


/*-----------------------------------------------------------------------------------------*
 * Get the point with the specified coordinates "x" from the pointpool
 *-----------------------------------------------------------------------------------------*/
Teuchos::RCP<GEO::CUT::Point> GEO::CUT::OctTreeNode::CreatePoint( unsigned newid, const double * x, Edge * cut_edge, Side * cut_side, double tolerance )
{
  if ( not IsLeaf() )
  {
    // call recursively CreatePoint for the child where the Point shall lie in
    Teuchos::RCP<Point> p = Leaf( x )->CreatePoint( newid, x, cut_edge, cut_side, tolerance );
    // add the pointer not only in the leaf but also on the current level
    AddPoint( x, p );
    return p;
  }
  else
  {
    // create a new point and add the point at the lowest level
    Teuchos::RCP<Point> p = GEO::CUT::CreatePoint( newid, x, cut_edge, cut_side, tolerance );
    AddPoint( x, p );
    return p;
  }
}


/*-----------------------------------------------------------------------------------------*
 * Simply insert p into the pointpool and correspondingly modify the boundingbox size
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::AddPoint( const double * x, Teuchos::RCP<Point> p )
{
  points_.insert( p ); // insert the point in the pointpool
  bb_->AddPoint( x );   // modify the boundingbox size
}


/*-----------------------------------------------------------------------------------------*
 * get the leaf where the point with the given coordinates lies in
 *-----------------------------------------------------------------------------------------*/
GEO::CUT::OctTreeNode* GEO::CUT::OctTreeNode::Leaf( const double * x )
{
  // navigate to the right one of the 8 children nodes
  //
  //    z <0            1  |  3         z > 0        5  |  7
  //        ____ y     ____|____                    ____|____
  //       |               |                            |
  //       |            2  |  4                      6  |  8
  //       x
  //

  int idx = 0;
  if ( x[0] > splitpoint_( 0 ) ) // add an index of one to move in x direction
  {
    idx += 1;
  }
  if ( x[1] > splitpoint_( 1 ) ) // add an index of two to move in y direction
  {
    idx += 2;
  }
  if ( x[2] > splitpoint_( 2 ) ) // add an index of four to move in z direction
  {
    idx += 4;
  }
  return &*nodes_[idx];
}


/*-----------------------------------------------------------------------------------------*
 * split the current boounding box (tree-node)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::Split( int level )
{
  // We must not end up with a OctTreeNode that holds just nodes from the
  // cutter mesh. However, there is no real way to test this right now.

  if ( points_.size()>125 ) /// 125 = 1/8 *1000 -> see NewPoint
  {
    LINALG::Matrix<3,1> x;
    bool first = true;

    for ( RCPPointSet::iterator i=points_.begin(); i!=points_.end(); ++i )
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
      nodes_[i] = Teuchos::rcp( new OctTreeNode( norm_ ) );
    }

    // avoid empty room (room not covered by boundary boxes)
    for ( int i=0; i<8; ++i )
    {
      // always have the split point in all boxes
      nodes_[i]->bb_->AddPoint( splitpoint_ );

      // always have the outmost point in each box
      double x[3];
      bb_->CornerPoint( i, x );
      Leaf( x )->bb_->AddPoint( x );
    }

    for ( RCPPointSet::iterator i=points_.begin(); i!=points_.end(); ++i )
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


/*-----------------------------------------------------------------------------------------*
 * collect all edges
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::CollectEdges( const BoundingBox & edgebox, plain_edge_set & edges )
{
  if ( not IsLeaf() )
  {
    if ( edgebox.Within( norm_, *bb_ ) )
    {
      for ( int i=0; i<8; ++i )
      {
        nodes_[i]->CollectEdges( edgebox, edges );
      }
    }
  }
  else
  {
    Teuchos::RCP<BoundingBox> sbox = Teuchos::rcp( BoundingBox::Create() );
    for ( RCPPointSet::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;
      const plain_edge_set & sds = p->CutEdges();
      for ( plain_edge_set::const_iterator i=sds.begin(); i!=sds.end(); ++i )
      {
        Edge * s = *i;
        if ( edges.count( s )==0 )
        {
          sbox->Assign( *s );
          if ( sbox->Within( norm_, edgebox ) )
          {
            edges.insert( s );
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * collect all sides
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::CollectSides( const BoundingBox & sidebox, plain_side_set & sides )
{
  if ( not IsLeaf() )
  {
    if ( sidebox.Within( norm_, *bb_ ) )
    {
      for ( int i=0; i<8; ++i )
      {
        nodes_[i]->CollectSides( sidebox, sides );
      }
    }
  }
  else
  {
    Teuchos::RCP<BoundingBox> sbox = Teuchos::rcp( BoundingBox::Create() );
    for ( RCPPointSet::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;
      const plain_side_set & sds = p->CutSides();
      for ( plain_side_set::const_iterator i=sds.begin(); i!=sds.end(); ++i )
      {
        Side * s = *i;
        if ( sides.count( s )==0 )
        {
          sbox->Assign( *s );
          if ( sbox->Within( norm_, sidebox ) )
          {
            sides.insert( s );
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * collect all elements near the sidebox
 * (unused, does not work properly when there is no point adjacent to elements in a tree's leaf,
 * e.g. when side lies within an element)
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::CollectElements( const BoundingBox & sidebox, plain_element_set & elements )
{
  // see REMARK in cut_mesh.cpp
  dserror("collecting elements via the OctTreeNode does not find all possible element-side intersections");

  if ( not IsLeaf() )
  {
    if ( sidebox.Within( norm_, *bb_ ) ) // within check is a check of overlap between the 2 bounding boxes
    {
      for ( int i=0; i<8; ++i )
      {
        nodes_[i]->CollectElements( sidebox, elements );
      }
    }
  }
  else
  {
    Teuchos::RCP<BoundingBox> elementbox = Teuchos::rcp( BoundingBox::Create() );
    for ( RCPPointSet::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;
      const plain_element_set & els = p->Elements();

      // add all elements adjacent to the current point
      // REMARK: this does not find all elements that have an overlap with the sidebox!!!
      for ( plain_element_set::const_iterator i=els.begin(); i!=els.end(); ++i )
      {
        Element * e = *i;
        if ( elements.count( e )==0 )
        {
          elementbox->Assign( *e );
          if ( elementbox->Within( norm_, sidebox ) )
          {
            elements.insert( e );
          }
        }
      }
    }
  }
}


/*-----------------------------------------------------------------------------------------*
 * reset the Point::Position of outside points
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::ResetOutsidePoints()
{
  for ( RCPPointSet::iterator i=points_.begin();
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


/*-----------------------------------------------------------------------------------------*
 * print the tree at a given level
 *-----------------------------------------------------------------------------------------*/
void GEO::CUT::OctTreeNode::Print( int level, std::ostream & stream )
{
  if ( not IsLeaf() )
  {
    for ( int i=0; i<8; ++i )
    {
      nodes_[i]->Print( level+1, stream );
    }
  }
  else
  {
    for ( RCPPointSet::iterator i=points_.begin(); i!=points_.end(); ++i )
    {
      Point * p = &**i;
      p->Plot( stream );
    }
    stream << "\n";
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
GEO::CUT::PointPool::PointPool( double norm )
    : tree_( norm ),
      probdim_( DRT::Problem::Instance()->NDim() )
{ }
