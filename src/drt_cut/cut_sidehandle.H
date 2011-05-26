#ifndef CUT_SIDEHANDLE_H
#define CUT_SIDEHANDLE_H

#include <vector>

#include "cut_side.H"
#include "cut_boundarycell.H"
//#include "cut_elementhandle.H"

#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../../src/drt_fem_general/drt_utils_gausspoints.H"

namespace GEO
{
namespace CUT
{
  class Node;
  class Side;
  class Mesh;

class SideHandle
{
public:

  virtual ~SideHandle() {}

  virtual DRT::Element::DiscretizationType Shape() = 0;

  virtual void Coordinates( Epetra_SerialDenseMatrix & xyze ) = 0;

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rst ) = 0;

  virtual void CollectSides( std::set<Side*> & sides ) = 0;

  template <DRT::Element::DiscretizationType distype>
  Teuchos::RCP<DRT::UTILS::GaussPoints> CreateProjected( BoundaryCell * bc )
  {
    const unsigned nen = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;

//     LINALG::Matrix<3, nen> xyze;
    LINALG::Matrix<2, nen> xie;

    const std::vector<GEO::CUT::Point*> & cpoints = bc->Points();
    if ( cpoints.size() != nen )
      throw std::runtime_error( "non-matching number of points" );

    for ( unsigned i=0; i<nen; ++i )
    {
      GEO::CUT::Point * p = cpoints[i];
      const LINALG::Matrix<2,1> & xi = LocalCoordinates( p );
      std::copy( xi.A(), xi.A()+2, &xie( 0, i ) );
//       const LINALG::Matrix<3,1> & eta = e->LocalCoordinates( p );
//       std::copy( eta.A(), eta.A()+3, &xyze( 0, i ) );
    }

    Teuchos::RCP<DRT::UTILS::GaussPoints> gp =
      DRT::UTILS::GaussIntegration::CreateProjected<distype>( xie, bc->CubatureDegree( Shape() ) );
    return gp;
  }

  const LINALG::Matrix<2,1> & LocalCoordinates( Point * p )
  {
    std::map<Point*, LINALG::Matrix<2,1> >::iterator i = local_coordinates_.find( p );
    if ( i!=local_coordinates_.end() )
    {
      return i->second;
    }
    LINALG::Matrix<2,1> & rst = local_coordinates_[p];
    LINALG::Matrix<3,1> xyz;
    p->Coordinates( xyz.A() );
    LocalCoordinates( xyz, rst );
    return rst;
  }

private:
  std::map<Point*, LINALG::Matrix<2,1> > local_coordinates_;
};

class LinearSideHandle : public SideHandle
{
public:

  LinearSideHandle()
    : side_( NULL )
  {
  }

  explicit LinearSideHandle( Side * s )
    : side_( s )
  {
  }

  virtual DRT::Element::DiscretizationType Shape()
  {
    return side_->Shape();
  }

  virtual void Coordinates( Epetra_SerialDenseMatrix & xyze )
  {
    xyze.Reshape( 3, side_->Nodes().size() );
    side_->Coordinates( xyze.A() );
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rs )
  {
    LINALG::Matrix<3,1> rst;
    side_->LocalCoordinates( xyz, rst );
    rs( 0 ) = rst( 0 );
    rs( 1 ) = rst( 1 );
  }

  virtual void CollectSides( std::set<Side*> & sides )
  {
    sides.insert( side_ );
  }

private:
  Side * side_;
};

class QuadraticSideHandle : public SideHandle
{
public:

  virtual void Coordinates( Epetra_SerialDenseMatrix & xyze )
  {
    xyze.Reshape( 3, nodes_.size() );
    for ( std::vector<Node*>::iterator i=nodes_.begin();
          i!=nodes_.end();
          ++i )
    {
      Node * n = *i;
      n->Coordinates( &xyze( 0, i-nodes_.begin() ) );
    }
  }

  virtual void CollectSides( std::set<Side*> & sides )
  {
    std::copy( subsides_.begin(), subsides_.end(), std::inserter( sides, sides.begin() ) );
  }

protected:
  std::vector<Side*> subsides_;
  std::vector<Node*> nodes_;
};

class Tri6SideHandle : public QuadraticSideHandle
{
public:

  Tri6SideHandle( Mesh & mesh, int sid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::tri6;
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rst );
};

class Quad4SideHandle : public QuadraticSideHandle
{
public:

  Quad4SideHandle( Mesh & mesh, int sid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::quad4;
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rst );
};

class Quad8SideHandle : public QuadraticSideHandle
{
public:

  Quad8SideHandle( Mesh & mesh, int sid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::quad8;
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rst );
};

class Quad9SideHandle : public QuadraticSideHandle
{
public:

  Quad9SideHandle( Mesh & mesh, int sid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape()
  {
    return DRT::Element::quad9;
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<2,1> & rst );
};

}
}

#endif
