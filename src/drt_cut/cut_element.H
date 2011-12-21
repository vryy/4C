#ifndef CUT_ELEMENT_H
#define CUT_ELEMENT_H

#include "cut_node.H"
#include "cut_edge.H"
#include "cut_side.H"

namespace GEO
{

namespace CUT
{
  class VolumeCell;
  class IntegrationCell;
  class BoundaryCell;
  class Creator;

/// Base of all elements.
class Element
{
public:

  Element( int eid,
           const std::vector<Side*> & sides,
           const std::vector<Node*> & nodes,
           bool active )
    : eid_( eid ),
      active_( active ),
      sides_( sides ),
      nodes_( nodes )
  {
    for ( std::vector<Side*>::const_iterator i=sides.begin();
          i!=sides.end();
          ++i )
    {
      Side * s = *i;
      s->Register( this );
    }
    for ( std::vector<Node*>::const_iterator i=nodes.begin();
          i!=nodes.end();
          ++i )
    {
      Node * n = *i;
      n->Register( this );
    }
  }

  virtual ~Element() {}

  int Id() const { return eid_; }

  virtual DRT::Element::DiscretizationType Shape() const = 0;

  virtual bool PointInside( Point* p ) = 0;

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz ) = 0;

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst ) = 0;

  bool Cut( Mesh & mesh, Side & side, int recursion );

  void MakeCutLines( Mesh & mesh, Creator & creator );

  void MakeFacets( Mesh & mesh );

  void MakeVolumeCells( Mesh & mesh );

  void CreateIntegrationCells( Mesh & mesh, int count, bool levelset );

  void MomentFitGaussWeights( Mesh & mesh, bool include_inner, std::string Bcellgausstype );
  /// After the integration and boundary cells have been created, there might
  /// still be some empty volume cells. Those are artificial cells with zero
  /// volume and need to be removed from their elements, before we count nodal
  /// dof sets.
  void RemoveEmptyVolumeCells();

  void FindNodePositions();

  bool IsCut();

  bool OnSide( Facet * f );

  bool OnSide( const std::vector<Point*> & facet_points );

  void GetIntegrationCells( plain_integrationcell_set & cells );

  void GetBoundaryCells( plain_boundarycell_set & bcells );

  const std::vector<Side*> & Sides() { return sides_; }

  bool OwnedSide( Side * side )
  {
    return std::find( sides_.begin(), sides_.end(), side )!=sides_.end();
  }

  void Coordinates( double * xyze )
  {
    double * x = xyze;
    for ( std::vector<Node*>::iterator i=nodes_.begin(); i!=nodes_.end(); ++i )
    {
      Node & n = **i;
      n.Coordinates( x );
      x += 3;
    }
  }

  const std::vector<Node*> & Nodes() { return nodes_; }

  const plain_side_set & CutSides() const { return cut_faces_; }

  void GetCutPoints( PointSet & cut_points );

  const plain_volumecell_set & VolumeCells() const { return cells_; }

  const plain_facet_set & Facets() const { return facets_; }

  int NumGaussPoints( DRT::Element::DiscretizationType shape );

  void DebugDump();

  void GnuplotDump();

  void DumpFacets();

  void AssignOtherVolumeCell( VolumeCell * vc ) { cells_.insert( vc ); }

protected:

  //const std::vector<Side*> & Sides() { return sides_; }

private:

  bool FindCutPoints( Mesh & mesh, Side & side, Side & other, int recursion );

  bool FindCutLines( Mesh & mesh, Side & side, Side & other );

  int eid_;

  bool active_;

  std::vector<Side*> sides_;

  std::vector<Node*> nodes_;

  plain_side_set cut_faces_;

  plain_facet_set facets_;

  //Teuchos::RCP<GEO::CELL::Cell> cell_;
  plain_volumecell_set cells_;
};

/// tet4 element
template <>
class ConcreteElement<DRT::Element::tet4> : public Element
{
public:

  ConcreteElement<DRT::Element::tet4>( int eid,
                                       const std::vector<Side*> & sides,
                                       const std::vector<Node*> & nodes,
                                       bool active )
    : Element( eid, sides, nodes, active )
  {}

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tet4; }

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::tet4 );

    xyz = 0;

    const std::vector<Node*> & nodes = Nodes();
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x( nodes[i]->point()->X() );
      xyz.Update( funct( i ), x, 1 );
    }
  }

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::tet4>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::tet4 );

    LINALG::Matrix<numNodesElement,1> scalar( &ns[0] );
    LINALG::Matrix<1,1> res;
    res.MultiplyTN( funct, scalar );
    return res( 0 );
  }

  void Coordinates( LINALG::Matrix<3, 4> & xyze )
  { Element::Coordinates( xyze.A() ); }

};

/// wedge6 element
template <>
class ConcreteElement<DRT::Element::wedge6> : public Element
{
public:

  ConcreteElement<DRT::Element::wedge6>( int eid,
                                         const std::vector<Side*> & sides,
                                         const std::vector<Node*> & nodes,
                                         bool active )
    : Element( eid, sides, nodes, active )
  {}

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::wedge6; }

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge6>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::wedge6 );

    xyz = 0;

    const std::vector<Node*> & nodes = Nodes();
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x( nodes[i]->point()->X() );
      xyz.Update( funct( i ), x, 1 );
    }
  }

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::wedge6>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::wedge6 );

    LINALG::Matrix<numNodesElement,1> scalar( &ns[0] );
    LINALG::Matrix<1,1> res;
    res.MultiplyTN( funct, scalar );
    return res( 0 );
  }

  void Coordinates( LINALG::Matrix<3, 6> & xyze )
  { Element::Coordinates( xyze.A() ); }

};

/// pyramid5 element
template <>
class ConcreteElement<DRT::Element::pyramid5> : public Element
{
public:

  ConcreteElement<DRT::Element::pyramid5>( int eid,
                                           const std::vector<Side*> & sides,
                                           const std::vector<Node*> & nodes,
                                           bool active )
    : Element( eid, sides, nodes, active )
  {}

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::pyramid5; }

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::pyramid5>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::pyramid5 );

    xyz = 0;

    const std::vector<Node*> & nodes = Nodes();
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x( nodes[i]->point()->X() );
      xyz.Update( funct( i ), x, 1 );
    }
  }

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::pyramid5>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::pyramid5 );

    LINALG::Matrix<numNodesElement,1> scalar( &ns[0] );
    LINALG::Matrix<1,1> res;
    res.MultiplyTN( funct, scalar );
    return res( 0 );
  }
  void Coordinates( LINALG::Matrix<3, 5> & xyze )
  { Element::Coordinates( xyze.A() ); }

};

/// hex8 element
template <>
class ConcreteElement<DRT::Element::hex8> : public Element
{
public:

  ConcreteElement<DRT::Element::hex8>( int eid,
                                       const std::vector<Side*> & sides,
                                       const std::vector<Node*> & nodes,
                                       bool active )
    : Element( eid, sides, nodes, active )
  {}

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::hex8; }

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  virtual void GlobalCoordinates( const LINALG::Matrix<3,1> & rst, LINALG::Matrix<3,1> & xyz )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );

    xyz = 0;

    const std::vector<Node*> & nodes = Nodes();
    for ( int i=0; i<numNodesElement; ++i )
    {
      LINALG::Matrix<3,1> x( nodes[i]->point()->X() );
      xyz.Update( funct( i ), x, 1 );
    }
  }

  virtual double Scalar( const std::vector<double> & ns, const LINALG::Matrix<3,1> & rst )
  {
    const int numNodesElement = DRT::UTILS::DisTypeToNumNodePerEle<DRT::Element::hex8>::numNodePerElement;
    LINALG::Matrix<numNodesElement,1> funct;
    DRT::UTILS::shape_function_3D( funct, rst( 0 ), rst( 1 ), rst( 2 ), DRT::Element::hex8 );

    LINALG::Matrix<numNodesElement,1> scalar( &ns[0] );
    LINALG::Matrix<1,1> res;
    res.MultiplyTN( funct, scalar );
    return res( 0 );
  }

  void Coordinates( LINALG::Matrix<3, 8> & xyze )
  { Element::Coordinates( xyze.A() ); }

};


}
}

#endif
