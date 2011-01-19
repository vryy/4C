#ifndef CUT_ELEMENT_H
#define CUT_ELEMENT_H

#include "cut_node.H"
#include "cut_edge.H"
#include "cut_side.H"
#include "cut_tetgen.H"

namespace GEO
{
  namespace CELL
  {
    class Cell;
  }

namespace CUT
{
  class VolumeCell;
  class IntegrationCell;
  class BoundaryCell;

class Element
{
public:

  Element( int eid,
           const std::vector<Side*> & sides,
           const std::vector<Node*> & nodes )
    : eid_( eid ),
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

  virtual bool PointInside( Point* p ) = 0;

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;

  bool Cut( Mesh & mesh, Side & side );

  void MakeFacets( Mesh & mesh );

  void MakeVolumeCells( Mesh & mesh );

  void Delaunay( Mesh & mesh, std::vector<std::vector<Point*> > & tets );

  void FindNodePositions();

  bool IsCut();

#ifdef QHULL
  virtual void FillTetgen( tetgenio & out )
  { throw std::runtime_error( "not supported" ); }
#endif

  bool OnSide( const std::vector<Point*> & facet_points );

  void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells );

  void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells );

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

  void GetCutPoints( std::set<Point*> & cut_points );

  void GetAllPoints( std::set<Point*> & cut_points );

  const std::set<VolumeCell*> & VolumeCells() const { return cells_; }

protected:

  //const std::vector<Side*> & Sides() { return sides_; }

private:

  bool FindCutPoints( Mesh & mesh, Side & side, Side & other );

  bool FindCutLines( Mesh & mesh, Side & side, Side & other );

  int eid_;

  std::vector<Side*> sides_;

  std::vector<Node*> nodes_;

  std::set<Side*> cut_faces_;

  std::set<Facet*> facets_;

  //Teuchos::RCP<GEO::CELL::Cell> cell_;
  std::set<VolumeCell*> cells_;
};


template <>
class ConcreteElement<DRT::Element::tet4> : public Element
{
public:

  ConcreteElement<DRT::Element::tet4>( int eid,
                                       const std::vector<Side*> & sides,
                                       const std::vector<Node*> & nodes )
    : Element( eid, sides, nodes )
  {}

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 4> & xyze )
  { Element::Coordinates( xyze.A() ); }

#ifdef QHULL
  virtual void FillTetgen( tetgenio & out );
#endif

};


template <>
class ConcreteElement<DRT::Element::wedge6> : public Element
{
public:

  ConcreteElement<DRT::Element::wedge6>( int eid,
                                         const std::vector<Side*> & sides,
                                         const std::vector<Node*> & nodes )
    : Element( eid, sides, nodes )
  {}

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 6> & xyze )
  { Element::Coordinates( xyze.A() ); }

};


template <>
class ConcreteElement<DRT::Element::pyramid5> : public Element
{
public:

  ConcreteElement<DRT::Element::pyramid5>( int eid,
                                           const std::vector<Side*> & sides,
                                           const std::vector<Node*> & nodes )
    : Element( eid, sides, nodes )
  {}

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 5> & xyze )
  { Element::Coordinates( xyze.A() ); }

};


template <>
class ConcreteElement<DRT::Element::hex8> : public Element
{
public:

  ConcreteElement<DRT::Element::hex8>( int eid,
                                       const std::vector<Side*> & sides,
                                       const std::vector<Node*> & nodes )
    : Element( eid, sides, nodes )
  {}

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 8> & xyze )
  { Element::Coordinates( xyze.A() ); }

};


}
}

#endif
