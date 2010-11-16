#ifndef CUT_ELEMENT_H
#define CUT_ELEMENT_H

#include "cut_node.H"
#include "cut_edge.H"
#include "cut_side.H"
#include "cut_tetgen.H"

namespace GEO
{
namespace CUT
{
  class CellGenerator;

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
  }

  virtual ~Element() {}

  int Id() const { return eid_; }

  virtual bool PointInside( Point* p ) = 0;

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;

  virtual void FillComplete( Mesh & mesh );

  virtual void Cut( Mesh & mesh, LinearSide & side ) = 0;

  virtual void MakeFacets( Mesh & mesh ) = 0;

  virtual void FindNodePositions() = 0;

  virtual bool IsCut() = 0;

  virtual void GenerateTetgen( Mesh & mesh, CellGenerator * generator ) = 0;

#ifdef QHULL
  virtual void FillTetgen( tetgenio & out )
  { throw std::runtime_error( "not supported" ); }
#endif

  virtual bool OnSide( const std::vector<Point*> & facet_points ) = 0;

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

protected:

  //const std::vector<Side*> & Sides() { return sides_; }

private:

  int eid_;

  std::vector<Side*> sides_;

  std::vector<Node*> nodes_;
};


class LinearElement : public Element
{
public:

  LinearElement( int eid,
                 const std::vector<Side*> & sides,
                 const std::vector<Node*> & nodes )
    : Element( eid, sides, nodes )
  {
    for ( std::vector<Node*>::const_iterator i=nodes.begin(); i!=nodes.end(); ++i )
    {
      Node * n = *i;
      n->point()->AddElement( this );
    }
  }

  virtual void Cut( Mesh & mesh, LinearSide & side );

  virtual void MakeFacets( Mesh & mesh );

  virtual void FindNodePositions();

  virtual bool IsCut();

  virtual void GenerateTetgen( Mesh & mesh, CellGenerator * generator );

  void GenerateTetgen( Mesh & mesh, Element * parent, CellGenerator * generator );

  virtual bool OnSide( const std::vector<Point*> & facet_points );

private:

  bool FindCutPoints( Mesh & mesh, LinearSide & side, LinearSide & other );

  bool FindCutLines( Mesh & mesh, LinearSide & side, LinearSide & other );

  std::set<Side*> cut_faces_;

  std::set<Facet*> facets_;
};


class QuadraticElement : public Element
{
public:

  QuadraticElement( int eid,
                    const std::vector<Side*> & sides,
                    const std::vector<Node*> & nodes )
    : Element( eid, sides, nodes )
  {}

  virtual void Cut( Mesh & mesh, LinearSide & side );

  virtual void MakeFacets( Mesh & mesh );

  virtual void FindNodePositions();

  virtual bool IsCut();

  virtual void GenerateTetgen( Mesh & mesh, CellGenerator * generator );

  virtual bool OnSide( const std::vector<Point*> & facet_points );

protected:

  std::vector<Element*> subelements_;
};


template <>
class ConcreteElement<DRT::Element::tet4> : public LinearElement
{
public:

  ConcreteElement<DRT::Element::tet4>( int eid,
                                       const std::vector<Side*> & sides,
                                       const std::vector<Node*> & nodes )
    : LinearElement( eid, sides, nodes )
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
class ConcreteElement<DRT::Element::wedge6> : public LinearElement
{
public:

  ConcreteElement<DRT::Element::wedge6>( int eid,
                                         const std::vector<Side*> & sides,
                                         const std::vector<Node*> & nodes )
    : LinearElement( eid, sides, nodes )
  {}

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 6> & xyze )
  { Element::Coordinates( xyze.A() ); }

};


template <>
class ConcreteElement<DRT::Element::pyramid5> : public LinearElement
{
public:

  ConcreteElement<DRT::Element::pyramid5>( int eid,
                                           const std::vector<Side*> & sides,
                                           const std::vector<Node*> & nodes )
    : LinearElement( eid, sides, nodes )
  {}

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 5> & xyze )
  { Element::Coordinates( xyze.A() ); }

};


template <>
class ConcreteElement<DRT::Element::hex8> : public LinearElement
{
public:

  ConcreteElement<DRT::Element::hex8>( int eid,
                                       const std::vector<Side*> & sides,
                                       const std::vector<Node*> & nodes )
    : LinearElement( eid, sides, nodes )
  {}

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 8> & xyze )
  { Element::Coordinates( xyze.A() ); }

};


template <>
class ConcreteElement<DRT::Element::hex20> : public QuadraticElement
{
public:

  ConcreteElement<DRT::Element::hex20>( int eid,
                                        const std::vector<Side*> & sides,
                                        const std::vector<Node*> & nodes )
    : QuadraticElement( eid, sides, nodes )
  {}

  virtual void FillComplete( Mesh & mesh );

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 20> & xyze )
  { Element::Coordinates( xyze.A() ); }

private:

};

template <>
class ConcreteElement<DRT::Element::hex27> : public QuadraticElement
{
public:

  ConcreteElement<DRT::Element::hex27>( int eid,
                                        const std::vector<Side*> & sides,
                                        const std::vector<Node*> & nodes )
    : QuadraticElement( eid, sides, nodes )
  {}

  virtual void FillComplete( Mesh & mesh );

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 27> & xyze )
  { Element::Coordinates( xyze.A() ); }

private:

};

template <>
class ConcreteElement<DRT::Element::tet10> : public QuadraticElement
{
public:

  ConcreteElement<DRT::Element::tet10>( int eid,
                                        const std::vector<Side*> & sides,
                                        const std::vector<Node*> & nodes )
    : QuadraticElement( eid, sides, nodes )
  {}

  virtual void FillComplete( Mesh & mesh );

  virtual bool PointInside( Point* p );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

  void Coordinates( LINALG::Matrix<3, 10> & xyze )
  { Element::Coordinates( xyze.A() ); }

private:

};

#if 0

template <>
class ConcreteElement<DRT::Element::wedge15> : public QuadraticElement
{
public:

  ConcreteElement<DRT::Element::wedge15>( int eid,
                                          const std::vector<Side*> & sides,
                                          const std::vector<Node*> & nodes )
    : QuadraticElement( eid, sides, nodes )
  {}

};

#endif


}
}

#endif
