#ifndef CUT_ELEMENTHANDLE_H
#define CUT_ELEMENTHANDLE_H

#include <vector>

#include "cut_element.H"

namespace GEO
{
namespace CUT
{
  class Node;
  class Element;
  class Mesh;

class ElementHandle
{
public:

  virtual ~ElementHandle() {}

  virtual bool IsCut() = 0;

  virtual void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells ) = 0;

  virtual void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells ) = 0;

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;
};

class LinearElementHandle : public ElementHandle
{
public:

  LinearElementHandle()
    : element_( NULL )
  {
  }

  explicit LinearElementHandle( Element * e )
    : element_( e )
  {
  }

  virtual bool IsCut() { return element_->IsCut(); }

  virtual void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells )
  {
    element_->GetIntegrationCells( cells );
  }

  virtual void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells )
  {
    element_->GetBoundaryCells( bcells );
  }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst )
  {
    element_->LocalCoordinates( xyz, rst );
  }

private:
  Element * element_;
};

class QuadraticElementHandle : public ElementHandle
{
public:

  virtual bool IsCut();

  virtual void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells );

  virtual void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells );

protected:
  std::vector<Element*> subelements_;
  std::vector<Node*> nodes_;
};

class Hex20ElementHandle : public QuadraticElementHandle
{
public:

  Hex20ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

class Hex27ElementHandle : public QuadraticElementHandle
{
public:

  Hex27ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

class Tet10ElementHandle : public QuadraticElementHandle
{
public:

  Tet10ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

}
}

#endif
