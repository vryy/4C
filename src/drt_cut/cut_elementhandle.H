#ifndef CUT_ELEMENTHANDLE_H
#define CUT_ELEMENTHANDLE_H

#include <vector>
#include <map>

#include "../linalg/linalg_fixedsizematrix.H"

#include "../../src/drt_fem_general/drt_utils_gausspoints.H"

#include "cut_element.H"

namespace GEO
{
namespace CUT
{
  class Point;
  class Node;
  class Element;
  class Mesh;

class ElementHandle
{
public:

  virtual ~ElementHandle() {}

  virtual bool IsCut() = 0;

  virtual DRT::Element::DiscretizationType Shape() = 0;

  virtual void GetVolumeCells( std::set<GEO::CUT::VolumeCell*> & cells ) = 0;

  virtual void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells ) = 0;

  virtual void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells ) = 0;

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst ) = 0;

  virtual void VolumeCells( std::set<VolumeCell*> & cells ) const = 0;

  virtual const std::vector<Node*> & Nodes() = 0;

  const LINALG::Matrix<3,1> & LocalCoordinates( Point * p )
  {
    std::map<Point*, LINALG::Matrix<3,1> >::iterator i = local_coordinates_.find( p );
    if ( i!=local_coordinates_.end() )
    {
      return i->second;
    }
    LINALG::Matrix<3,1> & rst = local_coordinates_[p];
    LINALG::Matrix<3,1> xyz;
    p->Coordinates( xyz.A() );
    LocalCoordinates( xyz, rst );
    return rst;
  }

  void VolumeCellGaussPoints( std::vector<DRT::UTILS::GaussIntegration> & intpoints );

private:

  template <DRT::Element::DiscretizationType distype>
  Teuchos::RCP<DRT::UTILS::GaussPoints> CreateProjected( GEO::CUT::IntegrationCell * ic );

  std::map<Point*, LINALG::Matrix<3,1> > local_coordinates_;
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

  virtual DRT::Element::DiscretizationType Shape() { return element_->Shape(); }

  virtual void GetVolumeCells( std::set<GEO::CUT::VolumeCell*> & cells );

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

  virtual void VolumeCells( std::set<VolumeCell*> & cells ) const
  {
    cells = element_->VolumeCells();
  }

  virtual const std::vector<Node*> & Nodes()
  {
    return element_->Nodes();
  }

private:
  Element * element_;
};

class QuadraticElementHandle : public ElementHandle
{
public:

  virtual bool IsCut();

  virtual void GetVolumeCells( std::set<GEO::CUT::VolumeCell*> & cells );

  virtual void GetIntegrationCells( std::set<GEO::CUT::IntegrationCell*> & cells );

  virtual void GetBoundaryCells( std::set<GEO::CUT::BoundaryCell*> & bcells );

  virtual void VolumeCells( std::set<VolumeCell*> & cells ) const;

  virtual const std::vector<Node*> & Nodes()
  {
    return nodes_;
  }

protected:
  std::vector<Element*> subelements_;
  std::vector<Node*> nodes_;
};

class Hex20ElementHandle : public QuadraticElementHandle
{
public:

  Hex20ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::hex20; }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

class Hex27ElementHandle : public QuadraticElementHandle
{
public:

  Hex27ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::hex27; }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

class Tet10ElementHandle : public QuadraticElementHandle
{
public:

  Tet10ElementHandle( Mesh & mesh, int eid, const std::vector<int> & nids );

  virtual DRT::Element::DiscretizationType Shape() { return DRT::Element::tet10; }

  virtual void LocalCoordinates( const LINALG::Matrix<3,1> & xyz, LINALG::Matrix<3,1> & rst );

};

}
}

#endif
