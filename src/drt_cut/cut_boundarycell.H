#ifndef CUT_BOUNDARYCELL_H
#define CUT_BOUNDARYCELL_H

#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cmath>

#include <Epetra_SerialDenseMatrix.h>

#include "../drt_lib/drt_element.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"

namespace GEO
{
namespace CUT
{
  class Facet;
  class VolumeCell;
  class Mesh;
  class Point;

class BoundaryCell
{
public:

  BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet, const std::vector<Point*> & points )
    : xyz_( xyz ),
      facet_( facet ),
      points_( points )
  {
  }

  virtual ~BoundaryCell() {}

  virtual DRT::Element::DiscretizationType Shape() const = 0;

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const = 0;

  virtual void DumpGmsh( std::ofstream & file ) = 0;

  virtual double Area() = 0;

  virtual void Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal ) const = 0;

  const std::vector<Point*> & Points() const { return points_; }

  const Epetra_SerialDenseMatrix & Coordinates() const { return xyz_; }

  Facet * GetFacet() { return facet_; }

protected:

  virtual DRT::UTILS::GaussRule2D MySimpleGaussRule() = 0;

  template <DRT::Element::DiscretizationType distype>
  double MyArea()
  {
    const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<distype>::numNodePerElement;
    LINALG::Matrix<3,numnodes> xyze( this->xyz_, true );
    LINALG::Matrix<numnodes,1> funct;
    LINALG::Matrix<2,numnodes> deriv;
    LINALG::Matrix<2,2> metrictensor;

    DRT::UTILS::GaussRule2D gaussrule = this->MySimpleGaussRule();
    DRT::UTILS::IntegrationPoints2D intpoints( gaussrule );

    double area = 0;
    double drs;
    for ( int i=0; i<intpoints.nquad; ++i )
    {
      double * eta = intpoints.qxg[i];
      DRT::UTILS::shape_function_2D( funct, eta[0], eta[1], distype );
      DRT::UTILS::shape_function_2D_deriv1( deriv, eta[0], eta[1], distype );
      DRT::UTILS::ComputeMetricTensorForBoundaryEle<distype>( xyze, deriv, metrictensor, drs, NULL );
      if ( not std::isnan( drs ) )
        area += intpoints.qwgt[i] * drs;
    }
    return area;
  }

  Epetra_SerialDenseMatrix xyz_;
  Facet * facet_;
  std::vector<Point*> points_;
};

class Tri3BoundaryCell : public BoundaryCell
{
public:

  Tri3BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet, const std::vector<Point*> & points )
    : BoundaryCell( xyz, facet, points )
  {
  }

//   static void CreateCell( Mesh & mesh, VolumeCell * cell, Facet * facet, const std::vector<Point*> & side );

//   static void CollectCoordinates( const std::vector<Point*> & side,
//                                   std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz );

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tri3; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const { return 20; }

  virtual void DumpGmsh( std::ofstream & file );

  virtual double Area() { return MyArea<DRT::Element::tri3>(); }

  virtual void Normal( const LINALG::Matrix<2,1> & pos, LINALG::Matrix<3,1> & normal ) const;

protected:

  virtual DRT::UTILS::GaussRule2D MySimpleGaussRule() { return DRT::UTILS::intrule_tri_3point; }
};

class Quad4BoundaryCell : public BoundaryCell
{
public:

  Quad4BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet, const std::vector<Point*> & points )
    : BoundaryCell( xyz, facet, points )
  {
  }

//   static void CreateCell( Mesh & mesh, VolumeCell * cell, Facet * facet, const std::vector<Point*> & side );

//   static void CollectCoordinates( const std::vector<Point*> & side,
//                                   std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz );

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::quad4; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const { return 20; }

  virtual void DumpGmsh( std::ofstream & file );

  virtual double Area() { return MyArea<DRT::Element::quad4>(); }

  virtual void Normal( const LINALG::Matrix<2,1> & pos, LINALG::Matrix<3,1> & normal ) const;

protected:

  virtual DRT::UTILS::GaussRule2D MySimpleGaussRule() { return DRT::UTILS::intrule_quad_4point; }
};

}
}

#endif
