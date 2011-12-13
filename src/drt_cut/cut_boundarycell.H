#ifndef CUT_BOUNDARYCELL_H
#define CUT_BOUNDARYCELL_H

#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cmath>

#include <Epetra_SerialDenseMatrix.h>

#include "cut_cycle.H"

#include "../drt_lib/drt_element.H"

#include "../drt_fem_general/drt_utils_integration.H"
#include "../drt_fem_general/drt_utils_fem_shapefunctions.H"
#include "../drt_fem_general/drt_utils_boundary_integration.H"
#include "../../src/drt_fem_general/drt_utils_gausspoints.H"

namespace GEO
{
namespace CUT
{
  class Facet;
  class VolumeCell;
  class Mesh;
  class Point;

/// Base of boundary cells.
/*!
 * Boundary cells are used to represent the cut surface. Each volume cell has
 * its own boundary cells at any cut surface with outward normals.
 */
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

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL ) = 0;

  virtual double Area() = 0;

  virtual void Normal( const LINALG::Matrix<2,1> & xsi, LINALG::Matrix<3,1> & normal ) const = 0;

  const Cycle & PointCycle() const { return points_; }

  const std::vector<Point*> & Points() const { return points_(); }

  const Epetra_SerialDenseMatrix & Coordinates() const { return xyz_; }

  Facet * GetFacet() { return facet_; }

  void Clear() { points_.clear(); }

  bool IsValid() const { return points_.size() > 0; }

  virtual DRT::UTILS::GaussIntegration gaussRule()=0;


  template <DRT::Element::DiscretizationType celldistype>
  void Transform( const LINALG::Matrix<2,1> & eta,
                  LINALG::Matrix<3,1>       & x_gp_lin,
                  LINALG::Matrix<3,1>       & normal,
                  double                    & drs)
  {
      const int numnodes = DRT::UTILS::DisTypeToNumNodePerEle<celldistype>::numNodePerElement;
      LINALG::Matrix<3,numnodes> xyze( this->xyz_, true );
      LINALG::Matrix<numnodes,1> funct;
      LINALG::Matrix<2,numnodes> deriv;
      LINALG::Matrix<2,2> metrictensor;

      DRT::UTILS::shape_function_2D( funct, eta(0), eta(1), celldistype );
      DRT::UTILS::shape_function_2D_deriv1( deriv, eta(0), eta(1), celldistype );
      DRT::UTILS::ComputeMetricTensorForBoundaryEle<celldistype>( xyze, deriv, metrictensor, drs, &normal );

      x_gp_lin.Multiply(xyze, funct);

      return;
  }


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
  Cycle points_;
};

/// tri3 boundary cell
class Tri3BoundaryCell : public BoundaryCell
{
public:

  Tri3BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet, const std::vector<Point*> & points )
    : BoundaryCell( xyz, facet, points )
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tri3; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const { return 20; }

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL );

  virtual double Area() { return MyArea<DRT::Element::tri3>(); }

  virtual void Normal( const LINALG::Matrix<2,1> & pos, LINALG::Matrix<3,1> & normal ) const;

  virtual DRT::UTILS::GaussIntegration gaussRule();

protected:

  virtual DRT::UTILS::GaussRule2D MySimpleGaussRule() { return DRT::UTILS::intrule_tri_3point; }
};

/// quad4 boundary cell
class Quad4BoundaryCell : public BoundaryCell
{
public:

  Quad4BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet, const std::vector<Point*> & points )
    : BoundaryCell( xyz, facet, points )
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::quad4; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const { return 20; }

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL );

  virtual double Area() { return MyArea<DRT::Element::quad4>(); }

  virtual void Normal( const LINALG::Matrix<2,1> & pos, LINALG::Matrix<3,1> & normal ) const;

  virtual DRT::UTILS::GaussIntegration gaussRule();

protected:

  virtual DRT::UTILS::GaussRule2D MySimpleGaussRule() { return DRT::UTILS::intrule_quad_4point; }
};

/// Irregular boundary cell generated during the cut
class ArbitraryBoundaryCell : public BoundaryCell
{
public:

  ArbitraryBoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet, const std::vector<Point*> & points,
      DRT::UTILS::GaussIntegration gaussRule)
    : BoundaryCell( xyz, facet, points ),
      gaussRule_(gaussRule)
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::dis_none; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const {  return 0; }

  virtual void DumpGmsh( std::ofstream & file, int * value=NULL );

  virtual double Area() { return 0.0;}

  virtual void Normal( const LINALG::Matrix<2,1> & pos, LINALG::Matrix<3,1> & normal ) const;

  virtual DRT::UTILS::GaussIntegration gaussRule();

protected:

  virtual DRT::UTILS::GaussRule2D MySimpleGaussRule() { return DRT::UTILS::intrule_quad_4point; }
private:
  DRT::UTILS::GaussIntegration gaussRule_;

};

}
}

#endif
