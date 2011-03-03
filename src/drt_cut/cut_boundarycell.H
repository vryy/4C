#ifndef CUT_BOUNDARYCELL_H
#define CUT_BOUNDARYCELL_H

#include <fstream>
#include <vector>
#include <set>
#include <map>

#include <Epetra_SerialDenseMatrix.h>

#include "../drt_lib/drt_element.H"

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

  virtual double Area() const = 0;

  const std::vector<Point*> & Points() const { return points_; }

  const Epetra_SerialDenseMatrix & Coordinates() const { return xyz_; }

  Facet * GetFacet() { return facet_; }

protected:
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

  static void CreateCell( Mesh & mesh, VolumeCell * cell, Facet * facet, const std::vector<Point*> & side );

//   static void CollectCoordinates( const std::vector<Point*> & side,
//                                   std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz );

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tri3; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const { return 7; }

  virtual void DumpGmsh( std::ofstream & file );

  virtual double Area() const;
};

class Quad4BoundaryCell : public BoundaryCell
{
public:

  Quad4BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet, const std::vector<Point*> & points )
    : BoundaryCell( xyz, facet, points )
  {
  }

  static void CreateCell( Mesh & mesh, VolumeCell * cell, Facet * facet, const std::vector<Point*> & side );

//   static void CollectCoordinates( const std::vector<Point*> & side,
//                                   std::map<Facet*, std::vector<Epetra_SerialDenseMatrix> > & sides_xyz );

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::quad4; }

  virtual int CubatureDegree( DRT::Element::DiscretizationType elementshape ) const { return 7; }

  virtual void DumpGmsh( std::ofstream & file );

  virtual double Area() const;
};

}
}

#endif
