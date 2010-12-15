#ifndef CUT_BOUNDARYCELL_H
#define CUT_BOUNDARYCELL_H

#include <fstream>
#include <vector>
#include <set>

#include <Epetra_SerialDenseMatrix.h>

#include "../drt_lib/drt_element.H"

namespace GEO
{
namespace CUT
{
  class Facet;

class BoundaryCell
{
public:

  BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet )
    : xyz_( xyz ),
      facet_( facet )
  {
  }

  virtual ~BoundaryCell() {}

  virtual DRT::Element::DiscretizationType Shape() const = 0;

  const Epetra_SerialDenseMatrix & Coordinates() const { return xyz_; }

  Facet * GetFacet() { return facet_; }

protected:
  Epetra_SerialDenseMatrix xyz_;
  Facet * facet_;
};

class Tri3BoundaryCell : public BoundaryCell
{
public:

  Tri3BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet )
    : BoundaryCell( xyz, facet )
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::tri3; }

};

class Quad4BoundaryCell : public BoundaryCell
{
public:

  Quad4BoundaryCell( const Epetra_SerialDenseMatrix & xyz, Facet * facet )
    : BoundaryCell( xyz, facet )
  {
  }

  virtual DRT::Element::DiscretizationType Shape() const { return DRT::Element::quad4; }

};

}
}

#endif
