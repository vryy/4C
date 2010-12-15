#ifndef CUT_LEVELSETINTERSECTION_H
#define CUT_LEVELSETINTERSECTION_H

#include <Epetra_SerialDenseMatrix.h>

#include "../drt_lib/drt_element.H"

#include "cut_mesh.H"
#include "cut_pointpool.H"
#include "cut_levelsetside.H"

namespace GEO
{
namespace CUT
{
  class Node;
  class Edge;
  class Side;
  class Element;

class LevelSetIntersection
{
public:

  LevelSetIntersection()
    : pp_( Teuchos::rcp( new PointPool ) ),
      mesh_( 1, pp_, false ),
      side_( 1 )
  {
  }

  void AddElement( int eid,
                   const std::vector<int> & nids,
                   const Epetra_SerialDenseMatrix & xyz,
                   const double * lsv,
                   DRT::Element::DiscretizationType distype );

  void Cut();

private:

  Teuchos::RCP<PointPool> pp_;

  Mesh mesh_;
  LevelSetSide side_;
};

}
}

#endif
