#ifndef CUT_MESHINTERSECTION_H
#define CUT_MESHINTERSECTION_H

#include <Epetra_SerialDenseMatrix.h>

#include "../drt_lib/drt_element.H"

#include "cut_mesh.H"
#include "cut_pointpool.H"

namespace GEO
{
namespace CUT
{
  class Node;
  class Edge;
  class Side;
  class Element;

class MeshIntersection
{
public:

  MeshIntersection()
    : pp_( Teuchos::rcp( new PointPool ) ),
      mesh_( pp_, false ),
      cut_mesh_( pp_, true )
  {}

  void AddElement( int eid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype );

  void AddCutSide( int sid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype );

  void Cut( CellGenerator * generator=NULL );

  GEO::CUT::Side * GetCutSides( int sid );

private:

  Teuchos::RCP<PointPool> pp_;

  Mesh mesh_;
  Mesh cut_mesh_;
};

}
}

#endif
