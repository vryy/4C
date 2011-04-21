#ifndef CUT_MESHINTERSECTION_H
#define CUT_MESHINTERSECTION_H

#include <Epetra_SerialDenseMatrix.h>

#include <Teuchos_RCP.hpp>

#include <vector>

#include "../drt_lib/drt_element.H"

#include "cut_meshhandle.H"
#include "cut_pointpool.H"

namespace GEO
{
namespace CUT
{
  class Node;
  class Edge;
  class Side;
  class Element;
  class ElementHandle;

class MeshIntersection
{
public:

  explicit MeshIntersection( int numcutmesh=1 )
    : pp_( Teuchos::rcp( new PointPool ) ),
      mesh_( 1, pp_, false )
  {
    cut_mesh_.reserve( numcutmesh );
    for ( int i=0; i<numcutmesh; ++i )
    {
      cut_mesh_.push_back( Teuchos::rcp( new MeshHandle( 1, pp_, true ) ) );
    }
  }

  ElementHandle * AddElement( int eid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype );

  SideHandle * AddCutSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype, int mi=0 );

  SideHandle * AddCutSide( int sid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype, int mi=0 );

  void Cut( bool include_inner );

  Node * GetNode( int nid ) const;

  ElementHandle * GetElement( int eid ) const;

  SideHandle * GetCutSide( int sid, int mi=0 ) const;

  void Status();

  Mesh & NormalMesh() { return mesh_.LinearMesh(); }

  Mesh & CutMesh( int i=0 ) { return cut_mesh_[i]->LinearMesh(); }

private:

  Teuchos::RCP<PointPool> pp_;

  MeshHandle mesh_;
  std::vector<Teuchos::RCP<MeshHandle> > cut_mesh_;
};

}
}

#endif
