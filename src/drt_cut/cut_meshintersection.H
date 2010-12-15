#ifndef CUT_MESHINTERSECTION_H
#define CUT_MESHINTERSECTION_H

#include <Epetra_SerialDenseMatrix.h>

#include <Teuchos_RCP.hpp>

#include <vector>

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

  explicit MeshIntersection( int numcutmesh=1 )
    : pp_( Teuchos::rcp( new PointPool ) ),
      mesh_( 1, pp_, false )
  {
    cut_mesh_.reserve( numcutmesh );
    for ( int i=0; i<numcutmesh; ++i )
    {
      cut_mesh_.push_back( Teuchos::rcp( new Mesh( 1, pp_, true ) ) );
    }
  }

  void AddElement( int eid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype );

  void AddCutSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype, int mi=0 );

  void AddCutSide( int sid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype, int mi=0 );

  void Cut( CellGenerator * generator=NULL );

  void SelfCut();

  Element * GetElement( int eid );

  Side * GetCutSide( int sid, int mi=0 );

  void Status();

  Mesh & NormalMesh() { return mesh_; }

  Mesh & CutMesh( int i=0 ) { return *cut_mesh_[i]; }

private:

  Teuchos::RCP<PointPool> pp_;

  Mesh mesh_;
  std::vector<Teuchos::RCP<Mesh> > cut_mesh_;
};

}
}

#endif
