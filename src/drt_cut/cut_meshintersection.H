#ifndef CUT_MESHINTERSECTION_H
#define CUT_MESHINTERSECTION_H

#include <Epetra_SerialDenseMatrix.h>

#include <Teuchos_RCP.hpp>

#include <vector>

#include "../drt_lib/drt_element.H"

#include "cut_meshhandle.H"
#include "cut_pointpool.H"
#include "cut_options.H"

namespace GEO
{
namespace CUT
{
  class Node;
  class Edge;
  class Side;
  class Element;
  class ElementHandle;

/// interface class for a surface mesh cut
class MeshIntersection
{
public:

  explicit MeshIntersection( int numcutmesh=1 )
    : pp_( Teuchos::rcp( new PointPool ) ),
      mesh_( options_, 1, pp_, false )
  {
    cut_mesh_.reserve( numcutmesh );
    for ( int i=0; i<numcutmesh; ++i )
    {
      cut_mesh_.push_back( Teuchos::rcp( new MeshHandle( options_, 1, pp_, true ) ) );
    }
  }

  void SetFindPositions( bool positions ) { options_.SetFindPositions( positions ); }

  ElementHandle * AddElement( int eid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype );

  SideHandle * AddCutSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype, int mi=0 );

  SideHandle * AddCutSide( int sid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype, int mi=0 );

  void Cut( bool include_inner );

  void CreateNodalDofSetNEW ( bool include_inner, DRT::Discretization & dis);

  void FindNodalCellSets( bool include_inner,
                          std::set<int> & eids,
                          std::map<int, ElementHandle*> & sourrounding_elements,
                          std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_inside,
                          std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_outside,
                          std::vector<plain_volumecell_set> & cell_sets_inside,
                          std::vector<plain_volumecell_set> & cell_sets_outside,
                          std::vector<plain_volumecell_set> & cell_sets );


  void ConnectNodalDOFSets( std::vector<Node *> & nodes,
                            bool include_inner,
                            const std::vector<plain_volumecell_set> & connected_vc_sets,
                            std::vector<std::vector<int> > &    nodaldofset_vc_sets);

  Node * GetNode( int nid ) const;

  ElementHandle * GetElement( int eid ) const;

  SideHandle * GetCutSide( int sid, int mi=0 ) const;

  void PrintCellStats();

  void Status();

  Mesh & NormalMesh() { return mesh_.LinearMesh(); }

  MeshHandle & GetMeshHandle() { return mesh_; }

  Mesh & CutMesh( int i=0 ) { return cut_mesh_[i]->LinearMesh(); }

  void DumpGmshNodalCellSet( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets, DRT::Discretization & dis );

  void DumpGmshCellSets( std::vector<plain_volumecell_set> & cell_sets,  DRT::Discretization & dis );

  void DumpGmshNumDOFSets( string filename, bool include_inner, DRT::Discretization & dis );

  void DumpGmshVolumeCells( std::string name, bool include_inner );

  void DumpGmshIntegrationCells( std::string name );

private:

  Teuchos::RCP<PointPool> pp_;

  Options options_;
  MeshHandle mesh_;
  std::vector<Teuchos::RCP<MeshHandle> > cut_mesh_;
};

}
}

#endif
