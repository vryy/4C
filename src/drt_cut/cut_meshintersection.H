/*!-----------------------------------------------------------------------------------------------*
\file cut_meshintersection.H

\brief provides the basic functionality for cutting a mesh

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_MESHINTERSECTION_H
#define CUT_MESHINTERSECTION_H

#include <vector>

#include "cut_meshhandle.H"
#include "cut_pointpool.H"
#include "cut_options.H"

namespace DRT
{
  class Discretization;
}

namespace GEO
{
namespace CUT
{
  class Node;
  class Edge;
  class Side;
  class Element;
  class ElementHandle;

/*!
\brief Interface class for the surface mesh cut. The surface mesh is in general triangulated.
*/
class MeshIntersection
{
public:


    /*!
    \brief This class holds data for all volumecells for that dofsets (dofset numbers) of non-row nodes have
       to be communicated between processors.
     */
    class DofSetData
    {
    public:

      //! constructor for creating volumecells dofset data for communication between processors (constructor during the Robin round)
      DofSetData(
          int                               set_index,         // set index for Volumecell
          bool                              inside_cell,       // cell inside or outside
          std::vector<LINALG::Matrix<3,1> > cut_points_coords, // coordinates of cut_points
          int peid,
          std::map<int,int> & node_dofsetnumber_map
      ) :
        set_index_(set_index),
        inside_cell_(inside_cell),
        peid_(peid)
      //node_dofsetnumber_map_(node_dofsetnumber_map)
      {
        std::copy(node_dofsetnumber_map.begin(), node_dofsetnumber_map.end(), std::inserter(node_dofsetnumber_map_, node_dofsetnumber_map_.begin()) );
        std::copy(cut_points_coords.begin(), cut_points_coords.end(), std::inserter(cut_points_coords_, cut_points_coords_.begin()) );
      }

      /*!
      \brief print the dofset data to screen
       */
      void print()
      {
        // print volumecell information
        std::cout << "Volumecell-DofSetData: " << endl;

        // print parent element Id
        std::cout << "\tparent element id: " << peid_ << endl;

        // print node_dofsetnumber_map
        std::cout << "\tnode_dofsetnumber_map" << endl;
        for(std::map<int,int>::iterator i=node_dofsetnumber_map_.begin(); i!= node_dofsetnumber_map_.end(); ++i)
        {
          std::cout << "\t\tnodeId \t" << i->first << "\t dofsetnumber \t" << i->second << endl;
        }
      }

      //! destructor
      ~DofSetData()
      {
        return;
      }

      int                                set_index_;              // set index for Volumecell
      bool                               inside_cell_;            // bool inside or outside cell
      std::vector<LINALG::Matrix<3,1> >  cut_points_coords_;      // coordinates for points of Volumecell
      int                                peid_;                   // parent element Id for volumecell
      std::map<int,int>                  node_dofsetnumber_map_;  // node Ids and dofset numbers that for that data has to be communicated


    private:

    }; // end class DofSetData



  explicit MeshIntersection( int numcutmesh=1, int myrank = -1)
    : pp_( Teuchos::rcp( new PointPool ) ),
      mesh_( options_, 1, pp_, false ),
      myrank_(myrank)
  {
    cut_mesh_.reserve( numcutmesh );
    for ( int i=0; i<numcutmesh; ++i )
    {
      cut_mesh_.push_back( Teuchos::rcp( new MeshHandle( options_, 1, pp_, true ) ) );
    }
  }


  void SetFindPositions( bool positions ) { options_.SetFindPositions( positions ); }

  void GetOptions(Options & options){options = options_;};

  /*!
  \brief Add this background element if it falls within the bounding box of cut mesh
  */
  ElementHandle * AddElement( int eid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype );

  SideHandle * AddCutSide( int sid, const std::vector<int> & nids, DRT::Element::DiscretizationType distype, int mi=0 );

  SideHandle * AddCutSide( int sid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype, int mi=0 );

  /// standard cut routine for non-parallel frameworks and cuttest
  void Cut( bool include_inner, std::string VCellgausstype="Tessellation", std::string BCellgausstype="Tessellation" );

  /*!
  \brief Routine for cutting the mesh. This creates lines, facets, volumecells and quadrature rules
   */
  void Cut_Mesh( bool include_inner );

  /*!
  \brief Routine for deciding the inside-outside position. This creates the dofset data, just serial
   */
  void Cut_Positions_Dofsets( bool include_inner );

  /*!
  \brief The routine which splits the volumecell into integrationcells by tessellation, or create Gaussian
  integration rules by moment fitting equations
   */
  void Cut_Finalize( bool include_inner, std::string VCellgausstype, std::string BCellgausstype);

  /*!
  \brief Create nodal dofset sets within the parallel cut framework
   */
  void CreateNodalDofSetNEW ( bool include_inner, DRT::Discretization & dis);

  /// fill parallel DofSetData with information that has to be communicated
  void FillParallelDofSetData(RCP<std::vector<DofSetData> > parallel_dofSetData, DRT::Discretization& dis);

  /// create parallel DofSetData for a volumecell that has to be communicated
  void CreateParallelDofSetDataVC( RCP<std::vector<DofSetData> > parallel_dofSetData , int eid, int set_index, bool inside, VolumeCell * cell, std::map<int,int>& node_dofset_map);

  /*!
  \brief find cell sets around each node (especially for quadratic elements)
  */
  void FindNodalCellSets( bool include_inner,
                          std::set<int> & eids,
                          std::map<int, ElementHandle*> & sourrounding_elements,
                          std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_inside,
                          std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_outside,
                          std::vector<plain_volumecell_set> & cell_sets_inside,
                          std::vector<plain_volumecell_set> & cell_sets_outside,
                          std::vector<plain_volumecell_set> & cell_sets );

  /// connect sets of volumecells for neighboring elements around a node
  void ConnectNodalDOFSets( std::vector<Node *> &                      nodes,
                            bool                                       include_inner,
                            DRT::Discretization&                       dis,
                            const std::vector<plain_volumecell_set> &  connected_vc_sets,
                            std::vector<std::vector<int> > &           nodaldofset_vc_sets,
                            std::vector<std::map<int,int> >&           vcsets_nid_dofsetnumber_map_toComm);

  Node * GetNode( int nid ) const;

  Side * GetSide( std::vector<int>& nodeids ) const;

  SideHandle * GetSide( int sid ) const;

  ElementHandle * GetElement( int eid ) const;

  SideHandle * GetCutSide( int sid, int mi=0 ) const;

  void PrintCellStats();

  void Status(std::string gausstype="Tessellation");

  Mesh & NormalMesh() { return mesh_.LinearMesh(); }

  MeshHandle & GetMeshHandle() { return mesh_; }

  Mesh & CutMesh( int i=0 ) { return cut_mesh_[i]->LinearMesh(); }

  void DumpGmshNodalCellSet( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets, DRT::Discretization & dis );

  void DumpGmshCellSets( std::vector<plain_volumecell_set> & cell_sets,  DRT::Discretization & dis );

  void DumpGmshNumDOFSets( string filename, bool include_inner, DRT::Discretization & dis );

  void DumpGmshVolumeCells( std::string name, bool include_inner );

  void DumpGmshIntegrationCells( std::string name );

  void DumpGmshVolumeCells( std::string name );






private:



  Teuchos::RCP<PointPool> pp_;

  Options options_;
  MeshHandle mesh_;
  std::vector<Teuchos::RCP<MeshHandle> > cut_mesh_;

  int myrank_;

};

}
}

#endif
