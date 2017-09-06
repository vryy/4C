/*-----------------------------------------------------------------------------------------------*/
/**
\file cut_parentintersection.H

\brief provides the basic functionality for cutting a mesh

<pre>
\maintainer  Ager Christoph
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249
</pre>

\level 2
*/
/*------------------------------------------------------------------------------------------------*/

#ifndef CUT_PARENTINTERSECTION_H
#define CUT_PARENTINTERSECTION_H

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
\brief Interface class for the general mesh cut. This class contains shared functionality between
       level set intersection and mesh intersection routines.
*/
class ParentIntersection
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
          int                                set_index,            ///< set index for Volumecell
          bool                               inside_cell,          ///< cell inside or outside
          std::vector<LINALG::Matrix<3,1> >& cut_points_coords,    ///< coordinates of cut_points
          int                                peid,                 ///< parent element id
          std::map<int,int> &                node_dofsetnumber_map ///< for the current volumecell in a parent element, for each node (nid) the current dofset number
      ) :
        set_index_(set_index),
        inside_cell_(inside_cell),
        peid_(peid)
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
        std::cout << "Volumecell-DofSetData: " << std::endl;

        // print parent element Id
        std::cout << "\tparent element id: " << peid_ << std::endl;

        // print node_dofsetnumber_map
        std::cout << "\tnode_dofsetnumber_map" << std::endl;
        for(std::map<int,int>::iterator i=node_dofsetnumber_map_.begin(); i!= node_dofsetnumber_map_.end(); ++i)
        {
          std::cout << "\t\tnodeId \t" << i->first << "\t dofsetnumber \t" << i->second << std::endl;
        }
      }

      //! destructor
      ~DofSetData()
      {
        return;
      }

      int                                set_index_;              ///< set index for Volumecell
      bool                               inside_cell_;            ///< bool inside or outside cell
      std::vector<LINALG::Matrix<3,1> >  cut_points_coords_;      ///< coordinates for points of Volumecell
      int                                peid_;                   ///< parent element Id for volumecell
      std::map<int,int>                  node_dofsetnumber_map_;  ///< node Ids and dofset numbers w.r.t volumecell for that data have to be communicated


    private:

    }; // end class DofSetData


  /// constructur for ParentIntersecton class
  ParentIntersection( int myrank = -1)
    : pp_( Teuchos::rcp( new PointPool ) ),
      mesh_( options_, 1, pp_, false,myrank ),
      myrank_(myrank)
  {
  }

  /// destructor
  virtual ~ParentIntersection(){};

  /*========================================================================*/
  //! @name set and get routines for options
  /*========================================================================*/

  /// set the option if positions have to be determined or not
  void SetFindPositions( bool positions ) { options_.SetFindPositions( positions ); }

  /// set the option if positions have to be determined or not
  void SetNodalDofSetStrategy( INPAR::CUT::NodalDofSetStrategy nodal_dofset_strategy) { options_.SetNodalDofSetStrategy( nodal_dofset_strategy ); }

  /// Set the position for the boundary cell creation
  void SetGenBoundaryCellPosition( INPAR::CUT::BoundaryCellPosition gen_bcell_position )
  { options_.SetGenBoundaryCellPosition( gen_bcell_position ); }

  /// get the options
  void GetOptions(Options & options){options = options_;};

  /*========================================================================*/
  //! @name Cut functionality routines
  /*========================================================================*/


  virtual void Cut_SelfCut(bool include_inner, bool screenoutput)
  {
    if(myrank_==0 and screenoutput) IO::cout << "\t * 2/6 Cut_SelfCut ... not performed";
    return;
  };

  virtual void Cut_CollisionDetection(bool include_inner, bool screenoutput){
    if(myrank_==0 and screenoutput) IO::cout << "\t * 3/6 Cut_CollisionDetection ... not performed";
    return;
  };

  virtual void Cut_Mesh( bool screenoutput )
  {
    if(myrank_==0 and screenoutput) IO::cout << "\t * 4/6 Cut_Mesh (LevelSet-Cut) ... not performed";
    return;
  };

  virtual void Cut_MeshIntersection( bool screenoutput )
  {
    if(myrank_==0 and screenoutput) IO::cout << "\t * 4/6 Cut_MeshIntersection (Mesh-Cut) ... not performed";
    return;
  };

  /*!
  \brief The routine which splits the volumecell into integrationcells by tessellation, or create Gaussian
  integration rules by moment fitting equations
   */
  void Cut_Finalize( bool include_inner,
                     INPAR::CUT::VCellGaussPts VCellgausstype,
                     INPAR::CUT::BCellGaussPts BCellgausstype,
                     bool tetcellsonly,
                     bool screenoutput);

  /*========================================================================*/
  //! @name nodal dofset routines
  /*========================================================================*/

  /// Create nodal dofset sets within the parallel cut framework
  void CreateNodalDofSet ( bool include_inner, const DRT::Discretization & dis);

  /// fill parallel DofSetData with information that has to be communicated
  void FillParallelDofSetData(
      std::vector<Teuchos::RCP<DofSetData> > & parallel_dofSetData,
      const DRT::Discretization& dis,
      bool include_inner
      );

  /// create parallel DofSetData for a volumecell that has to be communicated
  void CreateParallelDofSetDataVC(
      std::vector<Teuchos::RCP<DofSetData> >& parallel_dofSetData,
      int eid,
      int set_index,
      bool inside,
      VolumeCell * cell,
      std::map<int,int>& node_dofset_map
  );

  /// find cell sets around each node (especially for quadratic elements)
  void FindNodalCellSets( bool include_inner,
                          std::set<int> & eids,
                          std::vector<int> & sourrounding_elements,
                          std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_inside,
                          std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets_outside,
                          std::vector<plain_volumecell_set> & cell_sets_inside,
                          std::vector<plain_volumecell_set> & cell_sets_outside,
                          std::vector<plain_volumecell_set> & cell_sets );

  /// connect sets of volumecells for neighboring elements around a node
  void ConnectNodalDOFSets( std::vector<Node *> &                      nodes,
                            bool                                       include_inner,
                            const DRT::Discretization&                 dis,
                            const std::vector<plain_volumecell_set> &  connected_vc_sets,
                            std::vector<std::vector<int> > &           nodaldofset_vc_sets,
                            std::vector<std::map<int,int> >&           vcsets_nid_dofsetnumber_map_toComm);


  /*========================================================================*/
  //! @name get routines for nodes, elements, sides, mesh, meshhandles
  /*========================================================================*/

  /// get the node based on node id
  Node * GetNode( int nid ) const;

  /// get the mesh's side based on node ids and return the side
  SideHandle * GetSide( std::vector<int>& nodeids ) const;

  /// get the mesh's side based on side id
  SideHandle * GetSide( int sid ) const;

  /// get the mesh's element based on element id
  ElementHandle * GetElement( int eid ) const;

  /// get the linear mesh
  Mesh & NormalMesh() { return mesh_.LinearMesh(); }

  /// get the mesh handle
  MeshHandle & GetMeshHandle() { return mesh_; }

  /*========================================================================*/
  //! @name GMSH output
  /*========================================================================*/

  /// write gmsh debug output for nodal cell sets
  void DumpGmshNodalCellSet( std::map<Node*, std::vector<plain_volumecell_set> > & nodal_cell_sets, const DRT::Discretization & dis );

  /// write gmsh debug output for CellSets
  void DumpGmshCellSets( std::vector<plain_volumecell_set> & cell_sets, const DRT::Discretization & dis );

  /// write gmsh cut output for number of dofsets and the connected vc sets
  void DumpGmshNumDOFSets( std::string filename, bool include_inner, const DRT::Discretization & dis );

  /// write gmsh output for volumecells
  void DumpGmshVolumeCells( std::string name, bool include_inner );

  /// write gmsh output for volumecells
  void DumpGmshIntegrationCells( std::string name );

  /// write gmsh output for volumecells
  void DumpGmshVolumeCells( std::string name );


  /*========================================================================*/
  //! @name statistics, output
  /*========================================================================*/

  /// print cell statistics
  void PrintCellStats();

  /// status
  virtual void Status(INPAR::CUT::VCellGaussPts gausstype = INPAR::CUT::VCellGaussPts_Tessellation);

protected:

  /*!
  \brief Compares Tesselation and DirectDivergence Cut.
   */
  void DebugCut(Mesh& m);

  /*========================================================================*/
  //! @name private class variables
  /*========================================================================*/

  Teuchos::RCP<PointPool> pp_;                       ///< pointpool (octTree) whose nodes consist of bounding boxes, each bb contains a set of Teuchos::RCPs to points
  Options options_;                                  ///< options
  MeshHandle mesh_;                                  ///< the background mesh
  int myrank_;                                       ///< my processor Id

};

}
}

#endif
