#ifndef GEO_INTERSECTION_H
#define GEO_INTERSECTION_H

#ifdef QHULL

#include <map>
#include <vector>

#include <Teuchos_RCP.hpp>
#include <Epetra_CrsGraph.h>

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_dofset.H"
#include "geo_utils.H"


namespace DRT
{
  class Discretization;
}

namespace GEO
{
  namespace CUT
  {
    class MeshIntersection;
    class ElementHandle;
    class Node;
  }

//   class CutDofSet : public DRT::DofSet
//   {
//   public:

//     CutDofSet( GEO::CUT::MeshIntersection * mesh, bool include_inner )
//       : mesh_( mesh ),
//         include_inner_( include_inner )
//     {
//     }

//     int PlainNumDofPerNode( const DRT::Node & node, unsigned dspos ) const;

//     void Dof( DRT::Node & node, int nodaldofset, std::vector<int> & dofs ) const;

//   protected:

//     /// get number of nodal dofs for this element at this node
//     virtual int NumDofPerNode( const DRT::Node & node, unsigned dspos ) const;

//   private:

//     GEO::CUT::MeshIntersection * mesh_;
//     bool include_inner_;
//   };

  /// contains the cut
  class CutWizard
  {
  public:

    /*!
    \brief Constructor
    */
    CutWizard( DRT::Discretization & dis, bool include_inner, int numcutmesh=1 );

    /*!
    \brief Set the value of positions in the cut_option.H
    */
    void SetFindPositions( bool positions );

    /*!
    \brief Add this cut element to the cut libraries
    */
    void AddCutSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze );

    /*!
    \brief Add this background mesh element to the cut libraries
    */
    void AddElement( DRT::Element * ele );

    /*!
    \brief Get this element from the cut libraries
     */
    GEO::CUT::ElementHandle * GetElement( DRT::Element * ele );

    /*!
    \brief Get this node from the cut libraries
     */
    GEO::CUT::Node * GetNode( int nid );

    /*!
    \brief Cut routine for parallel framework in XFSI and XFLUIDFLUID
    */
    void CutParallel( bool include_inner, std::string VCellgausstype, std::string BCellgausstype );

    /*!
    \brief Routine for finding node positions and computing volumecell dofsets in a parallel way
    */
    void CutParallel_FindPositionDofSets(bool include_inner, bool parallel);

    /*!
    \brief Cut routine for standard non-parallel framework which is used only for cut_test
    */
    void Cut( bool include_inner, std::string VCellgausstype, std::string BCellgausstype );

    /*!
    \brief Cut routine for standard combustion and two phase flow framework
    */
    void Cut( std::map< int, DomainIntCells >& domainintcells,
              std::map< int, BoundaryIntCells >& boundaryintcells,
              std::string VCellgausstype,
              std::string BCellgausstype);

//     Teuchos::RCP<CutDofSet> CreateDofSet();

//     Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( const CutDofSet & dofset, const Epetra_Map & dbcmap );

    Teuchos::RCP<Epetra_Map> CutElementMap() const { return cutelementmap_; }

    GEO::CUT::MeshIntersection & Mesh() { return *mesh_; }

    /*!
    \brief Print the number of volumecells and boundarycells generated over the whole mesh during the cut
     */
    void PrintCellStats();

    /*!
    \brief Write the DOF details of the nodes
     */
    void DumpGmshNumDOFSets( bool include_inner );

    /*!
    \brief Write volumecell output in GMSH format throughout the domain
     */
    void DumpGmshVolumeCells( bool include_inner );

    /*!
    \brief Write the integrationcells and boundarycells in GMSH format throughout the domain
     */
    void DumpGmshIntegrationCells();

  private:

    DRT::Discretization & dis_;
    Teuchos::RCP<GEO::CUT::MeshIntersection> mesh_;
    Teuchos::RCP<Epetra_Map> cutelementmap_;
    bool include_inner_;
  };

  /*!
   * \brief  computes the intersection of two discretizations.
   *
   * This methods localizes the interface of the
   * cutter discretization within the xfem discretization. Each intersected xfem element
   * is tetrahedralized with help of a Constrained Delaunay Tetrahedralization.
   * The output consists of a list of qadratic tetrahedral volume integration cells and
   * a list of quadratic triangular boundary integration cells
   * representing the curved interface.
   */
  void computeIntersection(
    const Teuchos::RCP<DRT::Discretization>              xfemdis,                 ///< xfem discretization
    const Teuchos::RCP<DRT::Discretization>              cutterdis,               ///< cutter discretization
    const std::map<int,LINALG::Matrix<3,1> >&   currentcutterpositions,  ///< current positions of all cutter nodes
    const std::map<int,LINALG::Matrix<3,2> >&   currentXAABBs,           ///< XAABBs of current elements
    std::map< int, DomainIntCells >&                 domainintcells,          ///< domainintegrationcells for each intersected element
    std::map< int, BoundaryIntCells >&               boundaryintcells,        ///< boundaryintegrationcells for each intersected element
    const std::map<int,int>&                       labelPerElementId,
    const std::vector<int>&                        MovingFluideleGIDs,
    std::string VCellgausstype,                                          //volumecell gauss points generation
    std::string BCellgausstype                                          //boundarycell gauss points generation
    );

}

#endif
#endif
