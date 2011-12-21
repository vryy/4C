#ifndef GEO_INTERSECTION_H
#define GEO_INTERSECTION_H

#ifdef QHULL

#include <map>
#include <vector>

#include <Teuchos_RCP.hpp>

#include <Epetra_CrsGraph.h>

#include "../linalg/linalg_fixedsizematrix.H"
#include "../drt_lib/drt_dofset.H"

#include "integrationcell.H"

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

    CutWizard( DRT::Discretization & dis, bool include_inner, int numcutmesh=1 );

    void SetFindPositions( bool positions );

    void AddCutSide( int mi, DRT::Element * ele, const Epetra_SerialDenseMatrix & xyze );

    void AddElement( DRT::Element * ele );

    GEO::CUT::ElementHandle * GetElement( DRT::Element * ele );

    GEO::CUT::Node * GetNode( int nid );

    void Cut( bool include_inner, std::string VCellgausstype, std::string BCellgausstype );

    void CreateNodalDofSet ( bool include_inner, DRT::Discretization & dis );

    void Cut( std::map< int, DomainIntCells >& domainintcells,
              std::map< int, BoundaryIntCells >& boundaryintcells,
              std::string VCellgausstype,
              std::string BCellgausstype);

//     Teuchos::RCP<CutDofSet> CreateDofSet();

//     Teuchos::RCP<Epetra_CrsGraph> MatrixGraph( const CutDofSet & dofset, const Epetra_Map & dbcmap );

    Teuchos::RCP<Epetra_Map> CutElementMap() const { return cutelementmap_; }

    GEO::CUT::MeshIntersection & Mesh() { return *mesh_; }

    void PrintCellStats();

    void DumpGmshNumDOFSets( bool include_inner );

    void DumpGmshVolumeCells( bool include_inner );

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
