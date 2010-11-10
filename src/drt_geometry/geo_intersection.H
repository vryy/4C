#ifndef GEO_INTERSECTION_H
#define GEO_INTERSECTION_H

#ifdef QHULL

#include <map>
#include <vector>

#include <Teuchos_RCP.hpp>

#include "../linalg/linalg_fixedsizematrix.H"

#include "../drt_cut/cut_mesh.H"


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
  }

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
    const std::vector<int>&                        MovingFluideleGIDs
    );


  class CellGenerator : public GEO::CUT::CellGenerator
  {
  public:

    CellGenerator( DRT::Discretization & dis,
                   GEO::CUT::MeshIntersection & intersection,
                   std::map< int, DomainIntCells >   & domainintcells,
                   std::map< int, BoundaryIntCells > & boundaryintcells )
      : dis_( dis ),
        intersection_( intersection ),
        domainintcells_( domainintcells ),
        boundaryintcells_( boundaryintcells )
    {}

    virtual void Generate( GEO::CUT::Element* element, const tetgenio & out );

  private:

    /// xfem discretization
    DRT::Discretization & dis_;

    GEO::CUT::MeshIntersection & intersection_;

    /// domainintegrationcells for each intersected element
    std::map< int, DomainIntCells >   & domainintcells_;

    /// boundaryintegrationcells for each intersected element
    std::map< int, BoundaryIntCells > & boundaryintcells_;
  };
}

#endif
#endif
