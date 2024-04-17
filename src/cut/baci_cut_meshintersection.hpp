/*----------------------------------------------------------------------*/
/*! \file

\brief provides the specific functionality for cutting a mesh with other meshes

\level 3
 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_MESHINTERSECTION_HPP
#define FOUR_C_CUT_MESHINTERSECTION_HPP

#include "baci_config.hpp"

#include "baci_cut_parentintersection.hpp"

FOUR_C_NAMESPACE_OPEN

namespace DRT
{
  class Discretization;
}

namespace CORE::GEO
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
    class MeshIntersection : public virtual ParentIntersection
    {
      typedef ParentIntersection my;


     public:
      /// constructur for MeshIntersecton class
      explicit MeshIntersection(int numcutmesh = 1, int myrank = -1) : ParentIntersection(myrank)
      {
        cut_mesh_.reserve(numcutmesh);
        for (int i = 0; i < numcutmesh; ++i)
        {
          cut_mesh_.push_back(Teuchos::rcp(new MeshHandle(options_, 1, pp_, true, myrank)));
        }
      }

      /*========================================================================*/
      //! @name Add functionality for elements and sides
      /*========================================================================*/

      /// add this background element if it falls within the bounding box of cut mesh
      ElementHandle* AddElement(int eid, const std::vector<int>& nids,
          const CORE::LINALG::SerialDenseMatrix& xyz, CORE::FE::CellType distype,
          const double* lsv = nullptr);

      /// add a side of the cut mesh and return the sidehandle (e.g. quadratic sidehandle for
      /// quadratic sides)
      SideHandle* AddCutSide(
          int sid, const std::vector<int>& nids, CORE::FE::CellType distype, int mi = 0);

      /// add a side of the cut mesh and return the sidehandle (e.g. quadratic sidehandle for
      /// quadratic sides)
      SideHandle* AddCutSide(int sid, const std::vector<int>& nids,
          const CORE::LINALG::SerialDenseMatrix& xyz, CORE::FE::CellType distype, int mi = 0);

      /// build the static search tree for the collision detection in the self cut
      void BuildSelfCutTree();

      /// build the static search tree for the collision detection
      void BuildStaticSearchTree();

      /*========================================================================*/
      //! @name Cut functionality, routines
      /*========================================================================*/

      /// standard cut routine for non-parallel frameworks and cuttest
      void CutTest_Cut(bool include_inner,
          INPAR::CUT::VCellGaussPts VCellgausstype = INPAR::CUT::VCellGaussPts_Tessellation,
          INPAR::CUT::BCellGaussPts BCellgausstype = INPAR::CUT::BCellGaussPts_Tessellation,
          bool tetcellsonly = true, bool screenoutput = true,
          bool do_Cut_Positions_Dofsets = false);  // for cuttest with closed cutsides this option
                                                   // can be activated, otherwise this will fail!

      /// handles cut sides which cut each other
      void Cut_SelfCut(bool include_inner, bool screenoutput = true) override;

      /// detects if a side of the cut mesh possibly collides with an element of the background mesh
      void Cut_CollisionDetection(bool include_inner, bool screenoutput = true) override;

      /// Routine for cutting the mesh. This creates lines, facets, volumecells and quadrature rules
      void Cut_MeshIntersection(bool screenoutput = true) override;

      /// Routine for deciding the inside-outside position. This creates the dofset data, just
      /// serial
      void Cut_Positions_Dofsets(bool include_inner, bool screenoutput = true);

      /*========================================================================*/
      //! @name get routines for nodes, elements, sides, mesh, meshhandles
      /*========================================================================*/

      /// get the cut mesh's side based on side id
      SideHandle* GetCutSide(int sid, int mi = 0) const;

      /// get the cut mesh based on mesh id
      Mesh& CutMesh(int i = 0) { return cut_mesh_[i]->LinearMesh(); }

     private:
      /*========================================================================*/
      //! @name private class variables
      /*========================================================================*/

      std::vector<Teuchos::RCP<MeshHandle>> cut_mesh_;  ///< a vector of cut_meshes
    };

  }  // namespace CUT
}  // namespace CORE::GEO

FOUR_C_NAMESPACE_CLOSE

#endif
