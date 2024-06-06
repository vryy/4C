/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief provides the basic functionality for cutting a mesh with a level set function and/or a
       mesh


\level 2
*/
/*------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_COMBINTERSECTION_HPP
#define FOUR_C_CUT_COMBINTERSECTION_HPP

#include "4C_config.hpp"

#include "4C_cut_levelsetintersection.hpp"
#include "4C_cut_meshintersection.hpp"
#include "4C_cut_parentintersection.hpp"

FOUR_C_NAMESPACE_OPEN


namespace LinAlg
{
  class SerialDenseMatrix;
}

namespace Core::Geo
{
  namespace Cut
  {
    class Node;
    class Edge;
    class Side;
    class Element;

    /*!
    \brief Interface class for the level set cut.
    */
    class CombIntersection : public LevelSetIntersection, public MeshIntersection
    {
     public:
      /// constructur for LevelSetIntersecton class
      CombIntersection(int myrank);

      void Cut(bool screenoutput);

      void FindNodePositions();

      void add_element(int eid, const std::vector<int>& nids,
          const Core::LinAlg::SerialDenseMatrix& xyz, Core::FE::CellType distype, const double* lsv,
          const bool lsv_only_plus_domain);

      void AddLevelSetSide(int levelset_sid);

      void add_mesh_cutting_side(int sid, const std::vector<int>& nids,
          const Core::LinAlg::SerialDenseMatrix& xyz, Core::FE::CellType distype, int mi);
    };

  }  // namespace Cut
}  // namespace Core::Geo
FOUR_C_NAMESPACE_CLOSE

#endif
