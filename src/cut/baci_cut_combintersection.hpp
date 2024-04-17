/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief provides the basic functionality for cutting a mesh with a level set function and/or a
       mesh


\level 2
*/
/*------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_CUT_COMBINTERSECTION_HPP
#define FOUR_C_CUT_COMBINTERSECTION_HPP

#include "baci_config.hpp"

#include "baci_cut_levelsetintersection.hpp"
#include "baci_cut_meshintersection.hpp"
#include "baci_cut_parentintersection.hpp"

FOUR_C_NAMESPACE_OPEN


namespace LINALG
{
  class SerialDenseMatrix;
}

namespace CORE::GEO
{
  namespace CUT
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

      void AddElement(int eid, const std::vector<int>& nids,
          const CORE::LINALG::SerialDenseMatrix& xyz, CORE::FE::CellType distype, const double* lsv,
          const bool lsv_only_plus_domain);

      void AddLevelSetSide(int levelset_sid);

      void AddMeshCuttingSide(int sid, const std::vector<int>& nids,
          const CORE::LINALG::SerialDenseMatrix& xyz, CORE::FE::CellType distype, int mi);
    };

  }  // namespace CUT
}  // namespace CORE::GEO
FOUR_C_NAMESPACE_CLOSE

#endif
