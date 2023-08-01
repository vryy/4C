/*-----------------------------------------------------------------------------------------------*/
/*! \file

\brief provides the basic functionality for cutting a mesh with a level set function and/or a
       mesh


\level 2
*/
/*------------------------------------------------------------------------------------------------*/

#ifndef BACI_CUT_COMBINTERSECTION_H
#define BACI_CUT_COMBINTERSECTION_H

#include "baci_cut_parentintersection.H"
#include "baci_cut_meshintersection.H"
#include "baci_cut_levelsetintersection.H"


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
          const CORE::LINALG::SerialDenseMatrix& xyz, ::DRT::Element::DiscretizationType distype,
          const double* lsv, const bool lsv_only_plus_domain);

      void AddLevelSetSide(int levelset_sid);

      void AddMeshCuttingSide(int sid, const std::vector<int>& nids,
          const CORE::LINALG::SerialDenseMatrix& xyz, ::DRT::Element::DiscretizationType distype,
          int mi);
    };

  }  // namespace CUT
}  // namespace CORE::GEO
#endif
