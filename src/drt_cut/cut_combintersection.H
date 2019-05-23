/*-----------------------------------------------------------------------------------------------*/
/*!
\file cut_combintersection.H

\brief provides the basic functionality for cutting a mesh with a level set function and/or a
       mesh

\maintainer  Christoph Ager
             ager@lnm.mw.tum.de
             http://www.lnm.mw.tum.de
             089 - 289-15249

\level 2
*/
/*------------------------------------------------------------------------------------------------*/

#ifndef CUT_COMBINTERSECTION_H
#define CUT_COMBINTERSECTION_H

#include "cut_parentintersection.H"
#include "cut_meshintersection.H"
#include "cut_levelsetintersection.H"


class Epetra_SerialDenseMatrix;

namespace GEO
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

      void AddElement(int eid, const std::vector<int>& nids, const Epetra_SerialDenseMatrix& xyz,
          DRT::Element::DiscretizationType distype, const double* lsv,
          const bool lsv_only_plus_domain);

      void AddLevelSetSide(int levelset_sid);

      void AddMeshCuttingSide(int sid, const std::vector<int>& nids,
          const Epetra_SerialDenseMatrix& xyz, DRT::Element::DiscretizationType distype, int mi);
    };

  }  // namespace CUT
}  // namespace GEO
#endif
