/*!-----------------------------------------------------------------------------------------------*
\file cut_combintersection.H

\brief provides the basic functionality for cutting a mesh with a level set function

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

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
  class LevelSetSide;

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

  void AddElement(
      int eid,
      const std::vector<int> & nids,
      const Epetra_SerialDenseMatrix & xyz,
      DRT::Element::DiscretizationType distype,
      const double * lsv,
      const bool lsv_only_plus_domain
  );

  void AddLevelSetSide(int levelset_sid);

  void AddMeshCuttingSide( int sid, const std::vector<int> & nids, const Epetra_SerialDenseMatrix & xyz, DRT::Element::DiscretizationType distype, int mi );

};

}
}
#endif

