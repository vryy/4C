/*!-----------------------------------------------------------------------------------------------*
\file cut_levelsetintersection.H

\brief provides the basic functionality for cutting a mesh with a level set function

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef CUT_LEVELSETINTERSECTION_H
#define CUT_LEVELSETINTERSECTION_H

#include "cut_parentintersection.H"


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
class LevelSetIntersection : public virtual ParentIntersection
{
  typedef ParentIntersection my;


public:

  /// constructur for LevelSetIntersecton class
  LevelSetIntersection(int myrank =-1, bool create_side = true);

  ///
  void AddCutSide( int levelset_sid );

  ///
  bool HasLSCuttingSide(int sid){ return true; /*return sid == side_->Id();*/ };

  /*========================================================================*/
  //! @name Cut functionality, routines
  /*========================================================================*/

  /*!
  \brief Performs the cut of the mesh with the level set
  */
  void Cut_Mesh( bool screenoutput = false );

  /*!
  \brief Performs all the level set cut operations including find positions and triangulation. (Used for the test cases)
   */
  void Cut( bool include_inner=true, bool screenoutput = false, INPAR::CUT::VCellGaussPts VCellGP=INPAR::CUT::VCellGaussPts_Tessellation);



  /*========================================================================*/
  //! @name Add functionality for elements
  /*========================================================================*/

  /// add this background element if it is cut (determined by level set)
  GEO::CUT::ElementHandle * AddElement(
      int eid,
      const std::vector<int> & nids,
      const Epetra_SerialDenseMatrix & xyz,
      DRT::Element::DiscretizationType distype,
      const double * lsv,
      const bool lsv_only_plus_domain = false
  );

protected:

  /*========================================================================*/
  //! @name private class variables
  /*========================================================================*/

  Teuchos::RCP<LevelSetSide> side_;
};

}
}

#endif
