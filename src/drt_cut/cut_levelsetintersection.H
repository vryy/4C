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

#include <Epetra_SerialDenseMatrix.h>

#include "cut_parentintersection.H"

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
class LevelSetIntersection : public ParentIntersection
{
public:

  /// constructur for LevelSetIntersecton class
  LevelSetIntersection(int myrank =-1);

  /*========================================================================*/
  //! @name Cut functionality, routines
  /*========================================================================*/
  /*!
  \brief Performs all the level set cut operations.(Used for the test cases)
   */
  void Cut( bool include_inner=true, bool screenoutput = false );

  /*!
  \brief Performs the cut of the mesh with the level set
  */
  void Cut_Mesh( bool include_inner, bool screenoutput = false );

  /*========================================================================*/
  //! @name Add functionality for elements
  /*========================================================================*/

  /// add this background element if it is cut (determined by level set)
  void AddElement( int eid,
                   const std::vector<int> & nids,
                   const Epetra_SerialDenseMatrix & xyz,
                   const double * lsv,
                   DRT::Element::DiscretizationType distype );

  /*========================================================================*/
  //! @name statistics, output
  /*========================================================================*/

  /// status
  void Status(INPAR::CUT::VCellGaussPts gausstype = INPAR::CUT::VCellGaussPts_Tessellation);

private:

  /*========================================================================*/
  //! @name private class variables
  /*========================================================================*/

  Teuchos::RCP<LevelSetSide> side_;
};

}
}

#endif
