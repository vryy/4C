/*----------------------------------------------------------------------*/
/*!
\file crackUtils.H

\brief Utility functions for crack propagation problem

<pre>
Maintainer: Sudhakar
            sudhakar@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/
/*----------------------------------------------------------------------*/

#ifndef CRACKUTILS_H_
#define CRACKUTILS_H_

#include "../drt_lib/drt_discret.H"

namespace DRT {
namespace CRACK {
namespace UTILS{


  /*!
  \brief Update the given EpetraVector for the new map of discretization
   */
  void UpdateThisEpetraVectorCrack( const Teuchos::RCP<DRT::Discretization>& discret, Teuchos::RCP<Epetra_Vector>& vec,
                                  const std::map<int,int>& oldnewIds );

  /*!
  \brief Return the angle of the vector that points from a to b
   */
  double FindAngle( const double a[3], const double b[3] );

  /*!
  \brief convert the angle calculated in the range [-pi,pi] to [0,2pi]
   */
  void convertAngleTo_02PI_range( double &ang );

  /*!
  \brief Returns true when the given element has given nodeid
   */
  bool ElementHasThisNodeId( const DRT::Element* ele, int nodeid );

} /* namespace CRACK */
} /* namespace CRACK */
} /* namespace DRT */
#endif
