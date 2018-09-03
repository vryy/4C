/*!----------------------------------------------------------------------
\file contact_analytical.H

\brief A set of analytical solutions for convergence analysis of contact/meshtying methods

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*-----------------------------------------------------------------------*/
#ifndef CONTACT_ANALYTICAL_H
#define CONTACT_ANALYTICAL_H

#include "../linalg/linalg_fixedsizematrix.H"

/*!
\brief CONTACT: namespace of the BACI contact and meshtying algorithms
*/
namespace CONTACT
{
  /*!
  \brief Analytical solutions for 2D elasticity problems

  \param pos (in)       : position where analytical solution is evaluated
  \param uanalyt (out)  : analytical displacement solution
  \param epsanalyt (out): analytical strain solution
  */
  void AnalyticalSolutions2D(const LINALG::Matrix<2, 1>& pos, LINALG::Matrix<2, 1>& uanalyt,
      LINALG::Matrix<4, 1>& epsanalyt, LINALG::Matrix<2, 2>& derivanalyt);

  /*!
  \brief Analytical solutions for 3D elasticity problems

  \param pos (in)       : position where analytical solution is evaluated
  \param uanalyt (out)  : analytical displacement solution
  \param epsanalyt (out): analytical strain solution
  */
  void AnalyticalSolutions3D(const LINALG::Matrix<3, 1>& pos, LINALG::Matrix<3, 1>& uanalyt,
      LINALG::Matrix<6, 1>& epsanalyt, LINALG::Matrix<3, 3>& derivanalyt);

}  // namespace CONTACT

#endif  // #ifndef CONTACT_ANALYTICAL_H
