/*!----------------------------------------------------------------------
\file contact_analytical.H

\brief A set of analytical solutions for convergence analysis of contact/meshtying methods

<pre>
-------------------------------------------------------------------------
                        BACI Contact library
            Copyright (2008) Technical University of Munich

Under terms of contract T004.008.000 there is a non-exclusive license for use
of this work by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library is proprietary software. It must not be published, distributed,
copied or altered in any form or any media without written permission
of the copyright holder. It may be used under terms and conditions of the
above mentioned license by or on behalf of Rolls-Royce Ltd & Co KG, Germany.

This library contains and makes use of software copyrighted by Sandia Corporation
and distributed under LGPL licence. Licensing does not apply to this or any
other third party software used here.

Questions? Contact Dr. Michael W. Gee (gee@lnm.mw.tum.de)
                   or
                   Prof. Dr. Wolfgang A. Wall (wall@lnm.mw.tum.de)

http://www.lnm.mw.tum.de

-------------------------------------------------------------------------
</pre>

<pre>
Maintainer: Alexander Popp
            popp@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15238
</pre>

*----------------------------------------------------------------------*/
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
  void AnalyticalSolutions2D(const LINALG::Matrix<2,1>& pos,
                             LINALG::Matrix<2,1>& uanalyt,
                             LINALG::Matrix<4,1>& epsanalyt,
                             LINALG::Matrix<2,2>& derivanalyt);

  /*!
  \brief Analytical solutions for 3D elasticity problems

  \param pos (in)       : position where analytical solution is evaluated
  \param uanalyt (out)  : analytical displacement solution
  \param epsanalyt (out): analytical strain solution
  */
  void AnalyticalSolutions3D(const LINALG::Matrix<3,1>& pos,
                             LINALG::Matrix<3,1>& uanalyt,
                             LINALG::Matrix<6,1>& epsanalyt,
                             LINALG::Matrix<3,3>& derivanalyt);

} // namespace CONTACT

#endif  // #ifndef CONTACT_ANALYTICAL_H
