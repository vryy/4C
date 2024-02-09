/*-----------------------------------------------------------------------*/
/*! \file
\brief A set of analytical solutions for convergence analysis of contact/meshtying methods

\level 2

*/
/*-----------------------------------------------------------------------*/
#ifndef BACI_CONTACT_ANALYTICAL_HPP
#define BACI_CONTACT_ANALYTICAL_HPP

#include "baci_config.hpp"

#include "baci_linalg_fixedsizematrix.hpp"

BACI_NAMESPACE_OPEN

namespace CONTACT
{
  /*!
  \brief Analytical solutions for 2D elasticity problems

  \param pos (in)       : position where analytical solution is evaluated
  \param uanalyt (out)  : analytical displacement solution
  \param epsanalyt (out): analytical strain solution
  */
  void AnalyticalSolutions2D(const CORE::LINALG::Matrix<2, 1>& pos,
      CORE::LINALG::Matrix<2, 1>& uanalyt, CORE::LINALG::Matrix<4, 1>& epsanalyt,
      CORE::LINALG::Matrix<2, 2>& derivanalyt);

  /*!
  \brief Analytical solutions for 3D elasticity problems

  \param pos (in)       : position where analytical solution is evaluated
  \param uanalyt (out)  : analytical displacement solution
  \param epsanalyt (out): analytical strain solution
  */
  void AnalyticalSolutions3D(const CORE::LINALG::Matrix<3, 1>& pos,
      CORE::LINALG::Matrix<3, 1>& uanalyt, CORE::LINALG::Matrix<6, 1>& epsanalyt,
      CORE::LINALG::Matrix<3, 3>& derivanalyt);

}  // namespace CONTACT

BACI_NAMESPACE_CLOSE

#endif  // CONTACT_ANALYTICAL_H
