/*----------------------------------------------------------------------*/
/*! \file

\brief A class that wraps Teuchos::SerialDenseVector

\level 0
*/
/*----------------------------------------------------------------------*/

#include "baci_linalg_serialdensevector.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Compute vector 2-norm                                               |
 *----------------------------------------------------------------------*/
double CORE::LINALG::Norm2(const CORE::LINALG::SerialDenseVector& v) { return v.normFrobenius(); }

/*----------------------------------------------------------------------*
 |  b = alpha*a + beta*b                                                |
 *----------------------------------------------------------------------*/
void CORE::LINALG::Update(double alpha, const CORE::LINALG::SerialDenseVector& a, double beta,
    CORE::LINALG::SerialDenseVector& b)
{
  b.scale(beta);
  CORE::LINALG::SerialDenseVector acopy(a);
  acopy.scale(alpha);
  b += acopy;
}

FOUR_C_NAMESPACE_CLOSE
