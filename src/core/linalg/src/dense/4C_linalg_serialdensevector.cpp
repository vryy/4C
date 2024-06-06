/*----------------------------------------------------------------------*/
/*! \file

\brief A class that wraps Teuchos::SerialDenseVector

\level 0
*/
/*----------------------------------------------------------------------*/

#include "4C_linalg_serialdensevector.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Compute vector 2-norm                                               |
 *----------------------------------------------------------------------*/
double Core::LinAlg::Norm2(const Core::LinAlg::SerialDenseVector& v) { return v.normFrobenius(); }

/*----------------------------------------------------------------------*
 |  b = alpha*a + beta*b                                                |
 *----------------------------------------------------------------------*/
void Core::LinAlg::Update(double alpha, const Core::LinAlg::SerialDenseVector& a, double beta,
    Core::LinAlg::SerialDenseVector& b)
{
  b.scale(beta);
  Core::LinAlg::SerialDenseVector acopy(a);
  acopy.scale(alpha);
  b += acopy;
}

FOUR_C_NAMESPACE_CLOSE
