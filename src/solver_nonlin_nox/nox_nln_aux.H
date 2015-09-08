/*-----------------------------------------------------------*/
/*!
\file nox_nln_aux.H

\maintainer Michael Hiermeier

\date Jul 31, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_AUX_H_
#define NOX_NLN_AUX_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_statustest_factory.H"

// forward declaration
class Epetra_Vector;

namespace NOX {
namespace StatusTest {
  class Generic;
} // namespace StatusTest
namespace NLN {
namespace AUX {
  /*! \brief Calculate the root mean square for the NOX status test

      $ \delta_{rms} = \sqrt{\frac{1}{N} \sum\limits_{i=1}^{N} \left( \frac{x_{i}^{k} - x_{i}^{k-1}}{\mathcal{RTOL} | x_{i}^{k-1} | + \mathcal{ATOL}} \right)} $

      \param atol  : absolute tolerance
      \param rtol  : relative tolerance
      \param xnew  : new / current iterate $x_{i}^{k}$
      \param xincr : current step increment $x_{i}^{k} - x_{i}^{k-1}$
  */
  double RootMeanSquareNorm(const double& atol, const double& rtol,
      Teuchos::RCP<const Epetra_Vector> xnew, Teuchos::RCP<const Epetra_Vector> xincr,
      const bool& disable_implicit_weighting=false);

  //! Do a recursive search for a nox_nln_statustest_normf test in the StatusTest object
  //! and return the class variable value of the desired quantity.
  double GetNormFClassVariable(
      const NOX::StatusTest::Generic& test,
      const NOX::NLN::StatusTest::QuantityType& qType,
      const std::string& classVariableName);

} // namepsace AUX
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_AUX_H_ */
