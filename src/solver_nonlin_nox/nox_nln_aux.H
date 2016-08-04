/*-----------------------------------------------------------*/
/*!
\file nox_nln_aux.H

\brief Auxiliary methods.

\maintainer Michael Hiermeier

\date Jul 31, 2015

\level 3

*/
/*-----------------------------------------------------------*/

#ifndef NOX_NLN_AUX_H_
#define NOX_NLN_AUX_H_

#include <Teuchos_RCP.hpp>
#include "nox_nln_statustest_factory.H"
#include "nox_nln_enum_lists.H"

#include <NOX_Abstract_Vector.H>

// forward declaration
class Epetra_Vector;
class Epetra_Comm;
namespace LINALG {
  class Solver;
} // namespace LINALG
namespace NOX {
namespace StatusTest {
  class Generic;
} // namespace StatusTest
namespace NLN {
namespace AUX {

  /*! Set printing parameters
   *
   *  Note: The Yes/No tuples are translated to booleans! */
  void SetPrintingParameters(Teuchos::ParameterList& p_nox,
      const Epetra_Comm& comm);

  NOX::NLN::LinSystem::LinearSystemType GetLinearSystemType(
      const std::map<enum NOX::NLN::SolutionType,Teuchos::RCP<LINALG::Solver> >& linsolvers);

  /*! \brief Calculate the root mean square for the NOX status test
   *  \f[
   *    \delta_{rms} = \sqrt{\frac{1}{N} \sum\limits_{i=1}^{N} \left( \frac{x_{i}^{k} - x_{i}^{k-1}}{\mathcal{RTOL} | x_{i}^{k-1} | + \mathcal{ATOL}} \right)}
   *  \f]
   *
   *  \param atol  : absolute tolerance
   *  \param rtol  : relative tolerance
   *  \param xnew  : new / current iterate $x_{i}^{k}$
   *  \param xincr : current step increment $x_{i}^{k} - x_{i}^{k-1}$
   */
  double RootMeanSquareNorm(const double& atol, const double& rtol,
      Teuchos::RCP<const Epetra_Vector> xnew, Teuchos::RCP<const Epetra_Vector> xincr,
      const bool& disable_implicit_weighting=false);

  /*! \brief Do a recursive search for a NOX::NLN::StatusTest::NormWRMS object in the StatusTest object list
   *  and return the class variable value of the desired quantity.
   *
   * \param test              : StatusTest object which will be scanned.
   * \param qType             : Quantity type of the NormWRMS test which we are looking for.
   * \param classVariableName : Name of the class variable which will be returned. (Type: double) */
  double GetNormWRMSClassVariable(
      const NOX::StatusTest::Generic& test,
      const NOX::NLN::StatusTest::QuantityType& qType,
      const std::string& classVariableName);

  /*! \brief Do a recursive search for a NOX::NLN::StatusTest::NormF object in the StatusTest object list
   *  and return the class variable value of the desired quantity.
   *
   * \param test              : StatusTest object which will be scanned.
   * \param qType             : Quantity type of the NormF test which we are looking for.
   * \param classVariableName : Name of the class variable which will be returned. (Type: double) */
  double GetNormFClassVariable(
      const NOX::StatusTest::Generic& test,
      const NOX::NLN::StatusTest::QuantityType& qType,
      const std::string& classVariableName);

  /*! Do a recursive search for a <T> status test and the given quantity
   *
   *  True is returned as soon as a status test of type <T> is found, which
   *  holds the given quantity. */
  template <class T>
  bool IsQuantity(const NOX::StatusTest::Generic& test,
      const NOX::NLN::StatusTest::QuantityType& qtype);

  /*! \brief Do a recursive search for a <T> status test class and return the NormType of the given quantity.
   *
   *  If there are more than one status tests of the type <T> which hold the given quantity, the normtype of
   *  the first we can find, will be returned! */
  template <class T>
  int GetNormType(const NOX::StatusTest::Generic& test,
      const NOX::NLN::StatusTest::QuantityType& qtype);

  /*! \brief Do a recursive search for a <T> status test class and return its status.
   *
   * If more than one of the given status test objects is combined in a combination list,
   * the AND combination of the different status is returned. I.e. if one of the status
   * is unconverged, the return status is unconverged.
   * If we cannot find the given status test class, a default value of -100 is returned.
   *
   * \param test : StatusTest object which will be scanned.
   */
  template <class T>
  int GetOuterStatus(const NOX::StatusTest::Generic& test);

  /*! \brief Convert the quantity type to a solution type
   *
   * \param qtype : Quantity type which has to be converted.
   */
  enum NOX::NLN::SolutionType ConvertQuantityType2SolutionType(
      const enum NOX::NLN::StatusTest::QuantityType& qtype);

  /*! \brief Map norm type stl_string to norm type enum
   *
   * \param name : Name of the vector norm type.
   */
  enum NOX::Abstract::Vector::NormType String2NormType(
      const std::string& name);

} // namepsace AUX
} // namespace NLN
} // namespace NOX


#endif /* NOX_NLN_AUX_H_ */
