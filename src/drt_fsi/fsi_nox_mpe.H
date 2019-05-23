/*----------------------------------------------------------------------*/
/*!

\brief Solve FSI problems using minimal polynomial vector extrapolation

\maintainer Matthias Mayr

\level 1

*/
/*----------------------------------------------------------------------*/

#ifndef FSI_NOX_MPE_H
#define FSI_NOX_MPE_H

#include <NOX_Direction_Generic.H>  // base class
#include <NOX_Utils.H>
#include <NOX_GlobalData.H>
#include <Teuchos_RCP.hpp>
#include <NOX_Direction_UserDefinedFactory.H>

namespace NOX
{
  namespace FSI
  {
    //! %Calculates MPE.
    /*!

      Calculates the update direction by means of the Minimal Polynomial
      vector extrapolation method as proposed by Avram Sidi. That is we
      calcualte a krylov subspace based on the residuals here and
      approximate the solution inside this subspace. The size of the
      subspace is given by the user. In each cycle there is one operator
      evaluation, so this is an expensive direction search.

      This implements the inner cycle of the MPE that Prof Sidi
      proposes. Repeated calls by the outer nonlinear NOX loop are
      required.

      <h2>Possible extensions</h2>

      Just ideas. This may or may not work.

      - reuse the krylov subspace in the next loop

      - adjust the tolerance of the inner loop not to calculate
        digits that are not used later on

      - Use Aitken line search in addition to this extrapolation.

      \note It seems that there is no point in applying MPE to a vector
      sequence already accelerated with Irons&Tuck's Aitken version.

    <h2>Parameters</h2>

      - "kmax" - size of krylov subspace, equals number of inner
                 iterations (defaults to 10)

      - "omega" - fixed relaxation parameter, required to avoid divergence
                  (defaults to 1.0)

      - "Tolerance" - Tolerance test for residual in krylov subspace
                      (defaults to 1e-1)

    <h2>References</h2>

    Avram Sidi: Efficient implementation of minimal polynomial and reduced
    rank extrapolation methods. Journal of Computational and Applied
    Mathematics 36(3): 305-337, 1991.
    doi:10.1016/0377-0427(91)90013-A
    */
    class MinimalPolynomial : public NOX::Direction::Generic
    {
     public:
      //! Constructor
      MinimalPolynomial(const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params);

      //! Destructor
      ~MinimalPolynomial();

      // derived
      bool reset(const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params);

      // derived
      bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
          const NOX::Solver::Generic& solver);

      // derived
      virtual bool compute(NOX::Abstract::Vector& dir, NOX::Abstract::Group& grp,
          const NOX::Solver::LineSearchBased& solver);

     private:
      //! Print error message and throw error
      void throwError(const std::string& functionName, const std::string& errorMsg);

     private:
      //! Printing Utils
      Teuchos::RCP<NOX::Utils> utils_;

      int kmax_;
      double omega_;
      double eps_;
      bool mpe_;
    };

    /// simple factory that creates MPE direction object
    class MinimalPolynomialFactory : public NOX::Direction::UserDefinedFactory
    {
     public:
      virtual Teuchos::RCP<NOX::Direction::Generic> buildDirection(
          const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
      {
        return Teuchos::rcp(new MinimalPolynomial(gd->getUtils(), params));
      }
    };

  }  // namespace FSI
}  // namespace NOX

#endif
