/*!----------------------------------------------------------------------
\file fsi_nox_fixpoint.H

\brief Calculates the fix point direction.

\maintainer Andreas Rauch

\level 1

*----------------------------------------------------------------------*/

#ifndef FSI_NOX_FIXPOINT_H
#define FSI_NOX_FIXPOINT_H

#include <NOX_Direction_Generic.H>  // base class
#include <NOX_Utils.H>
#include <NOX_GlobalData.H>
#include <Teuchos_RCP.hpp>
#include <NOX_Direction_UserDefinedFactory.H>


/// NOX: the nonlinear object solver from trilinos with local extensions
namespace NOX
{
  /// NOX::FSI: local extensions to NOX, some of them special to FSI
  namespace FSI
  {
    //! %Calculates the fix point direction.
    /*!

      Calculates the direction
      \f[
      d = f(x)
      \f]

      This is the old (established, trusted, whatever you want) way to
      solve the nonlinear FSI interface equations. The residuum \f$ f(x)
      \f$ is
      \f[
      f(x) = g = d_{\Gamma,i+1} - d_{\Gamma,i}
      \f]
      in the FSI context.

      To be used with relaxation (line search).

    <B>Parameters</B>

    None.

    */
    class FixPoint : public NOX::Direction::Generic
    {
     public:
      //! Constructor
      FixPoint(const Teuchos::RCP<NOX::Utils>& utils, Teuchos::ParameterList& params);

      //! Destructor
      ~FixPoint();

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
    };

    /// simple factory that creates fix point direction object
    class FixPointFactory : public NOX::Direction::UserDefinedFactory
    {
     public:
      virtual Teuchos::RCP<NOX::Direction::Generic> buildDirection(
          const Teuchos::RCP<NOX::GlobalData>& gd, Teuchos::ParameterList& params) const
      {
        return Teuchos::rcp(new FixPoint(gd->getUtils(), params));
      }
    };

  }  // namespace FSI
}  // namespace NOX

#endif
