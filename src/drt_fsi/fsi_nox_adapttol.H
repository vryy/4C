#ifdef CCADISCRET

#ifndef FSI_NOX_ADAPTTOL_H
#define FSI_NOX_ADAPTTOL_H

#include <NOX_Abstract_PrePostOperator.H>
#include <NOX_Utils.H>

#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

namespace NOX
{
  namespace FSI
  {

    /// calculate adaptive tolerance for linear solver in Newton loop
    /*!

      This is a NOX version of the M. Gee hack.

     */
    class AdaptiveTolerance : public NOX::Abstract::PrePostOperator
    {
    public:
      AdaptiveTolerance(Teuchos::RCP<NOX::Utils> utils,
                        Teuchos::ParameterList& lsParams);

      virtual void runPreIterate(const NOX::Solver::Generic& solver);

      virtual void runPostIterate(const NOX::Solver::Generic& solver);

    private:

      /// utils for screen output
      Teuchos::RCP<NOX::Utils> utils_;

      /// linear solver parameters
      Teuchos::ParameterList& lsParams_;

      /// nonlinear tolerance we strive to achieve
      double desirednlnres_;

      /// current nonlinear tolerance (what we gained after the last linear solve)
      double currentnlnres_;

      /// basic (unmodified) linear solver (AZ_r0) tolerance
      double plaintol_;

      /// improvement factor
      double better_;
    };
  }
}

#endif
#endif
