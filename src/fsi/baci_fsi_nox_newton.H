/*----------------------------------------------------------------------*/
/*! \file

\brief NOX Newton direction with adaptive linear solver tolerance for FSI

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef BACI_FSI_NOX_NEWTON_H
#define BACI_FSI_NOX_NEWTON_H

#include "baci_config.H"

#include "baci_inpar_fsi.H"

#include <NOX_Direction_Newton.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

#include <vector>

BACI_NAMESPACE_OPEN

namespace NOX
{
  namespace FSI
  {
    /// NOX Newton direction with adaptive linear solver tolerance
    class Newton : public ::NOX::Direction::Newton
    {
     public:
      Newton(const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params);


      // derived
      bool reset(
          const Teuchos::RCP<::NOX::GlobalData>& gd, Teuchos::ParameterList& params) override;

      // derived
      bool compute(::NOX::Abstract::Vector& dir, ::NOX::Abstract::Group& grp,
          const ::NOX::Solver::Generic& solver) override;

      void Residual(double current, double desired);

     private:
      /// Printing Utilities
      Teuchos::RCP<::NOX::Utils> utils_;

      //! "Direction" sublist with parameters for the direction vector
      Teuchos::ParameterList* paramsPtr;

      /// nonlinear tolerance we strive to achieve
      double desirednlnres_;

      /// current nonlinear tolerance (what we gained after the last linear solve)
      double currentnlnres_;

      /// basic (unmodified) linear solver (AZ_r0) tolerance
      double plaintol_;

      /// improvement factor
      double better_;

      /// verbosity level of FSI algorithm
      INPAR::FSI::Verbosity verbosity_;

      std::vector<double> cresiduals_;
      std::vector<double> dresiduals_;
    };
  }  // namespace FSI
}  // namespace NOX

BACI_NAMESPACE_CLOSE

#endif
