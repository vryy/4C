/*----------------------------------------------------------------------*/
/*! \file

\brief Linear system for solving FSI problems with NOX

\level 2

*/
/*----------------------------------------------------------------------*/

#ifndef FOUR_C_FSI_MONOLITHIC_LINEARSYSTEM_HPP
#define FOUR_C_FSI_MONOLITHIC_LINEARSYSTEM_HPP


#include "baci_config.hpp"

#include <AztecOO_StatusTestCombo.h>
#include <AztecOO_StatusTestMaxIters.h>
#include <AztecOO_StatusTestResNorm.h>
#include <NOX_Epetra_LinearSystem_AztecOO.H>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

BACI_NAMESPACE_OPEN


namespace FSI
{
  /// monolithic linear system
  /*!
    Implements ::NOX::Epetra::LinearSystemAztecOO.
    Currently, the only purpose is to introduce a combo status test
    to the aztec solver that allows for
    L2-norm, Linf-norm and MaxIters

    \author m.gee
    \date 03/08
   */
  class MonolithicLinearSystem : public ::NOX::Epetra::LinearSystemAztecOO
  {
   public:
    /// create
    explicit MonolithicLinearSystem(Teuchos::ParameterList& printingParams,
        Teuchos::ParameterList& linearSolverParams,
        const Teuchos::RCP<::NOX::Epetra::Interface::Jacobian>& iJac,
        const Teuchos::RCP<Epetra_Operator>& J,
        const Teuchos::RCP<::NOX::Epetra::Interface::Preconditioner>& iPrec,
        const Teuchos::RCP<Epetra_Operator>& M, const ::NOX::Epetra::Vector& cloneVector,
        const Teuchos::RCP<::NOX::Epetra::Scaling> scalingObject = Teuchos::null);


    /// overload the solve to incorporate our own convergence test
    bool applyJacobianInverse(Teuchos::ParameterList& linearSolverParams,
        const ::NOX::Epetra::Vector& input, ::NOX::Epetra::Vector& result) override;

   protected:
    void softreset(Teuchos::ParameterList& linearSolverParams);

    //! copy of the linear solver parameters passed in constructor
    Teuchos::ParameterList lsparams_;

    Teuchos::RCP<AztecOO_StatusTestMaxIters> aztest_maxiter_;
    //! an aztec 2-norm convergence test
    Teuchos::RCP<AztecOO_StatusTestResNorm> aztest_norm2_;
    //! an aztec inf-norm convergence test
    Teuchos::RCP<AztecOO_StatusTestResNorm> aztest_norminf_;
    // ! combination of the above tests
    Teuchos::RCP<AztecOO_StatusTestCombo> aztest_combo1_;
    Teuchos::RCP<AztecOO_StatusTestCombo> aztest_combo2_;

   private:
  };

}  // namespace FSI


BACI_NAMESPACE_CLOSE

#endif
