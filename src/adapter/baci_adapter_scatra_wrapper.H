/*----------------------------------------------------------------------*/
/*! \file
\brief Wrapper for the scatra time integrator.
\level 1
 */
/*----------------------------------------------------------------------*/


#ifndef BACI_ADAPTER_SCATRA_WRAPPER_H
#define BACI_ADAPTER_SCATRA_WRAPPER_H


#include "baci_adapter_scatra_interface.H"
#include "baci_scatra_timint_meshtying_strategy_base.H"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>


// forward declarations
namespace CORE::LINALG
{
  class SparseOperator;
}

namespace ADAPTER
{
  class AdapterScatraWrapper : public ScatraInterface
  {
   public:
    /// constructor
    explicit AdapterScatraWrapper(Teuchos::RCP<ScatraInterface> scatra);

    /// compute contribution of mechanical state to eq. system
    virtual void EvaluateAdditionalSolutionDependingModels(
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,  //!< system matrix
        Teuchos::RCP<Epetra_Vector> rhs                           //!< rhs vector
    );


   protected:
    //! return discretization
    Teuchos::RCP<DRT::Discretization> Discretization() const override
    {
      return scatra_timint_->Discretization();
    };

    //! add parameters specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override
    {
      scatra_timint_->AddTimeIntegrationSpecificVectors(forcedincrementalsolver);
    };

    //! return number of dofset associated with displacement dofs
    virtual int NdsDisp() { return scatra_timint_->NdsDisp(); };

    /// return rcp ptr to neumann loads vector
    Teuchos::RCP<Epetra_Vector> GetNeumannLoadsPtr() override
    {
      return scatra_timint_->GetNeumannLoadsPtr();
    };

    //! return meshtying strategy (includes standard case without meshtying)
    const Teuchos::RCP<SCATRA::MeshtyingStrategyBase>& Strategy() const override
    {
      return scatra_timint_->Strategy();
    };

    //! return scalar field phi at time n
    Teuchos::RCP<Epetra_Vector> Phin() override { return scatra_timint_->Phin(); }

   private:
    Teuchos::RCP<ScatraInterface> scatra_timint_;  ///< underlying structural time integration

  };  // class AdapterScatraWrapper
}  // namespace ADAPTER


#endif  // ADAPTER_SCATRA_WRAPPER_H
