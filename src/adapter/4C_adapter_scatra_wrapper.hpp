/*----------------------------------------------------------------------*/
/*! \file
\brief Wrapper for the scatra time integrator.
\level 1
 */
/*----------------------------------------------------------------------*/


#ifndef FOUR_C_ADAPTER_SCATRA_WRAPPER_HPP
#define FOUR_C_ADAPTER_SCATRA_WRAPPER_HPP


#include "4C_config.hpp"

#include "4C_adapter_scatra_interface.hpp"
#include "4C_scatra_timint_meshtying_strategy_base.hpp"

#include <Epetra_Vector.h>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


// forward declarations
namespace Core::LinAlg
{
  class SparseOperator;
}

namespace Adapter
{
  class AdapterScatraWrapper : public ScatraInterface
  {
   public:
    /// constructor
    explicit AdapterScatraWrapper(Teuchos::RCP<ScatraInterface> scatra);

    /// compute contribution of mechanical state to eq. system
    virtual void evaluate_additional_solution_depending_models(
        Teuchos::RCP<Core::LinAlg::SparseOperator> systemmatrix,  //!< system matrix
        Teuchos::RCP<Epetra_Vector> rhs                           //!< rhs vector
    );


   protected:
    //! return discretization
    Teuchos::RCP<Core::FE::Discretization> discretization() const override
    {
      return scatra_timint_->discretization();
    };

    //! add parameters specific for time-integration scheme
    void add_time_integration_specific_vectors(bool forcedincrementalsolver = false) override
    {
      scatra_timint_->add_time_integration_specific_vectors(forcedincrementalsolver);
    };

    //! return number of dofset associated with displacement dofs
    virtual int nds_disp() { return scatra_timint_->nds_disp(); };

    /// return rcp ptr to neumann loads vector
    Teuchos::RCP<Epetra_Vector> get_neumann_loads_ptr() override
    {
      return scatra_timint_->get_neumann_loads_ptr();
    };

    //! return meshtying strategy (includes standard case without meshtying)
    const Teuchos::RCP<ScaTra::MeshtyingStrategyBase>& strategy() const override
    {
      return scatra_timint_->strategy();
    };

    //! return scalar field phi at time n
    Teuchos::RCP<Epetra_Vector> phin() override { return scatra_timint_->phin(); }

   private:
    Teuchos::RCP<ScatraInterface> scatra_timint_;  ///< underlying structural time integration

  };  // class AdapterScatraWrapper
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
