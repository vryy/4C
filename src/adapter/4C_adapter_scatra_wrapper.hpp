// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_ADAPTER_SCATRA_WRAPPER_HPP
#define FOUR_C_ADAPTER_SCATRA_WRAPPER_HPP


#include "4C_config.hpp"

#include "4C_adapter_scatra_interface.hpp"
#include "4C_linalg_vector.hpp"
#include "4C_scatra_timint_meshtying_strategy_base.hpp"

#include <memory>

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
    explicit AdapterScatraWrapper(std::shared_ptr<ScatraInterface> scatra);

    /// compute contribution of mechanical state to eq. system
    virtual void evaluate_additional_solution_depending_models(
        std::shared_ptr<Core::LinAlg::SparseOperator> systemmatrix,  //!< system matrix
        std::shared_ptr<Core::LinAlg::Vector<double>> rhs            //!< rhs vector
    );


   protected:
    //! return discretization
    std::shared_ptr<Core::FE::Discretization> discretization() const override
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
    std::shared_ptr<Core::LinAlg::Vector<double>> get_neumann_loads_ptr() override
    {
      return scatra_timint_->get_neumann_loads_ptr();
    };

    //! return meshtying strategy (includes standard case without meshtying)
    const std::shared_ptr<ScaTra::MeshtyingStrategyBase>& strategy() const override
    {
      return scatra_timint_->strategy();
    };

    //! return scalar field phi at time n
    std::shared_ptr<Core::LinAlg::Vector<double>> phin() override { return scatra_timint_->phin(); }

   private:
    std::shared_ptr<ScatraInterface> scatra_timint_;  ///< underlying structural time integration

  };  // class AdapterScatraWrapper
}  // namespace Adapter


FOUR_C_NAMESPACE_CLOSE

#endif
