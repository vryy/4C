// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#ifndef FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_HPP
#define FOUR_C_SCATRA_TIMINT_CARDIAC_MONODOMAIN_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN


/*==========================================================================*/
// forward declarations
/*==========================================================================*/


namespace ScaTra
{
  class TimIntCardiacMonodomain : public virtual ScaTraTimIntImpl
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimIntCardiacMonodomain(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);


    //! setup algorithm
    void setup() override;

    //! time update of time-dependent materials
    virtual void element_material_time_update();

    void collect_runtime_output_data() override;

    //! Set ep-specific parameters
    void set_element_specific_scatra_parameters(Teuchos::ParameterList& eleparams) const override;

    void write_restart() const override;

    /*========================================================================*/
    //! @name electrophysiology
    /*========================================================================*/

    //! activation_time at times n+1
    Teuchos::RCP<Core::LinAlg::Vector<double>> activation_time_np_;

    //! activation threshold for postprocessing
    double activation_threshold_;

    //! maximum expected number of material internal state variables
    int nb_max_mat_int_state_vars_;

    //! material internal state at times n+1
    Teuchos::RCP<Core::LinAlg::MultiVector<double>> material_internal_state_np_;

    //! one component of the material internal state at times n+1 (for separated postprocessing)
    Teuchos::RCP<Core::LinAlg::Vector<double>> material_internal_state_np_component_;

    //! maximum expected number of material ionic currents
    int nb_max_mat_ionic_currents_;

    //! material ionic currents at times n+1
    Teuchos::RCP<Core::LinAlg::MultiVector<double>> material_ionic_currents_np_;

    //! one component of the material ionic currents at times n+1 (for separated postprocessing)
    Teuchos::RCP<Core::LinAlg::Vector<double>> material_ionic_currents_np_component_;

    //! parameter list
    const Teuchos::RCP<Teuchos::ParameterList> ep_params_;
  };

};  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
