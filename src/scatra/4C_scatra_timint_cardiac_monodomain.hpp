/*----------------------------------------------------------------------*/
/*! \file

\brief scatra time integration for cardiac monodomain

\level 2


 *------------------------------------------------------------------------------------------------*/

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

    void output_state() override;

    //! Set ep-specific parameters
    void set_element_specific_sca_tra_parameters(Teuchos::ParameterList& eleparams) const override;

    /*========================================================================*/
    //! @name electrophysiology
    /*========================================================================*/

    //! activation_time at times n+1
    Teuchos::RCP<Epetra_Vector> activation_time_np_;

    //! activation threshold for postprocessing
    double activation_threshold_;

    //! maximum expected number of material internal state variables
    int nb_max_mat_int_state_vars_;

    //! material internal state at times n+1
    Teuchos::RCP<Epetra_MultiVector> material_internal_state_np_;

    //! one component of the material internal state at times n+1 (for separated postprocessing)
    Teuchos::RCP<Epetra_Vector> material_internal_state_np_component_;

    //! maximum expected number of material ionic currents
    int nb_max_mat_ionic_currents_;

    //! material ionic currents at times n+1
    Teuchos::RCP<Epetra_MultiVector> material_ionic_currents_np_;

    //! one component of the material ionic currents at times n+1 (for separated postprocessing)
    Teuchos::RCP<Epetra_Vector> material_ionic_currents_np_component_;

    //! parameter list
    const Teuchos::RCP<Teuchos::ParameterList> ep_params_;
  };

};  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
