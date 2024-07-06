/*----------------------------------------------------------------------*/
/*! \file
\brief scatra time integration for loma
\level 2
 *------------------------------------------------------------------------------------------------*/

#ifndef FOUR_C_SCATRA_TIMINT_LOMA_HPP
#define FOUR_C_SCATRA_TIMINT_LOMA_HPP

#include "4C_config.hpp"

#include "4C_scatra_timint_implicit.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_ParameterList.hpp>
#include <Teuchos_RCP.hpp>

FOUR_C_NAMESPACE_OPEN


namespace ScaTra
{
  class ScaTraTimIntLoma : public virtual ScaTraTimIntImpl
  {
   public:
    /// Standard Constructor
    ScaTraTimIntLoma(Teuchos::RCP<Core::FE::Discretization> dis,
        Teuchos::RCP<Core::LinAlg::Solver> solver, Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> sctratimintparams,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<Core::IO::DiscretizationWriter> output);

    /*========================================================================*/
    //! @name Preconditioning
    /*========================================================================*/

    void setup_splitter() override;

    // -----------------------------------------------------------------
    // general methods
    // -----------------------------------------------------------------

    /// initialize algorithm
    void init() override;

    /// initialize algorithm
    void setup() override;

    //! set initial thermodynamic pressure
    void set_initial_therm_pressure();

    //! predict thermodynamic pressure and time derivative
    virtual void predict_therm_pressure() = 0;

    //! compute initial total mass in domain
    void compute_initial_mass();

    //! compute thermodynamic pressure and time derivative
    virtual void compute_therm_pressure() = 0;

    //! compute thermodyn. press. from mass cons. in domain
    void compute_therm_pressure_from_mass_cons();

    //! compute values of thermodynamic pressure at intermediate time steps
    //! (required for generalized-alpha)
    virtual void compute_therm_pressure_intermediate_values() = 0;

    //!  compute time derivative of thermodynamic pressure after solution
    virtual void compute_therm_pressure_time_derivative() = 0;

    //! update thermodynamic pressure and time derivative
    virtual void update_therm_pressure() = 0;

    //! return thermo. press. at time step n
    double therm_press_n() const { return thermpressn_; }

    //! return thermo. press. at time step n+1
    double therm_press_np() const { return thermpressnp_; }

    //! return thermo. press. at time step n+alpha_F
    virtual double therm_press_af() = 0;

    //! return thermo. press. at time step n+alpha_M
    virtual double therm_press_am() = 0;

    //! return time der. of thermo. press. at time step n+1
    double therm_press_dt_np() const { return thermpressdtnp_; }

    //! return time derivative of thermo. press. at time step n+alpha_F
    virtual double therm_press_dt_af() = 0;

    //! return time derivative of thermo. press. at time step n+alpha_M
    virtual double therm_press_dt_am() = 0;

   protected:
    /*!
     * @brief add parameters depending on the problem, i.e., loma, level-set, ...
     *
     * @param params parameter list
     */
    void add_problem_specific_parameters_and_vectors(Teuchos::ParameterList& params) override;

    virtual void add_therm_press_to_parameter_list(
        Teuchos::ParameterList& params  //!< parameter list
        ) = 0;

    //! the parameter list for loma problems
    Teuchos::RCP<Teuchos::ParameterList> lomaparams_;

    //! initial mass in domain
    double initialmass_;

    //! thermodynamic pressure at n
    double thermpressn_;
    //! thermodynamic pressure at n+1
    double thermpressnp_;

    //! time deriv. of thermodynamic pressure at n
    double thermpressdtn_;
    //! time deriv. of thermodynamic pressure at n+1
    double thermpressdtnp_;
  };
}  // namespace ScaTra

FOUR_C_NAMESPACE_CLOSE

#endif
