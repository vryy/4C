/*-----------------------------------------------------------*/
/*! \file

\brief Basic time integration driver for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_LOMA_HPP
#define FOUR_C_FLUID_TIMINT_LOMA_HPP


#include "4C_config.hpp"

#include "4C_fluid_implicit_integration.hpp"

FOUR_C_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLoma : public virtual FluidImplicitTimeInt
  {
    friend class TimIntGenAlpha;

   public:
    /// Standard Constructor
    TimIntLoma(const Teuchos::RCP<Core::FE::Discretization>& actdis,
        const Teuchos::RCP<Core::LinAlg::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<Core::IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void init() override;

    /*!
    \brief set scalar fields within outer iteration loop

    */
    void set_loma_iter_scalar_fields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
        const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
        Teuchos::RCP<Core::FE::Discretization> scatradis) override;

    /*!
    \brief set scalar fields

    */
    void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp, const double thermpressnp,
        Teuchos::RCP<const Epetra_Vector> scatraresidual,
        Teuchos::RCP<Core::FE::Discretization> scatradis, const int whichscalar = -1) override;

    /*!
    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element*/
    void set_element_custom_parameter();

    /*!
    \brief print turbulence model

    */
    void print_turbulence_model() override;

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void set_custom_ele_params_apply_nonlinear_boundary_conditions(
        Teuchos::ParameterList& eleparams) override;

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void set_custom_ele_params_assemble_mat_and_rhs(Teuchos::ParameterList& eleparams) override;

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void set_custom_ele_params_linear_relaxation_solve(Teuchos::ParameterList& eleparams) override;

    /*!
    \brief Call statistics manager (special case in TimIntLoma)

    */
    void call_statistics_manager() override;

    /// prepare AVM3-based scale separation
    void av_m3_preparation() override;

    /*!
    \brief return thermpressaf_ in TimIntLoma

    */
    double ReturnThermpressaf() override { return thermpressaf_; }


   protected:
   private:
    /// for low-Mach-number flow solver: thermodynamic pressure at n+alpha_F/n+1
    /// and at n+alpha_M/n as well as its time derivative at n+alpha_F/n+1 and n+alpha_M/n
    double thermpressaf_;
    double thermpressam_;
    double thermpressdtaf_;
    double thermpressdtam_;


  };  // class TimIntLoma

}  // namespace FLD


FOUR_C_NAMESPACE_CLOSE

#endif
