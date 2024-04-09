/*-----------------------------------------------------------*/
/*! \file

\brief Basic time integration driver for Low Mach number problems


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FOUR_C_FLUID_TIMINT_LOMA_HPP
#define FOUR_C_FLUID_TIMINT_LOMA_HPP


#include "baci_config.hpp"

#include "baci_fluid_implicit_integration.hpp"

BACI_NAMESPACE_OPEN


namespace FLD
{
  class TimIntLoma : public virtual FluidImplicitTimeInt
  {
    friend class TimIntGenAlpha;

   public:
    /// Standard Constructor
    TimIntLoma(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);


    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief set scalar fields within outer iteration loop

    */
    void SetLomaIterScalarFields(Teuchos::RCP<const Epetra_Vector> scalaraf,
        Teuchos::RCP<const Epetra_Vector> scalaram, Teuchos::RCP<const Epetra_Vector> scalardtam,
        Teuchos::RCP<const Epetra_Vector> fsscalaraf, const double thermpressaf,
        const double thermpressam, const double thermpressdtaf, const double thermpressdtam,
        Teuchos::RCP<DRT::Discretization> scatradis) override;

    /*!
    \brief set scalar fields

    */
    void SetScalarFields(Teuchos::RCP<const Epetra_Vector> scalarnp, const double thermpressnp,
        Teuchos::RCP<const Epetra_Vector> scatraresidual,
        Teuchos::RCP<DRT::Discretization> scatradis, const int whichscalar = -1) override;

    /*!
    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element*/
    void SetElementCustomParameter();

    /*!
    \brief print turbulence model

    */
    void PrintTurbulenceModel() override;

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void SetCustomEleParamsApplyNonlinearBoundaryConditions(
        Teuchos::ParameterList& eleparams) override;

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams) override;

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void SetCustomEleParamsLinearRelaxationSolve(Teuchos::ParameterList& eleparams) override;

    /*!
    \brief Call statistics manager (special case in TimIntLoma)

    */
    void CallStatisticsManager() override;

    /// prepare AVM3-based scale separation
    void AVM3Preparation() override;

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


BACI_NAMESPACE_CLOSE

#endif
