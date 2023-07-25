/*-----------------------------------------------------------*/
/*! \file

\brief base class of time integration schemes for porous fluid


\level 2

*/
/*-----------------------------------------------------------*/

#ifndef BACI_FLUID_TIMINT_PORO_H
#define BACI_FLUID_TIMINT_PORO_H


#include "baci_fluid_implicit_integration.H"

namespace FLD
{
  class TimIntPoro : public virtual FluidImplicitTimeInt
  {
   public:
    //! Standard Constructor
    TimIntPoro(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<CORE::LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /*!
    \brief initialization

    */
    void Init() override;

    /*!
    \brief call elements to calculate system matrix/rhs and assemble

    */
    void AssembleMatAndRHS() override;

    /*!
    \brief read restart data
    */
    void ReadRestart(int step) override;

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element*/
    virtual void SetElementCustomParameter();

    //! set the initial porosity field
    void SetInitialPorosityField(
        const INPAR::POROELAST::InitialField init,  //!< type of initial field
        const int startfuncno                       //!< number of spatial function
        ) override;

    /*!
    \brief update iterative increment

    */
    void UpdateIterIncrementally(
        Teuchos::RCP<const Epetra_Vector> vel  //!< input residual velocities
        ) override;

    /*!
    \brief update configuration and output to file/screen

    */
    void Output() override;

    /*!
    \brief Do some poro-specific stuff in AssembleMatAndRHS

    */
    virtual void PoroIntUpdate();

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams) override;

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void TimIntCalculateAcceleration() override;

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void SetElementGeneralFluidParameter() override;

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void SetElementTurbulenceParameters() override;


   protected:
    //! initial porosity (poroelasticity)
    Teuchos::RCP<Epetra_Vector> init_porosity_field_;

   private:
  };

}  // namespace FLD


#endif  // FLUID_TIMINT_PORO_H
