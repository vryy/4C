/*-----------------------------------------------------------*/
/*!

\brief base class of time integration schemes for porous fluid

\maintainer Christoph Ager

\level 2

*/
/*-----------------------------------------------------------*/

#ifndef FLD_TIMINT_PORO_H
#define FLD_TIMINT_PORO_H


#include "fluidimplicitintegration.H"

namespace FLD
{
  class TimIntPoro : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntPoro(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntPoro();

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief call elements to calculate system matrix/rhs and assemble

    */
    virtual void AssembleMatAndRHS();

    /*!
    \brief read restart data
    */
    virtual void ReadRestart(int step);

    //! @name Set general parameter in class f3Parameter
    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element*/
    virtual void SetElementCustomParameter();

    /// set the initial porosity field
    virtual void SetInitialPorosityField(
        const INPAR::POROELAST::InitialField,  ///< type of initial field
        const int startfuncno                  ///< number of spatial function
    );

    /*!
    \brief update iterative increment

    */
    virtual void UpdateIterIncrementally(
        Teuchos::RCP<const Epetra_Vector> vel  //!< input residual velocities
    );

    /*!
    \brief update configuration and output to file/screen

    */
    virtual void Output();

    /*!
    \brief Do some poro-specific stuff in AssembleMatAndRHS

    */
    virtual void PoroIntUpdate();

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    virtual void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams);

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    virtual void TimIntCalculateAcceleration();

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    virtual void SetElementGeneralFluidParameter();

    /*!

    \brief parameter (fix over all time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    virtual void SetElementTurbulenceParameters();


   protected:
    /// initial porosity (poroelasticity)
    Teuchos::RCP<Epetra_Vector> initporosityfield_;

   private:
  };  // class TimIntPoro

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_PORO_H
