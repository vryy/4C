/*-----------------------------------------------------------*/
/*! \file

\brief Basic fluid driver for stationary problems

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FLD_TIMINT_STAT_H
#define FLD_TIMINT_STAT_H


#include "fluidimplicitintegration.H"


namespace FLD
{
  class TimIntStationary : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntStationary(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntStationary();

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief Do time integration (time loop)

    */
    void TimeLoop();

    /*!
    \brief Set the part of the righthandside belonging to the last
           timestep for incompressible or low-Mach-number flow

       for low-Mach-number flow: distinguish momentum and continuity part
       (continuity part only meaningful for low-Mach-number flow)

       Stationary/af-generalized-alpha:

                     mom: hist_ = 0.0
                    (con: hist_ = 0.0)

       One-step-Theta:

                     mom: hist_ = veln_  + dt*(1-Theta)*accn_
                    (con: hist_ = densn_ + dt*(1-Theta)*densdtn_)

       BDF2: for constant time step:

                     mom: hist_ = 4/3 veln_  - 1/3 velnm_
                    (con: hist_ = 4/3 densn_ - 1/3 densnm_)


    */
    void SetOldPartOfRighthandside();

    /*!
    \brief Solve stationary problem

    */
    void SolveStationaryProblem();

    /*!
    \brief Set states in the time integration schemes: differs between GenAlpha and the others

    */
    void SetStateTimInt();

    /*!
    \brief Calculate time derivatives for
           stationary/one-step-theta/BDF2/af-generalized-alpha time integration
           for incompressible and low-Mach-number flow
    */
    void CalculateAcceleration(const Teuchos::RCP<const Epetra_Vector> velnp,  ///< velocity at n+1
        const Teuchos::RCP<const Epetra_Vector> veln,   ///< velocity at     n
        const Teuchos::RCP<const Epetra_Vector> velnm,  ///< velocity at     n-1
        const Teuchos::RCP<const Epetra_Vector> accn,   ///< acceleration at n-1
        const Teuchos::RCP<Epetra_Vector> accnp         ///< acceleration at n+1
    );

    /*!
    \brief Set gamma to a value

    */
    void SetGamma(Teuchos::ParameterList& eleparams);

    /*!
    \brief Scale separation

    */
    void Sep_Multiply();

    /*!
    \brief Output of filtered velocity

    */
    void OutputofFilteredVel(
        Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec);

    /*!

    \brief parameter (fix over a time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    void SetElementTimeParameter();

    /*!
    \brief return scheme-specific time integration parameter

    */
    virtual double TimIntParam() const;

    /*!
    \brief return scaling factor for the residual

    */
    double ResidualScaling() const { return 1.0; }

    /*!
    \brief velocity required for evaluation of related quantites required on element level

    */
    Teuchos::RCP<const Epetra_Vector> EvaluationVel() { return Teuchos::null; };

    /*!
    \brief treat turbulence models in AssembleMatAndRHS
    */
    virtual void TreatTurbulenceModels(Teuchos::ParameterList& eleparams);

    //! @name Time Step Size Adaptivity
    //@{

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const { return 1; }

    //! Give local order of accuracy of pressure part
    virtual int MethodOrderOfAccuracyPres() const { return 1; }

    //! Return linear error coefficient of velocity
    virtual double MethodLinErrCoeffVel() const { return 1.0; }

    //@}

   protected:
    /*!
    \brief break criterion for pseudo timeloop

    */
    bool NotFinished() { return step_ < stepmax_; }


   private:
  };  // class TimIntStationary

}  // namespace FLD

#endif  // #ifndef FLD_TIMINT_STAT_H
