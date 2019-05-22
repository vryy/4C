/*-----------------------------------------------------------*/
/*!

\brief Basic HDG time-integration scheme

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FLD_TIMINT_HDG_H
#define FLD_TIMINT_HDG_H


#include "fluid_timint_genalpha.H"


namespace FLD
{
  /*!
  \brief time integration for HDG fluid (only gen-alpha implemented)

  \author kronbichler (kronbichler@lnm.mw.tum.de)
  */
  class TimIntHDG : public TimIntGenAlpha
  {
   public:
    /// Standard Constructor
    TimIntHDG(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntHDG() {}

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief Set theta_ to its value, dependent on integration method for GenAlpha and BDF2

    */
    virtual void SetTheta();

    /*!
    \brief do explicit predictor step to start nonlinear iteration from
           a better initial value
    */
    virtual void ExplicitPredictor();

    /*!
    \brief Set custom parameters in the respective time integration class (Loma, RedModels...)

    */
    virtual void SetCustomEleParamsAssembleMatAndRHS(Teuchos::ParameterList& eleparams);

    /*!
    \brief Set states in the time integration schemes: additional vectors for HDG

    */
    virtual void SetStateTimInt();

    /*!
    \brief Call discret_->ClearState() after assembly (HDG needs to read from state vectors...)

    */
    virtual void ClearStateAssembleMatAndRHS();

    /*!
    \brief Set the part of the right hand side belonging to the last
           time step for incompressible or low-Mach-number flow

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
    virtual void SetOldPartOfRighthandside();

    /*!
    \brief Update the solution after convergence of the nonlinear
           iteration. Current solution becomes old solution of next
           time step.
    */
    virtual void TimeUpdate();

    /*!
    \brief set initial flow field for analytical test problems

    */
    virtual void SetInitialFlowField(
        const INPAR::FLUID::InitialField initfield, const int startfuncno);

    /*!
    \brief calculate error between a analytical solution and the
           numerical solution of a test problems

    */
    virtual Teuchos::RCP<std::vector<double>> EvaluateErrorComparedToAnalyticalSol();

    /*!
    \brief Reset state vectors
     */
    virtual void Reset(bool completeReset = false, int numsteps = 1, int iter = -1);

    /*!
    \brief update configuration and output to file/screen

    */
    virtual void Output();

    /*!
    \brief accessor to interior velocity

    */
    virtual Teuchos::RCP<Epetra_Vector> ReturnIntVelnp() { return intvelnp_; }
    virtual Teuchos::RCP<Epetra_Vector> ReturnIntVeln() { return intveln_; }
    virtual Teuchos::RCP<Epetra_Vector> ReturnIntVelnm() { return intvelnm_; }


   protected:
    /// copy constructor
    TimIntHDG(const TimIntHDG& old);

    /*!
    \brief update acceleration for generalized-alpha time integration

    */
    virtual void GenAlphaUpdateAcceleration();

    /*!
    \brief compute values at intermediate time steps for gen.-alpha

    */
    virtual void GenAlphaIntermediateValues();

    //! @name velocity gradient, velocity and pressure at time n+1, n, n-1
    //!  and n+alpha_F for element interior in HDG
    Teuchos::RCP<Epetra_Vector> intvelnp_;
    Teuchos::RCP<Epetra_Vector> intveln_;
    Teuchos::RCP<Epetra_Vector> intvelnm_;
    Teuchos::RCP<Epetra_Vector> intvelaf_;
    //@}

    //! @name acceleration/(scalar time derivative) at time n+1, n and n+alpha_M/(n+alpha_M/n)
    //!  and n-1 for element interior in HDG
    //@{
    Teuchos::RCP<Epetra_Vector> intaccnp_;  ///< acceleration at time \f$t^{n+1}\f$
    Teuchos::RCP<Epetra_Vector> intaccn_;   ///< acceleration at time \f$t^{n}\f$
    Teuchos::RCP<Epetra_Vector> intaccam_;  ///< acceleration at time \f$t^{n+\alpha_M}\f$
    Teuchos::RCP<Epetra_Vector> intaccnm_;  ///< acceleration at time \f$t^{n-1}\f$
    //@}

    //! @name other HDG-specific auxiliary vectors, like velocity for output
    Teuchos::RCP<Epetra_MultiVector> interpolatedVelocity_;
    Teuchos::RCP<Epetra_Vector> interpolatedPressure_;
    //@}

    /*!
    \brief calculate intermediate solution

    */
    virtual void CalcIntermediateSolution();

    /*!
    \brief Initialize forcing

    */
    virtual void InitForcing();

   private:
    ///< Print stabilization details to screen. Do nothing here because we do not use stabilization
    virtual void PrintStabilizationDetails() const {}

    ///< time algorithm flag actually set (we internally reset it)
    INPAR::FLUID::TimeIntegrationScheme timealgoset_;

    ///< Keep track of whether we do the first assembly of a time step because we reconstruct the
    ///< local HDG solution as part of assembly
    bool firstAssembly_;

  };  // class TimIntHDG

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_HDG_H
