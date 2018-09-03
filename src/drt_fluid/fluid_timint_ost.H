/*!----------------------------------------------------------------------
\file fluid_timint_ost.H

\brief One-Step-Theta time-integration scheme

\level 2
<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLUID_TIMINT_OST_H
#define FLUID_TIMINT_OST_H


#include "fluidimplicitintegration.H"


namespace FLD
{
  class TimIntOneStepTheta : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntOneStepTheta(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntOneStepTheta();

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief Print information about current time step to screen (reimplementation)

    */
    virtual void PrintTimeStepInfo();

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
    virtual void SetOldPartOfRighthandside();

    /*!
    \brief Set states in the time integration schemes: differs between GenAlpha and the others

    */
    virtual void SetStateTimInt();

    /*!
    \brief Calculate time derivatives for
           stationary/one-step-theta/BDF2/af-generalized-alpha time integration
           for incompressible and low-Mach-number flow
    */
    virtual void CalculateAcceleration(
        const Teuchos::RCP<const Epetra_Vector> velnp,  ///< velocity at     n+1
        const Teuchos::RCP<const Epetra_Vector> veln,   ///< velocity at     n
        const Teuchos::RCP<const Epetra_Vector> velnm,  ///< velocity at     n-1
        const Teuchos::RCP<const Epetra_Vector> accn,   ///< acceleration at n-1
        const Teuchos::RCP<Epetra_Vector> accnp         ///< acceleration at n+1
    );

    /*!
    \brief Set gamma to a value

    */
    virtual void SetGamma(Teuchos::ParameterList& eleparams);

    /*!
    \brief Scale separation

    */
    virtual void Sep_Multiply();

    /*!
    \brief Output of filtered velocity

    */
    virtual void OutputofFilteredVel(
        Teuchos::RCP<Epetra_Vector> outvec, Teuchos::RCP<Epetra_Vector> fsoutvec);

    /*!

    \brief parameter (fix over a time step) are set in this method.
           Therefore, these parameter are accessible in the fluid element
           and in the fluid boundary element

    */
    virtual void SetElementTimeParameter();

    /*!
    \brief set theta if starting algorithm is chosen.

    */
    virtual void SetTheta();

    /*!
    \brief return scheme-specific time integration parameter

    */
    virtual double TimIntParam() const { return 0.0; }

    /*!
    \brief return scaling factor for the residual

    */
    virtual double ResidualScaling() const
    {
      if (params_->get<bool>("ost new"))
        return 1.0 / (dta_);
      else
        return 1.0 / (theta_ * dta_);
    }

    /*!
    \brief velocity required for evaluation of related quantites required on element level

    */
    virtual Teuchos::RCP<const Epetra_Vector> EvaluationVel() { return velnp_; };

    /*!
    \ apply external forces to the fluid

    */
    virtual void ApplyExternalForces(Teuchos::RCP<Epetra_MultiVector> fext);

    /*!
    \output of external forces for restart

    */
    virtual void OutputExternalForces();

    /*!
    \brief read restart data

    */
    virtual void ReadRestart(int step);

    /*!
    \ Update of external forces

    */
    virtual void TimeUpdateExternalForces();

    /*!
    \brief treat turbulence models in AssembleMatAndRHS

    */
    virtual void TreatTurbulenceModels(Teuchos::ParameterList& eleparams);

    //! @name Time Step Size Adaptivity
    //@{

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const;

    //! Give local order of accuracy of pressure part
    virtual int MethodOrderOfAccuracyPres() const;

    /*! \brief Return linear error coefficient of velocity
     *
     *  The linear discretization error reads
     *  \f[
     *  e \approx \Delta t_n^2\left(\frac{1}{2}-\theta\right)\ddot{u}(t_n)
     *    + \Delta t_n^3\left(\frac{1}{6}-\frac{\theta}{2}\right)\dddot{u}(t_n)
     *    + HOT\left(\Delta t_n^4\right)
     *  \f]
     *
     *  \author mayr.mt \date 04/2015
     */
    virtual double MethodLinErrCoeffVel() const;

    //@}

   protected:
    /// copy constructor
    TimIntOneStepTheta(const TimIntOneStepTheta& old);

    //! @name time stepping variables
    bool startalgo_;  ///< flag for starting algorithm

    /// the vector containing external loads at t_n
    Teuchos::RCP<Epetra_Vector> external_loadsn_;

    /// the vector containing external loads at t_{n+1}
    Teuchos::RCP<Epetra_Vector> external_loadsnp_;

   private:
  };  // class TimIntOneStepTheta

}  // namespace FLD

#endif  // #ifndef FLD_TIMINT_OST_H
