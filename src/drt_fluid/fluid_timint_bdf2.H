/*-----------------------------------------------------------*/
/*!

\brief BDF-2 time-integration scheme

\maintainer Martin Kronbichler

\level 2

*/
/*-----------------------------------------------------------*/


#ifndef FLD_TIMINT_BDF2_H
#define FLD_TIMINT_BDF2_H


#include "fluidimplicitintegration.H"


namespace FLD
{
  class TimIntBDF2 : public virtual FluidImplicitTimeInt
  {
   public:
    /// Standard Constructor
    TimIntBDF2(const Teuchos::RCP<DRT::Discretization>& actdis,
        const Teuchos::RCP<LINALG::Solver>& solver,
        const Teuchos::RCP<Teuchos::ParameterList>& params,
        const Teuchos::RCP<IO::DiscretizationWriter>& output, bool alefluid = false);

    /// Destructor
    virtual ~TimIntBDF2();

    /*!
    \brief initialization

    */
    virtual void Init();

    /*!
    \brief Print information about current time step to screen

    */
    virtual void PrintTimeStepInfo();

    /*!
    \brief Set theta_ to its value, dependent on integration method for GenAlpha and BDF2

    */
    virtual void SetTheta();

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
    \brief return scheme-specific time integration parameter

    */
    virtual double TimIntParam() const { return 0.0; }

    /*!
    \brief return scaling factor for the residual

    */
    virtual double ResidualScaling() const { return 1.0 / (theta_ * dta_); }

    //! @name Time Step Size Adaptivity
    //@{

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const { return 2; };

    //! Give local order of accuracy of pressure part
    virtual int MethodOrderOfAccuracyPres() const { return 2; };

    /*! \brief Return linear error coefficient of velocity
     *
     *  The linear discretization error reads
     *  \f[
     *  e \approx -\frac{\left(\Delta t_n + \Delta t_{n-1}\right)^2}
     *                  {6\Delta t_n\left(2\Delta t_n+\Delta t_{n-1}\right)}
     *             \Delta t_n^3 \dddot{u}(t_n)
     *    + HOT\left(\Delta t_n^4\right)
     *  \f]
     *
     *  \author mayr.mt \date 04/2015
     */
    virtual double MethodLinErrCoeffVel() const;

    //@}

   protected:
    /// copy constructor
    TimIntBDF2(const TimIntBDF2& old);

    /*!
    \brief velocity required for evaluation of related quantities required on element level

    */
    virtual Teuchos::RCP<const Epetra_Vector> EvaluationVel() { return velnp_; };


   private:
    /// time factor for one-step-theta/BDF2 time integration
    double theta_;

  };  // class TimIntBDF2

}  // namespace FLD


#endif  // #ifndef FLD_TIMINT_BDF2_H
