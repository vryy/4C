/*----------------------------------------------------------------------*/
/*!
\file strtimint_genalpha.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRTIMINT_GENALPHA_H
#define STRTIMINT_GENALPHA_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_impl.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Generalised-alpha time integration
   *
   * References
   * - [1] NM Newmark, A method of computation for structural dynamics,
   *   Journal of Earthquake Mechanic Division, 85(EM3):67-94, 1959.
   * - [2] HM Hilber, TJR Hughes and RL Taylor, Improved numerical
   *   dissipation for time integration algorithms, Earthquake 
   *   Engineering and Structural Mechanics, 5:283-292, 1977.
   *     
   *
   * \author bborn
   * \date 06/08
   */
  class TimIntGenAlpha : public TimIntImpl
  {

    //! Zienkiewicz-Xie auxiliar scheme is friend
    friend class TimAdaZienXie;

  public:

    //! Mid-average type of internal forces
    enum MidAverageEnum
    {
      midavg_vague = 0,  //!< undefined mid-averaging type
      midavg_imrlike,  //!< alphaf-mid-averaging is done IMR-like, i.e. 
                       //!< \f$F_{int,m}\f$ 
                       //!< \f$= F_{int}(D_m)\f$
                       //!< \f$= F_{int}( (1-\alpha_f)*D_{n+1} + \alpha_f*D_n )\f$
                       //!< (IMR means implicit mid-point rule.)
      midavg_trlike    //!< alphaf-mid-averaging is done TR-like, i.e.
                       //!< \f$F_{int,m}\f$
                       //!< \f$= (1-\alpha_f)*F_{int,n+1} + \alpha_f*F_{int,n}\f$
                       //!<  \f$= (1-\alpha_f)*F_{int}(D_{n+1}) + \alpha_f*F_{int}(D_n)\f$
                       //!<  (TR means trapezoidal rule.)
    };

    //! convert input string to enum for mid-average type
    enum MidAverageEnum MapMidAvgStringToEnum
    (
      const std::string name  //!< input string
    );

    //! @name Construction
    //@{

    //! Constructor
    TimIntGenAlpha
    (
      const Teuchos::ParameterList& ioparams,  //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,  //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,  //!< the solver
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    // ....

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here
    void ResizeMStep() { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum TimInt::NameEnum MethodName()
    { 
      return TimInt::name_genalpha;
    }

    //! Provide number of steps, e.g. a single-step method returns 1,
    //! a m-multistep method returns m
    int MethodSteps()
    { 
      return 1;
    }

    //! Give linear order of accuracy of displacement part
    int MethodOrderOfAccuracyDis()
    { 
      return (fabs(MethodLinErrCoeffDis2()) < EPS6) ? 3 : 2;
    }

    //! Give linear order of accuracy of velocity part
    int MethodOrderOfAccuracyVel()
    { 
      return (fabs(MethodLinErrCoeffVel1()) < EPS6) ? 2 : 1;
    }

    //! Return linear error coefficient of displacements
    double MethodLinErrCoeffDis()
    {
      if (MethodOrderOfAccuracyDis() == 2)
        return MethodLinErrCoeffDis2();
      else 
        return MethodLinErrCoeffDis3();
    }

    //! 2nd order linear error coefficient of displacements
    double MethodLinErrCoeffDis2()
    {
      return 1./6. - beta_*(1.-alphaf_)/(1.-alpham_);
    }

    //! 3rd order linear error coefficient of displacements
    double MethodLinErrCoeffDis3()
    {
      return 1./12. - beta_*(1.-alphaf_)/(1.-alpham_);
    }

    //! Return linear error coefficient of velocities
    double MethodLinErrCoeffVel()
    {
      if (MethodOrderOfAccuracyVel() == 1)
        return MethodLinErrCoeffVel1();
      else 
        return MethodLinErrCoeffVel2();
    }

    //! 1st order linear error coefficient of velocities
    double MethodLinErrCoeffVel1()
    {
      return 1./2. - gamma_*(1.-alphaf_)/(1.-alpham_);
    }

    //! 2nd order linear error coefficient of velocities
    double MethodLinErrCoeffVel2()
    {
      return 1./3. - gamma_*(1.-alphaf_)/(1.-alpham_);
    }

    //! Consistent predictor with constant displacements
    //! and consistent velocities and displacements
    void PredictConstDisConsistVelAcc();

    //! Create force residual #fres_ and ist stiffness #stiff_
    void EvaluateForceStiffResidual();

    //! Determine characteristic norm for displacements
    //! \author lw (originally)
    double CalcRefNormDisplacement();

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce();

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #veln_ and #acc_
    //! from scratch by using the newly updated #disn_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors 
    //! (like #PredictConstDisConsistVelAcc)
    void UpdateIterIncrementally();

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #disn_, #veln_ and #accn_ by
    //! incrementing these vector proportional to the residual
    //! displacements #disi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual displacements #disi_ are blanked at these DOFs.
    void UpdateIterIteratively();

    //! Update step
    void UpdateStep();

    //@}

    //! @name System vectors
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext()
    {
      return fext_;
    }

    //@}

    //! @name Generalised-alpha specific methods
    //@{
    //! Evaluate mid-state vectors by averaging end-point vectors
    void EvaluateMidState();
    //@}

  protected:

    //! equal operator is NOT wanted
    TimIntGenAlpha operator=(const TimIntGenAlpha& old);

    //! copy constructor is NOT wanted
    TimIntGenAlpha(const TimIntGenAlpha& old);

    //! @name set-up
    //@{
    //! mid-average type more at #MidAverageEnum
    enum MidAverageEnum midavg_;
    //@}

    //! @name Key coefficients
    //! Please note, to obtain a second-order accurate scheme, you need
    //! to follow the following formulas in which \f$\rho_\infty\f$ is the
    //! spectral radius.
    //! \f[ \alpha_m = (2*\rho_\infty - 1)/(\rho_\infty + 1) \f]
    //! \f[ \alpha_f = \rho_\infty/(\rho_\infty + 1) \f]
    //! \f[ \beta = 1/4*(1 - \alpha_m + \alpha_f)^2 \mbox{(max. damp. of high-freq. modes)} \f]
    //! \f[ \gamma = 1/2 - \alpha_m + \alpha_f \f]
    //! The spectral radius is responsible for the magnitude of 
    //! numerical dissipation introduced. 
    //! For instance
    //! Without numerical dissipation at \f$\rho_\infty=1\f$
    //! \f[ \beta=0.25, \gamma=0.5, \alpha_m=0.5, \alpha_f=0.5 \f]
    //! Medium dissipation at \f$\rho_\infty=0.8\f$
    //! \f[ \beta=25/81, \gamma=11/18, \alpha_m=1/3, \alpha_f=4/9  \f]
    //! Strong numerical dissipation at \f$\rho_\infty=0.5\f$
    //! \f[ \beta=4/9, \gamma=10/12, \alpha_m=0, \alpha_f=1/3 \f]
    //@{
    double beta_;  //!< factor (0,1/2]
    double gamma_;  //!< factor (0,1]
    double alphaf_;  //!< factor [0,1)
    double alpham_;  //!< factor [0,1)
    //@}

    //! @name Global mid-state vectors
    //@{
    Teuchos::RCP<Epetra_Vector> dism_;  //!< mid-displacements
                                        //!< \f$D_m = D_{n+1-\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> velm_;  //!< mid-velocities
                                        //!< \f$V_m = V_{n+1-\alpha_f}\f$
    Teuchos::RCP<Epetra_Vector> accm_;  //!< mid-accelerations
                                        //!< \f$A_m = A_{n+1-\alpha_m}\f$
    //@}

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    Teuchos::RCP<Epetra_Vector> fint_;  //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintm_;  //!< internal mid-force
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;  //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextm_;  //!< external mid-force
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> finertm_;  //!< inertial force
    Teuchos::RCP<Epetra_Vector> fviscm_;  //!< viscous force

    Teuchos::RCP<Epetra_Vector> frobin_;  //!< Robin force due to FSI
    //@}

  };  // class TimIntGenAlpha

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_GENALPHA_H
#endif  // #ifdef CCADISCRET
