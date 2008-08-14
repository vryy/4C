/*----------------------------------------------------------------------*/
/*!
\file strtimint_gemm.H

<pre>
Maintainer: Burkhard Bornemann
            bornemann@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15237
</pre>
*/

/*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef STRTIMINT_GEMM_H
#define STRTIMINT_GEMM_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_impl.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief Generalised energy-momentum method (GEMM)
   *
   * <b>References</b>
   * - [1] D Kuhl and MA Crisfield, Energy-conserving and decaying algorithms
   *   in non-linear structural dynamics, International Journal for
   *   Numerical Methods in Engineering, 45(5):569-599, 1999.
   * - [2] JC Simo and N Tarnow, The discrete energy-momentum method, 
   *   Conserving algorithms fro nonlinear elastodynamics, Journal of
   *   Applied Mathematics and Physics, 43:757-792, 1992.
   *
   * \author bborn
   * \date 07/08
   */
  class TimIntGEMM : public TimIntImpl
  {

    //! Zienkiewicz-Xie auxiliar scheme is friend
    friend class TimAdaZienXie;

  public:

    //! @name Lifetime
    //@{

    //! Constructor
    TimIntGEMM
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
      return TimInt::name_gemm;
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
      // CHECK THIS
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

    //! @name Access methods
    //@{

    //! Return external force \f$F_{ext,n}\f$
    Teuchos::RCP<Epetra_Vector> Fext()
    {
      return fext_;
    }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew()
    {
      return fextn_;
    }

    //@}

    //! @name GEMM-specific methods
    //@{

    //! Evaluate mid-state vectors by averaging end-point vectors
    void EvaluateMidState();

    //! Evaluate ordinary internal force, its stiffness at mid-state
    void ApplyForceStiffInternalMid
    (
      const double time,  //!< evaluation time
      const double dt,  //!< step size
      const Teuchos::RCP<Epetra_Vector> dis,  //!< displacement state at t_n
      const Teuchos::RCP<Epetra_Vector> disn,  //!< displacement state at t_{n+1}
      const Teuchos::RCP<Epetra_Vector> disi,  //!< residual displacements
      const Teuchos::RCP<Epetra_Vector> vel,  // velocity state
      Teuchos::RCP<Epetra_Vector> fint,  //!< internal force
      Teuchos::RCP<LINALG::SparseMatrix> stiff  //!< stiffness matrix
    );
    //@}

  protected:

    //! equal operator is NOT wanted
    TimIntGEMM operator=(const TimIntGEMM& old);

    //! copy constructor is NOT wanted
    TimIntGEMM(const TimIntGEMM& old);

    //! @name Key coefficients
    //@{
    const double beta_;  //!< must be 0.25 
    const double gamma_;  //!< must be 0.5
    double alphaf_;  //!< factor [0,1)
    double alpham_;  //!< factor [0,1)
    double xi_;  //!< factor [0,1)
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
    //Teuchos::RCP<Epetra_Vector> fint_;  //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintm_;  //!< internal mid-force
    //Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;  //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextm_;  //!< external mid-force
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> finertm_;  //!< inertial force
    Teuchos::RCP<Epetra_Vector> fviscm_;  //!< viscous force

    Teuchos::RCP<Epetra_Vector> frobin_;  //!< Robin force due to FSI
    //@}

  };  // class TimIntGEMM

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_GEMM_H
#endif  // #ifdef CCADISCRET
