/*----------------------------------------------------------------------*/
/*!
\file strtimint_ost.H
\brief Structural time integration with one-step-theta
\level 1
\maintainer Alexander Popp
*/

/*----------------------------------------------------------------------*/
#ifndef STRTIMINT_OST_H
#define STRTIMINT_OST_H

/*----------------------------------------------------------------------*/
/* headers */
#include "strtimint_impl.H"
#include "../drt_lib/standardtypes_cpp.H"

/*----------------------------------------------------------------------*/
/* belongs to structural dynamics namespace */
namespace STR
{

  /*====================================================================*/
  /*!
   * \brief One-step-theta time integration (or Crank-Nicolson scheme or generalised trapezoidal rule)
   *
   * <h3> Background </h3>
   * One-step theta time integration is a finite difference method for 1st order ordinary differential equations (ODE) of the type
   * \f[
   *   F(y,\dot{y},t) = \dot{y}(t) - f(y(t),t) = 0
   * \f]
   *
   * The one-step-theta time integration method discretises this equation into the following reference formula
   * \f[
   *   \frac{y_{n+1} - y_n}{\Delta t}
   *   - \theta f(y_{n+1},t_{n+1}) + (1-\theta) f(y_n,t_n)
   *   = 0
   * \qquad
   * \f]
   * in which \f$\theta\in[0,1]\f$ is the key parameter. The method is implicit unless \f$\theta=0\f$, which is the forward Euler scheme. The method recovers the backward Euler method with \f$\theta=1\f$. The trapezoidal rule (TR, or average acceleration method) is obtained with \f$\theta=1/2\f$. Only the trapezoidal rule is second order accurate, all other schemes are only first order.
   *
   * This method is applied to the set of ODEs reflecting the second order degree of the governing equations in structural dynamics:
   * \f[\left\{\begin{array}{rcl}
   *   M \, \dot{D}(t) - M \, V(t) & = & 0
   * \\
   *   M \, \dot{V}(t) + C \, V(t) + F_{int}(D,t) - F_{ext}(t) & = & 0
   * \end{array}\right.\f]
   * in which \f$M\f$ is the global mass matrix (#mass_), \f$C\f$ a global damping matrix (in our algorithm the Rayleigh damping matrix, #damp_), \f$D(t)\f$ the displacements (#dis_), \f$V(t)\f$ the velocities (#vel_), \f$F_{int}\f$ the internal forces and \f$F_{ext}\f$ the external forces. One obtains
   * \f[\left\{\begin{array}{rcl}
   *   \frac{D_{n+1} - D_n}{\Delta t}
   *   - \theta V_{n+1} - (1-\theta) V_n
   *   & = & 0
   * \\
   *   M \frac{V_{n+1} + V_n}{\Delta t}
   *   + C \big( \theta V_{n+1} - (1-\theta) V_n \big)
   *   + F_{int,n+\theta}
   *   - F_{ext,n+\theta}
   *   & = & 0
   * \end{array}\right.\f]
   * with
   * \f[
   *   F_{int,n+\theta}
   *   = \theta F_{int}(D_{n+1},t_{n+1}) + (1-\theta) F_{int}(D_{n+1},t_{n+1})
   *   \quad\mbox{and}\quad
   *   F_{ext,n+\theta}
   *   = \theta F_{ext}(t_{n+1}) + (1-\theta) F_{ext}(t_{n+1})
   * \f]
   * These vector equations can be rewritten such that the unknown velocities \f$V_{n+1}\f$ can be suppressed or rather expressed by the unknown displacements \f$D_{n+1}\f$. The residual is achieved
   * \f[
   *   R_{n+\theta}(D_{n+1})
   *   = M A_{n+\theta}(D_{n+1}) + C V_{n+\theta}(D_{n+1})
   *   + F_{int,n+\theta} - F_{ext,n+\theta}
   *   = 0
   * \f]
   * in which
   * \f[\begin{array}{rclcl}
   *   A_{n+\theta}(D_{n+1})
   *   & = &
   *   \frac{1}{\theta \Delta t^2} ( D_{n+1} - D_n )
   *   - \frac{1}{\theta \Delta t} V_n
   *   & =: &
   *   \theta A_{n+1} + (1-\theta) A_n
   * \\
   *   V_{n+\theta}(D_{n+1})
   *   & = &
   *   \frac{1}{\Delta t} ( D_{n+1} - D_n )
   *   &&
   * \end{array}\f]
   *
   * <h3>Family members to be aware of</h3>
   * <table>
   *   <tr>
   *     <td>Name</td>
   *     <td>Abbreviation</td>
   *     <td>\f$\theta\f$</td>
   *     <td>Order of accuracy</td>
   *     <td>Stability</td>
   *   </tr>
   *   <tr>
   *     <td>Backward Euler</td>
   *     <td>BE</td>
   *     <td>\f$1\f$</td>
   *     <td>1</td>
   *     <td>A,L-stable</td>
   *   </tr>
   *   <tr>
   *     <td>Trapezoidal rule</td>
   *     <td>TR</td>
   *     <td>\f$\frac{1}{2}\f$</td>
   *     <td>2</td>
   *     <td>A-stable</td>
   *   </tr>
   * </table>
   *
   * <h3>References</h3>
   * - [1] HR Schwarz, Numerische Mathematik, Teubner, Stuttgart, 1997.
   * - [2] TJR Hughes, The finite element method, Dover, Mineola, 1987.
   * - [3] P Deuflhard and F Bornemann, Numerische Mathematik II: Integration gewohnlicher Differentialgleichungen, Walter de Gryter, Berlin, 1994.
   * - [4] ...
   *
   *
   * \author bborn
   * \date 06/08
   */
  class TimIntOneStepTheta : public TimIntImpl
  {

  public:

    //! verify if given coefficients are in admissible range;
    //! prints also info to STDOUT
    void VerifyCoeff();

    //! @name Construction
    //@{

    //! Constructor
    TimIntOneStepTheta
    (
      const Teuchos::ParameterList& timeparams,  //!< time parameters
      const Teuchos::ParameterList& ioparams,    //!< ioflags
      const Teuchos::ParameterList& sdynparams,  //!< input parameters
      const Teuchos::ParameterList& xparams,     //!< extra flags
      Teuchos::RCP<DRT::Discretization> actdis,  //!< current discretisation
      Teuchos::RCP<LINALG::Solver> solver,       //!< the solver
      Teuchos::RCP<LINALG::Solver> contactsolver,    //!< the solver for contact meshtying
      Teuchos::RCP<IO::DiscretizationWriter> output  //!< the output
    );

    //! Destructor
    // ....

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return bool
    \date 08/16
    \author rauch  */
    virtual void Init(
        const Teuchos::ParameterList& timeparams,
        const Teuchos::ParameterList& sdynparams,
        const Teuchos::ParameterList& xparams,
        Teuchos::RCP<DRT::Discretization> actdis,
        Teuchos::RCP<LINALG::Solver> solver );

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

    //! Resize #TimIntMStep<T> multi-step quantities
    //! Single-step method: nothing to do here
    void ResizeMStep() { ; }

    //@}

    //! @name Pure virtual methods which have to be implemented
    //@{

    //! Return name
    enum INPAR::STR::DynamicType MethodName() const
    {
      return INPAR::STR::dyna_onesteptheta;
    }

    //! Provide number of steps, a single-step method returns 1
    virtual int MethodSteps() const { return 1; }

    //! Give local order of accuracy of displacement part
    virtual int MethodOrderOfAccuracyDis() const
    {
      return fabs(MethodLinErrCoeff1()) < EPS6 ? 2 : 1;
    }

    //! Give local order of accuracy of velocity part
    virtual int MethodOrderOfAccuracyVel() const
    {
      return MethodOrderOfAccuracyDis();
    }

    //! Return linear error coefficient of displacements
    virtual double MethodLinErrCoeffDis() const
    {
      if (MethodOrderOfAccuracyDis() == 1)
        return MethodLinErrCoeff1();
      else
        return MethodLinErrCoeff2();
    }

    //! Return linear error coefficient of velocities
    virtual double MethodLinErrCoeffVel() const
    {
      return MethodLinErrCoeffDis();
    }

    //! Linear error coefficient if 1st order accurate
    virtual double MethodLinErrCoeff1() const
    {
      return 1./2. - theta_;
    }

    //! Linear error coefficient if 2nd order accurate
    virtual double MethodLinErrCoeff2() const
    {
      return 1./6. - theta_/2.;  // this is -1/12
    }

    //! return time integration factor
    virtual double TimIntParam() const {return 1.0 - theta_; }

    //! Consistent predictor with constant displacements
    //! and consistent velocities and displacements
    void PredictConstDisConsistVelAcc();

    //! Consistent predictor with constant velocities,
    //! extrapolated displacements and consistent accelerations
    //! \author mayr.mt
    void PredictConstVelConsistAcc();

    //! Consistent predictor with constant accelerations
    //! and extrapolated velocities and displacements
    //! \author mayr
    void PredictConstAcc();

    //! Create force residual #fres_ and its stiffness #stiff_
    void EvaluateForceStiffResidual(Teuchos::ParameterList& params);

    //! Evaluate/define the residual force vector #fres_ for
    //! relaxation solution with SolveRelaxationLinear
    void EvaluateForceStiffResidualRelax(Teuchos::ParameterList& params);

    //! Evaluate residual #fres_
    void EvaluateForceResidual();

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
    void UpdateStepState();

    //! Update element
    void UpdateStepElement();

    //! Read and set restart for forces
    void ReadRestartForce();

    //! Write internal and external forces for restart
    void WriteRestartForce(Teuchos::RCP<IO::DiscretizationWriter> output);

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

    //! @name One-Step-Theta specific methods
    //@{

    //! Evaluate mid-state vectors by averaging end-point vectors
    void EvaluateMidState();

    //@}

  protected:

    //! equal operator is NOT wanted
    TimIntOneStepTheta operator=(const TimIntOneStepTheta& old);

    //! copy constructor is NOT wanted
    TimIntOneStepTheta(const TimIntOneStepTheta& old);

    //! Clear mass matrix and evaluate mass matrix again.
    //! \note not implemented in base class.
    virtual void DetermineMass();

    //! @name Key coefficients
    //@{
    double theta_;  //!< factor (0,1]
    //@}

    //! @name Global mid-state vectors
    //@{
    //! mid-displacements \f$D_m = D_{n+\theta}\f$
    Teuchos::RCP<Epetra_Vector> dist_;
    //! mid-velocities \f$V_m = V_{n+1-\theta}\f$
    Teuchos::RCP<Epetra_Vector> velt_;
    //! mid-accelerations \f$A_m = A_{n+1-\theta}\f$
    Teuchos::RCP<Epetra_Vector> acct_;
    //@}

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    Teuchos::RCP<Epetra_Vector> fint_;  //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;  //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> finert_;  //!< inertia force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> finertt_;  //!< inertia mid-force
    Teuchos::RCP<Epetra_Vector> finertn_;  //!< inertia force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fvisct_;  //!< viscous force

    //@}

  };  // class TimIntOneStepTheta

}  // namespace STR

/*----------------------------------------------------------------------*/
#endif  // #ifndef STRTIMINT_OST_H
