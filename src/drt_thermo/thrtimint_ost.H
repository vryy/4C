/*----------------------------------------------------------------------*/
/*! \file
\brief one-step-theta time integrator for thermal field
\level 1
\maintainer Christoph Meier

*/


/*----------------------------------------------------------------------*
 | definitions                                               dano 08/09 |
 *----------------------------------------------------------------------*/
#ifndef THRTIMINT_OST_H
#define THRTIMINT_OST_H


/*----------------------------------------------------------------------*
 | headers                                                   dano 08/09 |
 *----------------------------------------------------------------------*/
#include "thrtimint_impl.H"


/*----------------------------------------------------------------------*
 | belongs to thermal dynamics namespace                     dano 08/09 |
 *----------------------------------------------------------------------*/
namespace THR
{
  /*====================================================================*/
  /*!
   * \brief One-step-theta time integration (or Crank-Nicholson scheme or generalised trapezoidal
   * rule)
   *
   * <h3> Background </h3>
   * One-step theta time integration is a finite difference method for 1st order ordinary
   * differential equations (ODE) of the type \f[ F(y,\dot{y},t) = \dot{y}(t) - f(y(t),t) = 0 \f]
   *
   * The one-step-theta time integration method discretises this equation into the following
   * reference formula \f[ \frac{y_{n+1} - y_n}{\Delta t}
   *   - \theta f(y_{n+1},t_{n+1}) + (1-\theta) f(y_n,t_n)
   *   = 0
   * \qquad
   * \f]
   * in which \f$\theta\in[0,1]\f$ is the key parameter. The method is implicit unless
   * \f$\theta=0\f$, which is the forward Euler scheme. The method recovers the backward Euler
   * method with \f$\theta=1\f$. The trapezoidal rule (TR, or average acceleration method) is
   * obtained with \f$\theta=1/2\f$. Only the trapezoidal rule is second order accurate, all other
   * schemes are only first order.
   *
   * This method is applied to the set of ODEs reflecting the first order degree of the governing
   * equations in thermal dynamics: \f[\left\{\begin{array}{rcl}
   *
   *   C \, R(t) + F_{int}(T,t) - F_{ext}(t) & = & 0
   * \end{array}\right.\f]
   * \f$C\f$ is a global capacity matrix, \f$T(t)\f$ the temperature (#temp_), \f$R(t)\f$ the
   * temperature rates (#rate_), \f$F_{int}\f$ the internal forces and \f$F_{ext}\f$ the external
   * forces. One obtains \f[\left\{\begin{array}{rcl} \frac{T_{n+1} - T_n}{\Delta t}
   *   - \theta R_{n+1} - (1-\theta) R_n
   *   & = & 0
   * \\
   *   C \frac{T_{n+1} - T_n}{\Delta t}
   *   + F_{int,n+\theta}
   *   - F_{ext,n+\theta}
   *   & = & 0
   * \end{array}\right.\f]
   * with
   * \f[
   *   F_{int,n+\theta}
   *   = \theta F_{int}(T_{n+1},t_{n+1}) + (1-\theta) F_{int}(T_{n+1},t_{n+1})
   *   \quad\mbox{and}\quad
   *   F_{ext,n+\theta}
   *   = \theta F_{ext}(t_{n+1}) + (1-\theta) F_{ext}(t_{n+1})
   * \f]
   * These vector equations can be rewritten such that the unknown temperature rates \f$r_{n+1}\f$
   * can be suppressed or rather expressed by the unknown temperatures \f$T_{n+1}\f$. The residual
   * is achieved \f[ R_{n+\theta}(T_{n+1}) = C R_{n+\theta}(T_{n+1})
   *   + F_{int,n+\theta} - F_{ext,n+\theta}
   *   = 0
   * \f]
   * in which
   * \f[\begin{array}{rclcl}
   *   R_{n+\theta}(T_{n+1})
   *   & = &
   *   \frac{1}{\Delta t} ( T_{n+1} - T_n )
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
   * - [3] P Deuflhard and F Bornemann, Numerische Mathematik II: Integration gewohnlicher
   * Differentialgleichungen, Walter de Gryter, Berlin, 1994.
   * - [4] ...
   *
   *
   * \author bborn
   * \date 06/08
   */
  class TimIntOneStepTheta : public TimIntImpl
  {
   public:
    //! verify if given coefficients are in admissable range
    //! prints also info to STDOUT
    void VerifyCoeff();

    //! @name Construction
    //@{

    //! Constructor
    TimIntOneStepTheta(const Teuchos::ParameterList& ioparams,  //!< ioflags
        const Teuchos::ParameterList& tdynparams,               //!< input parameters
        const Teuchos::ParameterList& xparams,                  //!< extra flags
        Teuchos::RCP<DRT::Discretization> actdis,               //!< current discretisation
        Teuchos::RCP<LINALG::Solver> solver,                    //!< the solver
        Teuchos::RCP<IO::DiscretizationWriter> output           //!< the output
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
    enum INPAR::THR::DynamicType MethodName() const { return INPAR::THR::dyna_onesteptheta; }

    //! Provide number of steps, a single-step method returns 1
    int MethodSteps() { return 1; }

    //! Give local order of accuracy of temperature part
    int MethodOrderOfAccuracy() { return fabs(1. / 2. - theta_) < EPS10 ? 2 : 1; }

    //! Return linear error coefficient
    virtual double MethodLinErrCoeff() { return 1. / 2. - theta_; }

    //! Consistent predictor with constant temperatures
    //! and consistent temperature rates and temperatures
    void PredictConstTempConsistRate();

    //! Evaluate ordinary internal force, its tangent at state
    void ApplyForceTangInternal(const double time,  //!< evaluation time
        const double dt,                            //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,     //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,    //!< residual temperatures
        Teuchos::RCP<Epetra_Vector> fcap,           //!< capacity force
        Teuchos::RCP<Epetra_Vector> fint,           //!< internal force
        Teuchos::RCP<LINALG::SparseMatrix> tang     //!< tangent matrix
    );

    //! Evaluate ordinary internal force
    void ApplyForceInternal(const double time,    //!< evaluation time
        const double dt,                          //!< step size
        const Teuchos::RCP<Epetra_Vector> temp,   //!< temperature state
        const Teuchos::RCP<Epetra_Vector> tempi,  //!< incremental temperatures
        Teuchos::RCP<Epetra_Vector> fint          //!< internal force
    );

    //! Evaluate a convective boundary condition
    // (nonlinear --> add term to tangent)
    void ApplyForceExternalConv(const double time,  //!< evaluation time
        const Teuchos::RCP<Epetra_Vector> tempn,    //!< temperature state T_n
        const Teuchos::RCP<Epetra_Vector> temp,     //!< temperature state T_n+1
        Teuchos::RCP<Epetra_Vector> fext,           //!< internal force
        Teuchos::RCP<LINALG::SparseMatrix> tang     //!< tangent matrix
    );

    //! Create force residual #fres_ and its tangent #tang_
    void EvaluateRhsTangResidual();

    //! Determine characteristic norm for temperatures
    //! \author lw (originally)
    double CalcRefNormTemperature();

    //! Determine characteristic norm for force
    //! \author lw (originally)
    double CalcRefNormForce();

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #tempn_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #PredictConstTempConsistRate)
    void UpdateIterIncrementally();

    //! Update iteration iteratively
    //!
    //! This is the ordinary update of #tempn_ and #raten_ by
    //! incrementing these vector proportional to the residual
    //! temperatures #tempi_
    //! The Dirichlet BCs are automatically respected, because the
    //! residual temperatures #tempi_ are blanked at these DOFs.
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
    Teuchos::RCP<Epetra_Vector> Fext() { return fext_; }

    //! Return external force \f$F_{ext,n+1}\f$
    Teuchos::RCP<Epetra_Vector> FextNew() { return fextn_; }

    //@}

    //! @name Generalised-alpha specific methods
    //@{

    //! Evaluate mid-state vectors by averaging end-point vectors
    void EvaluateMidState();

    //@}

   protected:
    //! equal operator is NOT wanted
    TimIntOneStepTheta operator=(const TimIntOneStepTheta& old);

    //! copy constructor is NOT wanted
    TimIntOneStepTheta(const TimIntOneStepTheta& old);

    //! @name Key coefficients
    //@{
    double theta_;  //!< factor (0,1]
    //@}

    //! @name Global mid-state vectors
    //@{
    //! mid-temperatures \f$T_m = T_{n+\theta}\f$
    Teuchos::RCP<Epetra_Vector> tempt_;

    //@}

    //! @name Global force vectors
    //! Residual \c fres_ exists already in base class
    //@{
    Teuchos::RCP<Epetra_Vector> fint_;   //!< internal force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fintn_;  //!< internal force at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fcap_;  //!< capacity force \f$C\cdot\Theta_n\f$ at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector>
        fcapn_;  //!< capacity force \f$C\cdot\Theta_{n+1}\f$ at \f$t_{n+1}\f$

    Teuchos::RCP<Epetra_Vector> fext_;   //!< external force at \f$t_n\f$
    Teuchos::RCP<Epetra_Vector> fextn_;  //!< external force at \f$t_{n+1}\f$

    //@}

  };  // class TimIntOneStepTheta

}  // namespace THR

/*----------------------------------------------------------------------*/
#endif  // #ifndef THRTIMINT_OST_H
