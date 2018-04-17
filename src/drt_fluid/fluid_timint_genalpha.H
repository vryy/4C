/*!----------------------------------------------------------------------
\file fluid_timint_genalpha.H

\brief Generalized-alpha time-integration scheme

\level 2
<pre>
\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifndef FLD_TIMINT_GENALPHA_H
#define FLD_TIMINT_GENALPHA_H


#include "fluidimplicitintegration.H"


namespace FLD
{

/*! \brief Generalized-Alpha time integration for fluid problems
 *
 *  <h3> References </h3>
 *  <ul>
 *  <li>[1] Jansen KE, Whiting CH, Hulbert GM: A generalized-\f$\alpha\f$ method
 *          for integrating the filered Navier--Stokes euqations with a
 *          stabilized finite element method, Comp. Meth. Appl. Mech. Engng.,
 *          190(3-4), pp. 305--319, 2000 </li>
 *  </ul>
 */
class TimIntGenAlpha : public virtual FluidImplicitTimeInt
{

public:

  /// Standard Constructor
  TimIntGenAlpha(      const Teuchos::RCP<DRT::Discretization>&      actdis,
      const Teuchos::RCP<LINALG::Solver>&           solver,
      const Teuchos::RCP<Teuchos::ParameterList>&   params,
      const Teuchos::RCP<IO::DiscretizationWriter>& output,
      bool                                          alefluid = false);

  /// Destructor
  virtual ~TimIntGenAlpha();

  /*!
  \brief initialization

  */
  virtual void Init();

  /*!
  \brief Print information about current time step to screen (reimplementation)

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
  \brief Set time factor in GenAlpha

  */
  virtual double SetTimeFac();

  /*!
  \brief Calculate time derivatives for
         stationary/one-step-theta/BDF2/af-generalized-alpha time integration
         for incompressible and low-Mach-number flow
  */
  virtual void CalculateAcceleration(
      const Teuchos::RCP<const Epetra_Vector>    velnp,    ///< velocity at     n+1
      const Teuchos::RCP<const Epetra_Vector>    veln,     ///< velocity at     n
      const Teuchos::RCP<const Epetra_Vector>    velnm,    ///< velocity at     n-1
      const Teuchos::RCP<const Epetra_Vector>    accn ,    ///< acceleration at n-1
      const Teuchos::RCP<Epetra_Vector>          accnp     ///< acceleration at n+1
  );

  /*!
  \brief compute values at intermediate time steps for gen.-alpha
         for given vectors and store result in given vectors
         Helper method which can be called from outside fluid (e.g. for coupled problems)

  */
  virtual void GenAlphaIntermediateValues(
      Teuchos::RCP<Epetra_Vector>& vecnp,
      Teuchos::RCP<Epetra_Vector>& vecn
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
  \brief Update velaf_ for GenAlpha

  */
  virtual void UpdateVelafGenAlpha();

  /*!
  \brief Output of filtered velocity

  */
  virtual void OutputofFilteredVel(Teuchos::RCP<Epetra_Vector> outvec,
                           Teuchos::RCP<Epetra_Vector> fsoutvec);

  /*!

  \brief parameter (fix over a time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  virtual void SetElementTimeParameter();

  /*!
  \brief return scheme-specific time integration parameter

  */
  virtual double TimIntParam() const { return (1.0 - alphaF_); }

  virtual double AlphaM() const {return alphaM_;}

  //! @name methods for fsi
  /// Extrapolation of vectors from mid-point to end-point t_{n+1}
  Teuchos::RCP<Epetra_Vector> ExtrapolateEndPoint(
      Teuchos::RCP<Epetra_Vector> vecn, ///< vector at time level t_n
      Teuchos::RCP<Epetra_Vector> vecm  ///< vector at time level of equilibrium
      );

  /*!
  \brief treat turbulence models in AssembleMatAndRHS

  */
  virtual void TreatTurbulenceModels(Teuchos::ParameterList& eleparams);

  //! @name Time Step Size Adaptivity
  //@{

  /*! \brief Give local order of accuracy of velocity part
   *
   *  Generalized-Alpha is second order accurate if
   *  \f$\gamma=\frac{1}{2}+\alpha_{m}-\alpha_{f}\f$
   *  (cf. eq. (17) in [Jansen, 2000]), otherwise only first order.
   *
   *  \author mayr.mt \date 04/2015
   */
  virtual int MethodOrderOfAccuracyVel() const;

  /*! \brief Give local order of accuracy of pressure part
   *
   *  Generalized-Alpha is second order accurate if
   *  \f$\gamma=\frac{1}{2}+\alpha_{m}-\alpha_{f}\f$
   *  (cf. eq. (17) in [Jansen, 2000]), otherwise only first order.
   *
   *  \author mayr.mt \date 04/2015
   */
  virtual int MethodOrderOfAccuracyPres() const;

  /*! \brief Return linear error coefficient of velocity
   *
   *  The linear discretization error reads
   *  \f[
   *  e \approx \Delta t_n^2\left(\frac{1}{2}-\gamma\right)\ddot{u}(t_n)
   *    + \Delta t_n^3\left(\frac{1}{6}-\frac{\gamma}{2}\right)\dddot{u}(t_n)
   *    + HOT\left(\Delta t_n^4\right)
   *  \f]
   *
   *  \author mayr.mt \date 04/2015
   */
  virtual double MethodLinErrCoeffVel() const;

  //@}

protected:

  /// copy constructor
  TimIntGenAlpha (const TimIntGenAlpha& old);

  /*!
  \brief update acceleration for generalized-alpha time integration

  */
  virtual void GenAlphaUpdateAcceleration();

  /*!
  \brief compute values at intermediate time steps for gen.-alpha

  */
  virtual void GenAlphaIntermediateValues() ;

  /*!
  \brief return scaling of the residual

  */
  virtual double ResidualScaling() const
  {

    return alphaM_/(gamma_*dta_);

  }

  /*!
  \brief velocity required for evaluation of related quantites required on element level

  */
  virtual Teuchos::RCP<const Epetra_Vector> EvaluationVel() {return velaf_;};

  /// time factors for generalized-alpha time integration
  double alphaM_;
  double alphaF_;
  double gamma_;

  //! @name time stepping variables
  bool    startalgo_;   ///< flag for starting algorithm


}; // class TimIntGenAlpha

} // namespace FLD

#endif  // #ifndef FLD_TIMINT_GENALPHA_H
