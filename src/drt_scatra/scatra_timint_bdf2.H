/*!----------------------------------------------------------------------
\file scatra_timint_bdf2.H
\brief BDF2 time-integration scheme

<pre>
\level 1

\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_BDF2_H
#define SCATRA_TIMINT_BDF2_H

#include "scatra_timint_implicit.H"

namespace SCATRA
{
class TimIntBDF2 : public virtual ScaTraTimIntImpl
{
public:
  /// Standard Constructor
  TimIntBDF2(Teuchos::RCP<DRT::Discretization> dis,
      Teuchos::RCP<LINALG::Solver>            solver,
      Teuchos::RCP<Teuchos::ParameterList>    params,
      Teuchos::RCP<Teuchos::ParameterList>    extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntBDF2();

  /// Setup time integration scheme
  virtual void Setup();

  /// compute values at intermediate time steps (required for generalized-alpha)
  void ComputeIntermediateValues(){return;};

  /// compute values at the interior of the elements (required for hdg)
  virtual void ComputeInteriorValues(){return;};

  ///  compute scalar time derivative
  void ComputeTimeDerivative();

  ///  compute scalar time derivate parameters of the input voltage to compute double layer current densities
  virtual void ComputeTimeDerivPot0(const bool init){return;};

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  virtual void Update(const int num=0);

  /// read restart data
  virtual void ReadRestart(const int step,Teuchos::RCP<IO::InputControl> input = Teuchos::null);

  // routine to return scalar field phi at time step n-1
  Teuchos::RCP<Epetra_Vector> Phinm(){return Teuchos::null;}

  /// routine to return scalar field phi at time step n+alpha_F
  Teuchos::RCP<Epetra_Vector> Phiaf(){return Teuchos::null;}

  /// routine to return scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phiam(){return Teuchos::null;}

  /// routine to return time derivative of scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phidtam(){return Teuchos::null;}

  /// routine to return fine-scale scalar field fsphi at time step n+1
  virtual Teuchos::RCP<Epetra_Vector> FsPhi() {
    if (Sep_ != Teuchos::null)
    Sep_->Multiply(false,*phinp_,*fsphinp_);
    return fsphinp_; };

  /// routine to return time integration specific parameters
  Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList()
  {
    Teuchos::RCP<Teuchos::ParameterList> timeparams;
    timeparams = Teuchos::rcp(new Teuchos::ParameterList());
    timeparams->set("using stationary formulation",false);
    timeparams->set("using generalized-alpha time integration",false);
    timeparams->set("total time",time_);
    timeparams->set("time factor",theta_*dta_);
    timeparams->set("alpha_F",1.0);
    return timeparams;
  }


protected:

  /// don't want = operator and cctor
  TimIntBDF2 operator = (const TimIntBDF2& old);

  /// copy constructor
  TimIntBDF2 (const TimIntBDF2& old);

  /// set time parameter for element evaluation
  virtual void SetElementTimeParameter(bool forcedincrementalsolver=false) const;

  //! set time for evaluation of Neumann boundary conditions
  virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

  //! calculate consistent initial conditions in compliance with initial scalar field
  //! this is not necessary for BDF2 time integration scheme
  void CalcInitialTimeDerivative(){return;};

  /// Set the part of the righthandside belonging to the last timestep.
  virtual void SetOldPartOfRighthandside();

  /// do explicit predictor step (-> better starting value for nonlinear solver)
  void ExplicitPredictor();

  /// add actual Neumann loads with time factor
  void AddNeumannToResidual();

  /// AVM3-based scale separation
  void AVM3Separation();

  /// dynamic Smagorinsky model
  void DynamicComputationOfCs();

  /// dynamic Vreman model
  void DynamicComputationOfCv();

  /// add parameters specific for time-integration scheme
  void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver=false);

  /// write additional data required for restart
  virtual void OutputRestart();

  /// return the right time-scaling-factor for the true residual
  double ResidualScaling() const { return 1.0/(dta_*theta_); }

  /// time factor for one-step-theta/BDF2 time integration
  double theta_;

  /// solution vector phi at time n-1
  Teuchos::RCP<Epetra_Vector>  phinm_;

  /// fine-scale solution vector at time n+1
  Teuchos::RCP<Epetra_Vector>  fsphinp_;
}; // class TimIntBDF2
} // namespace SCATRA
#endif  // #ifndef SCATRA_TIMINT_BDF2_H
