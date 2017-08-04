/*!----------------------------------------------------------------------
\file scatra_timint_genalpha.H
\brief Generalized-alpha time-integration scheme

<pre>
\level 1

\maintainer Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_GENALPHA_H
#define SCATRA_TIMINT_GENALPHA_H

#include "scatra_timint_implicit.H"

#include "../linalg/linalg_utils.H"

namespace SCATRA
{
class TimIntGenAlpha : public virtual ScaTraTimIntImpl
{

public:

  /// Standard Constructor
  TimIntGenAlpha(Teuchos::RCP<DRT::Discretization>       dis,
      Teuchos::RCP<LINALG::Solver>            solver,
      Teuchos::RCP<Teuchos::ParameterList>    params,
      Teuchos::RCP<Teuchos::ParameterList>    extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntGenAlpha();

  /// initialize time integration scheme
  virtual void Setup();

  /// Print information about current time step to screen (reimplementation for Gen. Alpha)
  virtual void PrintTimeStepInfo()
  {
    if (myrank_==0)
        printf("\nTIME: %11.4E/%11.4E  DT = %11.4E  %s(a_F=%3.2f | a_M=%3.2f | gamma=%3.2f) STEP = %4d/%4d\n",
               time_,maxtime_,dta_,MethodTitle().c_str(),alphaF_,alphaM_,gamma_,step_,stepmax_);
    return;
  }

  /// compute values at intermediate time steps
  void ComputeIntermediateValues();

  /// compute values at the interior of the elements (required for hdg)
  virtual void ComputeInteriorValues(){return;};

  ///  compute scalar time derivative
  void ComputeTimeDerivative();

  ///  compute scalar time derivate parameters of the input voltage
  virtual void ComputeTimeDerivPot0(const bool init){return;};

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  virtual void Update(const int num=0);

  /// read restart data
  virtual void ReadRestart(const int step,Teuchos::RCP<IO::InputControl> input = Teuchos::null);

  // routine to return scalar field phi at time step n-1
  Teuchos::RCP<Epetra_Vector> Phinm() { return Teuchos::null; }

  /// routine to return scalar field phi at time step n+alpha_F
  Teuchos::RCP<Epetra_Vector> Phiaf() { return phiaf_; }

  //! return scalar field phi at time n+alpha_F (gen-alpha) or n+1 (otherwise)
  Teuchos::RCP<Epetra_Vector> Phiafnp() { return phiaf_; }

  /// routine to return scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phiam() { return phiam_; }

  /// routine to return time derivative of scalar field phi at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Phidtam() { return phidtam_; }

  /// routine to return fine-scale scalar field fsphi at time step n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> FsPhi() {
    if (Sep_ != Teuchos::null)
    Sep_->Multiply(false,*phiaf_,*fsphiaf_);
    return fsphiaf_; };

  /// routine to return time integration specific parameters
  Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList()
  {
    Teuchos::RCP<Teuchos::ParameterList> timeparams;
    timeparams = Teuchos::rcp(new Teuchos::ParameterList());
    timeparams->set("using stationary formulation",false);
    timeparams->set("using generalized-alpha time integration",true);
    timeparams->set("total time",time_-(1-alphaF_)*dta_);
    timeparams->set("time factor",genalphafac_*dta_);
    timeparams->set("alpha_F",alphaF_);
    return timeparams;
  }


protected:

  /// don't want = operator and cctor
  TimIntGenAlpha operator = (const TimIntGenAlpha& old);

  /// copy constructor
  TimIntGenAlpha (const TimIntGenAlpha& old);

  /// set time parameter for element evaluation
  virtual void SetElementTimeParameter(bool forcedincrementalsolver=false) const;

  //! set time for evaluation of Neumann boundary conditions
  virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

  /// set time parameter for element evaluation applying an iterative backward Euler scheme
  virtual void SetElementTimeParameterBackwardEuler() const;

  //! calculate consistent initial scalar time derivatives in compliance with initial scalar field
  void CalcInitialTimeDerivative();

  /// Set the part of the righthandside belonging to the last timestep.
  virtual void SetOldPartOfRighthandside();

  /// do explicit predictor step (-> better starting value for nonlinear solver)
  void ExplicitPredictor() const;

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
  virtual void OutputRestart() const;

  /// return the right time-scaling-factor for the true residual
  double ResidualScaling() const { return 1.0/(dta_*genalphafac_); }

  /// scalar at time n+alpha_F and n+alpha_M
  Teuchos::RCP<Epetra_Vector>  phiaf_;
  Teuchos::RCP<Epetra_Vector>  phiam_;

  /// scalar time derivative at time n+alpha_M
  Teuchos::RCP<Epetra_Vector>  phidtam_;

  /// fine-scale part at time n+alpha_F
  Teuchos::RCP<Epetra_Vector>  fsphiaf_;

  /// time factors for generalized-alpha time integration
  double alphaM_;
  double alphaF_;
  double gamma_;
  double genalphafac_;
}; // class TimIntGenAlpha

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_GENALPHA_H
