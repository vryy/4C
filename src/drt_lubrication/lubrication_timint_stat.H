/*--------------------------------------------------------------------------*/
/*!
\file lubrication_timint_stat.H

\brief solution algorithm for stationary problems

<pre>
Maintainer: Andy Wirtz
            wirtz@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089-289-15270
</pre>
*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_STAT_H_
#define SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_STAT_H_


#include "lubrication_timint_implicit.H"

#include "../linalg/linalg_utils.H"

namespace LUBRICATION
{
class TimIntStationary : public TimIntImpl
{

public:

  /// Standard Constructor
  TimIntStationary(Teuchos::RCP<DRT::Discretization>      dis,
      Teuchos::RCP<LINALG::Solver>            solver,
      Teuchos::RCP<Teuchos::ParameterList>    params,
      Teuchos::RCP<Teuchos::ParameterList>    extraparams,
      Teuchos::RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~TimIntStationary();

  /// initialize time integration scheme
  virtual void Init();

  /// compute values at intermediate time steps (required for generalized-alpha)
  void ComputeIntermediateValues(){return;};

  ///  compute lubrication time derivative
  void ComputeTimeDerivative(){return;};

  ///  compute lubrication time derivate parameters of the input voltage
  void ComputeTimeDerivPot0(const bool init){return;};

  /// Update the solution after convergence of the nonlinear iteration.
  /// Current solution becomes old solution of next timestep.
  virtual void Update(const int num=0);

  /// read restart data
  void ReadRestart(int step);

  // routine to return pressure pre at time step n-1
  Teuchos::RCP<Epetra_Vector> Prenm() { return Teuchos::null; }

  /// routine to return pressure field pre at time step n+alpha_F
  Teuchos::RCP<Epetra_Vector> Preaf(){return Teuchos::null;}

  /// routine to return pressure field pre at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Pream(){return Teuchos::null;}

  /// routine to return time derivative of pressure field pre at time step n+alpha_M
  Teuchos::RCP<Epetra_Vector> Predtam(){return Teuchos::null;}

  /// routine to return fine-scale pressure field fspre
  Teuchos::RCP<Epetra_Vector> FsPre() {
    if (Sep_ != Teuchos::null)
    Sep_->Multiply(false,*prenp_,*fsprenp_);
    return fsprenp_; };

  /// routine to return time integration specific parameters
  Teuchos::RCP<Teuchos::ParameterList> LubricationTimeParameterList()
  {dserror("Not yet implemented!"); return Teuchos::null;}


protected:

  /// don't want = operator and cctor
  TimIntStationary operator = (const TimIntStationary& old);

  /// copy constructor
  TimIntStationary (const TimIntStationary& old);

  /// set time parameter for element evaluation
  void SetElementTimeParameter(bool forcedincrementalsolver=false) const;

  //! set time for evaluation of Neumann boundary conditions
  void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

  //! calculate consistent initial conditions in compliance with initial pressure field
  //! this is not necessary for stationary calculations
  void CalcInitialTimeDerivative(){return;};

  /// Set the part of the righthandside belonging to the last timestep.
  void SetOldPartOfRighthandside();

  /// do explicit predictor step (nothing to predict for stationary problems!)
  void ExplicitPredictor(){return;};

  /// add actual Neumann loads with time factor
  void AddNeumannToResidual();

  /// AVM3-based scale separation
  void AVM3Separation();

  /// add parameters specific for time-integration scheme
  void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver=false);

  /// dynamic Smagorinsky model
  void DynamicComputationOfCs() {dserror("no turbulence in stationary flows!"); return;};

  /// dynamic Vreman model
  void DynamicComputationOfCv() {dserror("no turbulence in stationary flows!"); return;};

  /// write additional data required for restart
  void OutputRestart();

  /// return the right time-scaling-factor for the true residual
  double ResidualScaling() const { return 1.0; }

private:

  /// fine-scale solution vector at time n+1
  Teuchos::RCP<Epetra_Vector>  fsprenp_;


}; // class TimIntStationary

} // namespace LUBRICATION

#endif /* SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_STAT_H_ */
