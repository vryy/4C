/*!-----------------------------------------------------------------------------------------------*
\file scatra_timint_tg.H

\brief  explicit Taylor Galerkin time integration scheme (TG) and implicit Characteristic Galerkin scheme (ICG)
        just for the pure 1st order transport equation, not yet implemented for convections-diffusion equation

<pre>
Maintainer: Benedikt Schott
            schott@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15241
</pre>
 *------------------------------------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_TG_H
#define SCATRA_TIMINT_TG_H

#include "scatra_timint_implicit.H"

#include "../linalg/linalg_utils.H"


namespace SCATRA
{

class TimIntTaylorGalerkin : public ScaTraTimIntImpl
{
public:

    /// Standard Constructor
    TimIntTaylorGalerkin(Teuchos::RCP<DRT::Discretization>       dis,         ///< discretization
                         Teuchos::RCP<LINALG::Solver>            solver,      ///< solver
                         Teuchos::RCP<Teuchos::ParameterList>    params,      ///< time integration parameters
                         Teuchos::RCP<Teuchos::ParameterList>    extraparams, ///< extra parameters
                         Teuchos::RCP<IO::DiscretizationWriter>  output);     ///< output discretization writer

    /// Destructor
    virtual ~TimIntTaylorGalerkin();

    /// Print information about current time step to screen (reimplementation for OST)
    virtual void PrintTimeStepInfo()
    {
      if (myrank_==0)
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d \n",
            time_,maxtime_,dta_,MethodTitle().c_str(),step_,stepmax_);
      return;
    }

    /// predict thermodynamic pressure and time derivative for low-Mach-number flow
    void PredictThermPressure();

    /// compute values at intermediate time steps (required for generalized-alpha)
    void ComputeIntermediateValues(){return;};

    /// compute values of thermodynamic pressure at intermediate time steps (required for generalized-alpha)
    void ComputeThermPressureIntermediateValues(){return;};

    /// compute thermodynamic pressure and time derivative for low-Mach-number flow
    void ComputeThermPressure();

    ///  compute scalar time derivative
    void ComputeTimeDerivative();

    ///  compute scalar time derivate parameters of the input voltage
    void ComputeTimeDerivPot0(const bool init)
    {
      dserror("Double layer charging is not implemented for BDF2! You have to use one-step-theta time integration scheme");
      return;
    };

    ///  compute time derivative of thermodynamic pressure
    void ComputeThermPressureTimeDerivative();

    /// Update the solution after convergence of the nonlinear iteration.
    /// Current solution becomes old solution of next timestep.
    void Update(const int num=0);

    /// update the solution after reinitialization
    void UpdateReinit();

    /// update thermodynamic pressure and time derivative for low-Mach-number flow
    void UpdateThermPressure();

    /// update density field for ELCH natural convection
    void UpdateDensityElch();

    /// read restart data
    void ReadRestart(int step);

    // routine to return scalar field phi at time step n-1
    Teuchos::RCP<Epetra_Vector> Phinm() { return phinm_; }

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

    /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
    double ThermPressAf(){return 1.0;}

    /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
    double ThermPressAm(){return 1.0;}

    /// routine to return time derivative of thermo. press. at time step n+alpha_F for low-Mach-number flow
    double ThermPressDtAf(){return 0.0;}

    /// routine to return time derivative of thermo. press. at time step n+alpha_M for low-Mach-number flow
    double ThermPressDtAm(){return 0.0;}

    /// redistribute the scatra discretization and vectors according to nodegraph
    void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph);

    /// routine to return time integration specific parameters
    Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList()
    {dserror("Not yet implemented!"); return Teuchos::null;}

protected:

    /// don't want = operator and cctor
    TimIntTaylorGalerkin operator = (const TimIntTaylorGalerkin& old);

    /// copy constructor
    TimIntTaylorGalerkin (const TimIntTaylorGalerkin& old);

    /// Initialization procedure before the first time step is done
    void PrepareFirstTimeStep();

    /// Set the part of the righthandside belonging to the last timestep.
    void SetOldPartOfRighthandside(){;}

    /// do explicit predictor step (-> better starting value for nonlinear solver)
    void ExplicitPredictor();

    /// set time for evaluation of Neumann boundary conditions
    void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params);

    /// add actual Neumann loads with time factor
    void AddNeumannToResidual();

    /// add actual Neumann inflow loads
    void ComputeNeumannInflowTG(
        Teuchos::RCP<LINALG::SparseOperator> matrix,
        Teuchos::RCP<Epetra_Vector>          rhs
        );

    /// AVM3-based scale separation
    void AVM3Separation();

    /// dynamic Smagorinsky model
    void DynamicComputationOfCs() {dserror("not implemented!"); return;};
    
    // dynamic Vreman model
    void DynamicComputationOfCv() {dserror("not implemented!"); return;};

    /// add parameters specific for time-integration scheme
    void AddSpecificTimeIntegrationParameters(Teuchos::ParameterList& params);

    /// write additional data required for restart
    void OutputRestart();

    /// return the right time-scaling-factor for the true residual
    double ResidualScaling() const { return 1.0; }

    /// update time-dependent electrode state variables at the end of an time step
    void ElectrodeKineticsTimeUpdate();

    /// set old part of RHS for galvanostatic equation
    void ElectrodeKineticsSetOldPartOfRHS();

private:

    /// fine-scale solution vector at time n+1
    Teuchos::RCP<Epetra_Vector>  fsphinp_;

    /// solution vector phi at time n-1
    Teuchos::RCP<Epetra_Vector>    phinm_;


}; // class TimIntTaylorGalerkin

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_TG_H
