/*!----------------------------------------------------------------------
\file scatra_timint_implicit.H
\brief Associated with control routine for con.-diff.(in)stat. solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Georg Bauer
            bauer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef SCATRA_TIMINT_IMPLICIT_H
#define SCATRA_TIMINT_IMPLICIT_H

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/linalg_solver.H"
#include "../drt_lib/linalg_utils.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_io/io.H"
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_inpar/inpar_scatra.H"

//forward declarations
class PeriodicBoundaryConditions;


using namespace std;
using namespace Teuchos;

/// Algorithm and time integration classes for scalar transport problems
namespace SCATRA{

/*!
\brief implicit time integration for scalar transport problems

\author gjb (bauer@lnm.mw.tum.de)
*/
class ScaTraTimIntImpl
{
  friend class ScaTraResultTest;

public:

  /// Standard Constructor
  ScaTraTimIntImpl(RCP<DRT::Discretization>       dis,
                   RCP<LINALG::Solver>            solver,
                   RCP<ParameterList>             params,
                   RCP<ParameterList>             xparams,
                   RCP<IO::DiscretizationWriter>  output);

  /// Destructor
  virtual ~ScaTraTimIntImpl();

  /// Provide enum of time integration scheme
  enum INPAR::SCATRA::TimeIntegrationScheme MethodName(){return timealgo_;}

  /// Provide title of time integration scheme
  std::string MethodTitle(){return MapTimIntEnumToString(MethodName());}

  /// Print information about current time step to screen
  virtual void PrintTimeStepInfo()
  {
    if (myrank_==0)
        printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d \n",
               time_,maxtime_,dta_,MethodTitle().c_str(),step_,stepmax_);
    return;
  }

  /// do time integration (time loop)
  void TimeLoop();

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /// do explicit predictor step (-> better starting value for nonlinear solver)
  virtual void ExplicitPredictor() = 0;

  /// contains the nonlinear iteration loop
  void NonlinearSolve();

  /// check convergence (or divergence) of nonlinear iteration
  bool AbortNonlinIter(const int itnum, const int itmax, const double ittol, double& actresidual);

  /// prepare Krylov space projection
  void PrepareKrylovSpaceProjection();

  /// solve linear system
  void Solve();

  /// call elements to calculate system matrix and rhs and assemble
  void AssembleMatAndRHS();

  /// potential residual scaling and potential addition of Neumann terms
  void ScalingAndNeumann();

  /// update the solution after convergence of the nonlinear iteration.
  virtual void Update() = 0;

  /// update configuration and output to file
  void Output();

  /// read restart data
  virtual void ReadRestart(int step)=0;

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  /// note: VM3 solver still needs an explicit toggle vector for construction
  const Teuchos::RCP<const Epetra_Vector> DirichletToggle();

  /// prepare AVM3-based scale separation
  void AVM3Preparation();

  /// AVM3-based scale separation
  virtual void AVM3Separation()=0;

  /// scaling of AVM3-based subgrid-diffusivity matrix
  void AVM3Scaling(ParameterList& eleparams);

  /// set the velocity field (zero or field by function)
  void SetVelocityField();

  /// set the convective velocity + subgrid velocity/viscosity field (Navier-Stokes solution)
  void SetVelocityField(
      Teuchos::RCP<const Epetra_Vector> fluidvel, ///< fluid velocity/press. vector
      Teuchos::RCP<const Epetra_Vector> fluidacc, ///< fluid acceleration vector
      Teuchos::RCP<const DRT::DofSet> dofset,
      Teuchos::RCP<DRT::Discretization> fluiddis);    ///< fluid discretization

 /// apply moving mesh data
 void ApplyMeshMovement(
     Teuchos::RCP<const Epetra_Vector> dispnp,     ///< displacement vector
     Teuchos::RCP<DRT::Discretization> fluiddis ); ///< fluid discretization

  /// set the initial scalar field phi
  void SetInitialField(
      const INPAR::SCATRA::InitialField, ///< type of initial field
      const int startfuncno              ///< number of spatial function
      );

  /// set initial thermodynamic pressure for low-Mach-number flow
  void SetInitialThermPressure();

  /// compute initial thermodyn. pressure time derivative for low-Mach-number flow
  void ComputeInitialThermPressureDeriv();

  /// compute initial total mass in domain for low-Mach-number flow
  void ComputeInitialMass();

  /// predict thermodynamic pressure and time derivative for low-Mach-number flow
  virtual void PredictThermPressure() = 0;

  /// compute values at intermediate time steps (required for generalized-alpha)
  virtual void ComputeIntermediateValues() = 0;

  /// compute values of thermodynamic pressure at intermediate time steps
  virtual void ComputeThermPressureIntermediateValues() = 0;

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  virtual void ComputeThermPressure()=0;

  /// compute thermodyn. press. from mass cons. in domain for low-Mach-number flow
  void ComputeThermPressureFromMassCons();

  /// compute density field from concentration(s) for ELCH natural convection
  void ComputeDensity(const double density);

  /// check convergence for low-Mach-number flow
  bool LomaConvergenceCheck(int itnum, int itmax, const double ittol);

  ///  update time derivative
  virtual void UpdateTimeDerivative() = 0;

  ///  update time derivative of thermodynamic pressure after solution
  virtual void UpdateThermPressureTimeDerivative() = 0;

  /// update thermodynamic pressure and time derivative for low-Mach-number flow
  virtual void UpdateThermPressure() = 0;

  /// update density field for ELCH natural convection
  virtual void UpdateDensityElch() = 0;

  /// output of mean values of scalar(s)
  void OutputMeanScalars();

  /// routine to return scalar field phi at time step n+1
  Teuchos::RCP<Epetra_Vector> Phinp() { return phinp_; }

  /// routine to return scalar field phi at time step n
  Teuchos::RCP<Epetra_Vector> Phin() { return phin_; }

  /// routine to return time derivative of scalar field phi at time step n+1
  Teuchos::RCP<Epetra_Vector> Phidtnp() { return phidtnp_; }

  /// routine to return scalar field phi at time step n+alpha_F
  virtual Teuchos::RCP<Epetra_Vector> Phiaf() = 0;

  /// routine to return scalar field phi at time step n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> Phiam() = 0;

  /// routine to return time derivative of scalar field phi at time step n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> Phidtam() = 0;

  /// routine to return scalar field convective velocity at time step n+1
  Teuchos::RCP<Epetra_MultiVector> ConVel() { return convel_; }

  /// output of all electrode information to screen and/or file (default: 2x true)
  void OutputElectrodeInfo(bool printtoscreen = true, bool printtofile = true);

  /// output of electrode information for a single boundary condition
  void OutputSingleElectrodeInfo(
      DRT::Condition* condition,///< pointer to the condition
      const int  condid,        ///< condition id
      const bool printtoscreen, ///< flag for output to screen
      const bool printtofile,   ///< flag for output to file
      double& currentsum,       ///< net current over all boundaries (out)
      double& pottangent,       ///< tangent of current w.r.t. electrode potential (out)
      double& currresidual      ///< neg. residual of current equation (out)
      );

//REINHARD
  /// preliminary: reinitialize level set function
  void Reinitialize();
//end REINHARD

//REINHARD
  /// preliminary: calculate mass loss for level set function
  void CalculateMassLoss();
//end REINHARD

//REINHARD
  /// preliminary: set velocity field for level set functon
  void SetReinitVelocityField();
//end REINHARD

  /// routine to return density field at time n+1 for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchNp() { return elchdensnp_; }

  /// routine to return density field at time n for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchN() { return elchdensn_; }

  /// routine to return density field at time n-1 for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchNm() { return elchdensnm_; }

  /// routine to return trueresidual at time step n+1 for low-Mach-number flow
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }

  /// routine to return scalar increment at time step n+1 for low-Mach-number flow
  Teuchos::RCP<Epetra_Vector> ScalIncNp() { return phiincnp_; }

  /// routine to return discretization
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

  /// routine to return thermo. press. at time step n+1 for low-Mach-number flow
  double ThermPressNp() { return thermpressnp_; }

  /// routine to return thermo. press. at time step n for low-Mach-number flow
  double ThermPressN() { return thermpressn_; }

  /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
  virtual double ThermPressAf() = 0;

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  virtual double ThermPressAm() = 0;

  /// routine to return time der. of thermo. press. at time step n+1 for low-Mach-number flow
  double ThermPressDtNp() { return thermpressdtnp_; }

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  virtual double ThermPressDtAm() = 0;

  /// return current time value
  double Time() const { return time_; }

  /// return current step number
  int Step() const { return step_; }

  /// total number of time steps
  int NStep() const { return stepmax_;}

  /// return number of transported scalars
  int NumScal() const { return numscal_; }

  /// return system matrix down-casted as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }

  /// return system matrix downcasted as block sparse matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }

  /// calculate error compared to analytical solution
  void EvaluateErrorComparedToAnalyticalSol();

  /// routine to calculate mass/heat flux vector field in computational domain
  Teuchos::RCP<Epetra_MultiVector> CalcFluxInDomain(
      const INPAR::SCATRA::FluxType fluxtype ///< fluxtype enum
  );

  /// routine to calculate mass/heat normal flux at specified boundaries
  Teuchos::RCP<Epetra_MultiVector> CalcFluxAtBoundary(std::vector<string>& condnames);

  /// returns map extractor used for convergence check either in ELCH or LOMA case
  LINALG::MapExtractor& Splitter() { return splitter_; }

  /// computes conductivity of an electrolyte solution (ELCH)
  double ComputeConductivity();

  /// apply galvanostatic control (update electrode potential) (ELCH)
  bool ApplyGalvanostaticControl();

protected:

  /// don't want = operator and cctor
  //ScaTraTimIntImpl operator = (const ScaTraTimIntImpl& old);
  /// copy constructor
  ScaTraTimIntImpl (const ScaTraTimIntImpl& old);

  /// setup of natural convection (ELCH)
  void SetupElchNatConv();

  /// helper function to get algorithm title
  std::string MapTimIntEnumToString
  (
    const enum INPAR::SCATRA::TimeIntegrationScheme term  ///< the enum
  );

  /// increment time and step value
  void IncrementTimeAndStep()
  {
    step_ += 1;
    time_ += dta_;
  }

  /// calculate initial value for phidtn_
  void CalcInitialPhidt();

  /// initialization procedure before the first time step is done
  virtual void PrepareFirstTimeStep() = 0;

  /// Set the part of the righthandside belonging to the last timestep.
  virtual void SetOldPartOfRighthandside() = 0;

  /// set time for evaluation of Neumann boundary conditions
  virtual void SetTimeForNeumannEvaluation(ParameterList& params) = 0;

  /// add actual Neumann loads multipl. with time factor to the residual
  virtual void AddNeumannToResidual() = 0;

  /// add parameters specific for time-integration scheme
  virtual void AddSpecificTimeIntegrationParameters(ParameterList& params) = 0;

  /// Apply Dirichlet boundary conditions on provided state vector
  void ApplyDirichletBC
  (
    const double time,  ///< evaluation time
    Teuchos::RCP<Epetra_Vector> phinp, ///< transported scalar(s) (may be = null)
    Teuchos::RCP<Epetra_Vector> phidt  ///< first time derivative (may be = null)
  );

  /// Apply Neumann boundary conditions
  void ApplyNeumannBC
  (
    const double time,  ///< evaluation time
    const Teuchos::RCP<Epetra_Vector> phinp,   ///< transported scalar(s) state
    Teuchos::RCP<Epetra_Vector> neumann_loads ///< Neumann loads
  );

  /// export node-based multivector to column map and add it to a parameter list
  void AddMultiVectorToParameterList(
      Teuchos::ParameterList& p,           ///< parameter list
      const std::string       name,        ///< naming string for added vector
      Teuchos::RCP<Epetra_MultiVector> vec ///< Epetra_MultiVector to be added
  );

  /// write state vectors (phinp and convective velocity) to BINIO
  void OutputState();

  /// write additional data required for restart
  virtual void OutputRestart() = 0;

  /// routine to calculate mass/heat flux vectors (in domain / at boundary)
  Teuchos::RCP<Epetra_MultiVector> CalcFlux();

  /// compute contribution of electrode kinetics to eq. system
  void EvaluateElectrodeKinetics(
      RCP<LINALG::SparseOperator> matrix,
      RCP<Epetra_Vector>          rhs);

  /// compute potential Neumann inflow
  void ComputeNeumannInflow(
      RCP<LINALG::SparseOperator> matrix,
      RCP<Epetra_Vector>          rhs);

  /// return the right time-scaling-factor for the true residual
  virtual double ResidualScaling() const = 0;

  /// write flux vectors to BINIO
  void OutputFlux();

  //! compute outward pointing unit normal vectors at given bc's
  RCP<Epetra_MultiVector> ComputeNormalVectors(const vector<string>& condnames);

  //! update time-dependent electrode state variables at the end of an time step
  virtual void ElectrodeKineticsTimeUpdate(const bool init = false) = 0;

  /// the scalar transport discretization
  RCP<DRT::Discretization>   discret_;
  /// the solver
  RCP<LINALG::Solver>        solver_;
  /// the parameter list
  RCP<ParameterList>         params_;
  /// the parameter list containing extra parameters (application dependent)
  RCP<ParameterList>         extraparams_;
  /// the discretization writer
  RCP<IO::DiscretizationWriter>  output_;
  /// processor id
  int                        myrank_;

  /// actual time
  double time_;

  /// actual step number
  int step_;

  /// flag defining type of scalar transport problem
  string prbtype_;

  /// solvertype and flags for nonlinear (always incremental) and
  /// (linear) incremental solver
  INPAR::SCATRA::SolverType solvtype_;
  bool nonlinear_;
  bool incremental_;

  /// flag for potential inclusion of reaction term
  string reaction_;

  /// flag for Eulerian or ALE formulation of equation(s)
  bool isale_;

  /// flag for type of scalar transport problem
  const INPAR::SCATRA::ScaTraType scatratype_;

  /// flag controlling level set reinitialization strategy
  INPAR::SCATRA::ReinitializationAction reinitaction_;

  /// flag controlling level set mass calculation strategy
  INPAR::SCATRA::MassCalculation masscalc_;

  bool reinitswitch_;

  /// maximum number of steps
  int stepmax_;

  /// maximum simulation time
  double maxtime_;

  /// used time integration scheme
  INPAR::SCATRA::TimeIntegrationScheme timealgo_;

  /// write results every upres_ steps
  int upres_;

  /// write restart data every uprestart_ steps
  int uprestart_;

  /// flag for writing the flux vector fields
  INPAR::SCATRA::FluxType writeflux_;

  /// flag for printing out mean values of scalars
  const bool outmean_;

  /// time step size
  double dta_;

  /// time step size
  double dtp_;

  /// type of velocity field
  const INPAR::SCATRA::VelocityField cdvel_;

  /// maximum number of nonlinear iterations
  int itemax_;

  /// parameters for sampling/dumping period
  int samstart_;
  int samstop_;
  int dumperiod_;

  /// form of convective term
  INPAR::SCATRA::ConvForm convform_;

  /// flag for potential Neumann inflow
  bool neumanninflow_;

  /// fine-scale subgrid-viscosity flag
  INPAR::SCATRA::FSSUGRDIFF fssgd_;

  /// ELCH-specific parameter F/RT (default value for the temperature is 298K)
  double frt_;

  /// LOMA-specific parameter: flag for turbulence model
  bool turbmodel_;

  /// LOMA-specific parameters: Smagorinsky constant and turbulent Prandtl number
  double Cs_;
  double tpn_;

  /// LOMA-specific parameter: thermodynamic pressure at n+1 and n
  double thermpressnp_;
  double thermpressn_;

  /// LOMA-specific parameter: time deriv. of thermodynamic pressure at n+1 and n
  double thermpressdtnp_;
  double thermpressdtn_;

  /// LOMA-specific parameter: initial mass in domain
  double initialmass_;

  /// system matrix (either sparse matrix or block sparse matrix)
  RCP<LINALG::SparseOperator> sysmat_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// maps for extracting Dirichlet and free DOF sets
  RCP<LINALG::MapExtractor> dbcmaps_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// The residual vector
  RCP<Epetra_Vector>    residual_;

  /// phi at times n+1 and n
  RCP<Epetra_Vector>    phinp_;
  RCP<Epetra_Vector>    phin_;

  /// time derivative of phi at times n+1 and n
  RCP<Epetra_Vector>  phidtnp_;
  RCP<Epetra_Vector>  phidtn_;

  /// true (rescaled) residual vector without zeros at Dirichlet conditions
  RCP<Epetra_Vector>    trueresidual_;

  /// density (nodebased!) at times n+1, n and n-1 for ELCH natural convection algorithm
  RCP<Epetra_Vector>    elchdensnp_;
  RCP<Epetra_Vector>    elchdensn_;
  RCP<Epetra_Vector>    elchdensnm_;

  /// scalar and velocity increment at time n+1
  RCP<Epetra_Vector>       phiincnp_;
  //RCP<Epetra_MultiVector>  velincnp_;

  /// histvector --- a linear combination of phinm, phin (BDF)
  ///                or phin, phidtn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  /// convective velocity
  RCP<Epetra_MultiVector>    convel_;

  /// acceleration and pressure required for computation of subgrid-scale velocity
  RCP<Epetra_MultiVector>    accpre_;

  /// displacements at time n+1 for ALE formulation
  RCP<Epetra_MultiVector>    dispnp_;

  /// Nonlinear iteration increment vector
  RCP<Epetra_Vector>    increment_;

  /// subgrid-diffusivity(-scaling) vector
  RCP<Epetra_Vector>    subgrdiff_;

  /// only necessary for AVM3: fine-scale subgrid-diffusivity matrix
  RCP<LINALG::SparseMatrix> sysmat_sd_;

  /// only necessary for AVM3: scale-separation matrix
  RCP<LINALG::SparseMatrix> Sep_;

  /// only necessary for AVM3: normalized fine-scale subgrid-viscosity matrix
  RCP<LINALG::SparseMatrix> Mnsv_;

  /// Extractor used for convergence check either in ELCH or LOMA case
  LINALG::MapExtractor splitter_;

  /// Periodic Boundary Conditions. Guess what.
  RCP<PeriodicBoundaryConditions> pbc_;

  /// connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > >     pbcmapmastertoslave_;

  /// time measurement
  double dtele_;
  double dtsolve_;

  /// error file handle
  FILE* errfile_;

  /// vector for statistical evaluation of normal fluxes
  RCP<Epetra_SerialDenseVector> sumnormfluxintegral_;

  /// number of transported scalars
  int numscal_;

  /// boolean indicating if initial velocity field has been set
  bool initialvelset_;

  /// the last step number when fluxes have been computed
  int lastfluxoutputstep_;

  /// number of iterations in galvanostatic mode (ELCH)
  int gstatnumite_;

  /// value of electric potential increment in galvanostatic mode (ELCH)
  double gstatincrement_;

  /// unit outer normal vector field for flux output
  RCP<Epetra_MultiVector> normals_;

  /// densification coefficients for ELCH natural convection
  vector<double> densific_;

  /// initial concentrations of the electrolyte solution (ELCH)
  vector<double> c0_;

  /// flag for Krylov space projection
  bool                             project_;
  /// integrated basis vector (for each undetermined component)
  Teuchos::RCP<Epetra_MultiVector> w_;
  /// basis vector for nullspace (vector of ones for each undetermined component)
  Teuchos::RCP<Epetra_MultiVector> c_;

}; // class ScaTraTimIntImpl

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_IMPLICIT_H
#endif  // #ifdef CCADISCRET
