/*!----------------------------------------------------------------------
\file scatra_timint_implicit.H
\brief Associated with control routine for con.-diff.(in)stat. solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o implicit characteristic Galerkin (ICG) time-integration scheme (level-set transport)

     o explicit taylor galerkin (TG) time-integration schemes (level-set transport)

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

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_CrsGraph.h>
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_fluid.H"


//forward declarations
class PeriodicBoundaryConditions;
namespace DRT
{
  class Discretization;
  class DofSet;
  class Condition;
}
namespace IO
{
  class DiscretizationWriter;
}
namespace LINALG
{
  class Solver;
  class SparseMatrix;
  class MapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
}

namespace FLD
{
  class Meshtying;
}


/// Algorithm and time integration classes for scalar transport problems
namespace SCATRA{

/*!
\brief implicit time integration for scalar transport problems

\author gjb (bauer@lnm.mw.tum.de)
*/
class ScaTraTimIntImpl
{

public:

  /// Standard Constructor
  ScaTraTimIntImpl(Teuchos::RCP<DRT::Discretization>       dis,
                   Teuchos::RCP<LINALG::Solver>            solver,
                   Teuchos::RCP<Teuchos::ParameterList>    params,
                   Teuchos::RCP<Teuchos::ParameterList>    xparams,
                   Teuchos::RCP<IO::DiscretizationWriter>  output);

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
  }

  /// print information about current pseudo time step to screen
  virtual void PrintPseudoTimeStepInfoReinit()
  {
    if (myrank_==0)
      printf("TIMESTEPFACTOR: %11.4E      %s          PSEUDOSTEP = %4d/%4d \n",
          dta_,MethodTitle().c_str(),step_,stepmax_);
  }

  /// increment time and step value
  void SetTimeStep(double time = 0.0, int step = 0)
  {
    time_ = time;
    step_ = step;
  }

  /// evaluate integrated error norm for reinitialization stop criterion
  double EvaluateGradientNormError();

  /// do time integration (time loop)
  void TimeLoop();

  /// do time integration for reinitialization (time loop)
  void TimeLoop_Reinit();

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /// calculate initial value for phidtn_ (only necessary for OnestepTheta and GenAlpha)
  void CalcInitialPhidt();

  /// calculate initial value for phidtn_ (assembly of system)
  void CalcInitialPhidtAssemble();

  /// calculate initial value for phidtn_ (solver part)
  void CalcInitialPhidtSolve();

  /// do explicit predictor step (-> better starting value for nonlinear solver)
  virtual void ExplicitPredictor() = 0;

  /// general solver call for coupled algorithms (decides if linear/nonlinear internally)
  void Solve();

  /// update the solution after convergence of the nonlinear iteration.
  virtual void Update() = 0;

  /// update the solution after reinitialization
  virtual void UpdateReinit() = 0;

  /// update configuration and output to file
  void Output();

  /// output to file for reinitialized scatra field
  void OutputReinit(const int globalstep, const double globaltime);

  /// output reinitialization steps
  void OutputReinitializationSteps();

  /// write state vectors to Gmsh postprocessing files
  void OutputToGmshReinitializationSteps(
      const int reinit_step,
      const double reinit_pseudo_time
  ) const;

  /// write state vectors (phinp) to Gmsh postprocessing files
  void OutputToGmshReinit(
      const int globalstep,
      const double globaltime,
      const int step,
      const double time
  ) const;

  /// Access output object
  IO::DiscretizationWriter& DiscWriter() { return *output_; }

  /// read restart data
  virtual void ReadRestart(int step)=0;

  /// set the velocity field (zero or field by function)
  void SetVelocityField();

  /// set the convective velocity + subgrid velocity/viscosity field (Navier-Stokes solution)
  void SetVelocityField(
      Teuchos::RCP<const Epetra_Vector> confluidvel, ///< convective fluid velocity/press. vector
      Teuchos::RCP<const Epetra_Vector> fluidacc, ///< fluid acceleration vector
      Teuchos::RCP<const Epetra_Vector> fluidvel,  ///< fluid velocity vector
      Teuchos::RCP<const Epetra_Vector> fluidfsvel, ///< fluid fine-scale velocity vector
      Teuchos::RCP<const DRT::DofSet> dofset,
      Teuchos::RCP<DRT::Discretization> fluiddis);    ///< fluid discretization

  /// call different reinitialization techniques based on solving the sussman equation
  bool CallReinitialization();

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

  /// compute thermodynamic pressure and time derivative for low-Mach-number flow
  virtual void ComputeThermPressure()=0;

  /// compute thermodyn. press. from mass cons. in domain for low-Mach-number flow
  void ComputeThermPressureFromMassCons();

  /// compute density field from concentration(s) for ELCH natural convection
  void ComputeDensity();

  /// check convergence (only for low-Mach-number flow, for the time being)
  bool ConvergenceCheck(int itnum, int itmax, const double ittol);

  ///  compute time derivative
  virtual void ComputeTimeDerivative() = 0;

  ///  compute time derivative of thermodynamic pressure after solution
  virtual void ComputeThermPressureTimeDerivative() = 0;

  /// update thermodynamic pressure and time derivative for low-Mach-number flow
  virtual void UpdateThermPressure() = 0;

  /// update density field for ELCH natural convection
  virtual void UpdateDensityElch() = 0;

  /// iterative update of concentrations
  void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc);

  /// output of mean values of scalar(s)
  void OutputMeanScalars();

  /// routine to return scalar field phi at time step n+1
  Teuchos::RCP<Epetra_Vector> Phinp() { return phinp_; }

  /// routine to return scalar field phi at time step n
  Teuchos::RCP<Epetra_Vector> Phin() { return phin_; }

  // routine to return scalar field phi at time step n-1
  virtual Teuchos::RCP<Epetra_Vector> Phinm() = 0;

  /// routine to return time derivative of scalar field phi at time step n+1
  Teuchos::RCP<Epetra_Vector> Phidtnp() { return phidtnp_; }

  /// routine to return time derivative of scalar field phi at time step n
  Teuchos::RCP<Epetra_Vector> Phidtn() { return phidtn_; }

  /// routine to return scalar field phi at time step n+alpha_F
  virtual Teuchos::RCP<Epetra_Vector> Phiaf() = 0;

  /// routine to return scalar field phi at time step n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> Phiam() = 0;

  /// routine to return time derivative of scalar field phi at time step n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> Phidtam() = 0;

  /// routine to return scalar field convective velocity at time step n+1
  Teuchos::RCP<Epetra_MultiVector> ConVel() { return convel_; }

  /// routine to reset phi vector due to reinitialization for OST scheme
  virtual void SetPhin(Teuchos::RCP<Epetra_Vector> phireinitn) = 0;

  /// routine to set phi vectors due to reinitialization
  void SetPhinp(Teuchos::RCP<Epetra_Vector> phinp);

  /// routine to set phi vectors due to reinitialization
  void SetPhiReinit(Teuchos::RCP<Epetra_Vector> phi);

  /// routine to check for steady state due to reinitialization
  void CheckSteadyState(double&    Gradient_Error_old,
  		                bool&      STOP);

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
      double& currresidual,      ///< neg. residual of current equation (out)
      double& electrodesurface   ///< surface of the electrode
      );

  /// routine to return density field at time n+1 for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchNp() { return elchdensnp_; }

  /// routine to return density field at time n for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchN() { return elchdensn_; }

  /// routine to return density field at time n-1 for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchNm() { return elchdensnm_; }

  /// routine to return trueresidual at time step n+1 for low-Mach-number flow
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }

  /// routine to return discretization
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

  /// routine to return information on periodic bondary conditions (row map)
  Teuchos::RCP<std::map<int,std::vector<int> > > PBCmap() { return pbcmapmastertoslave_; };

  /// routine to return the periodic boundary condition member
  Teuchos::RCP<PeriodicBoundaryConditions> PBC() { return pbc_; }

  /// routine to return thermo. press. at time step n+1 for low-Mach-number flow
  double ThermPressNp() const { return thermpressnp_; }

  /// routine to return thermo. press. at time step n for low-Mach-number flow
  double ThermPressN() const { return thermpressn_; }

  /// routine to return thermo. press. at time step n+alpha_F for low-Mach-number flow
  virtual double ThermPressAf() = 0;

  /// routine to return thermo. press. at time step n+alpha_M for low-Mach-number flow
  virtual double ThermPressAm() = 0;

  /// routine to return time der. of thermo. press. at time step n+1 for low-Mach-number flow
  double ThermPressDtNp() const { return thermpressdtnp_; }

  /// routine to return time derivative of thermo. press. at time step n+alpha_F for low-Mach-number flow
  virtual double ThermPressDtAf() = 0;

  /// routine to return time derivative of thermo. press. at time step n+alpha_M for low-Mach-number flow
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
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

  /// return system matrix downcasted as block sparse matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

  /// return dof map of vector of unknowns
  Teuchos::RCP<const Epetra_Map> DofRowMap();

  /// return dof map of vector of unknowns for multiple dofsets
  Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds);

  /// return residual
  Teuchos::RCP<Epetra_Vector> Residual();

  /// calculate error compared to analytical solution
  void EvaluateErrorComparedToAnalyticalSol();

  /// routine to calculate mass/heat flux vector field in computational domain
  Teuchos::RCP<Epetra_MultiVector> CalcFluxInDomain(
      const INPAR::SCATRA::FluxType fluxtype ///< fluxtype enum
  );

  /// routine to calculate mass/heat normal flux at specified boundaries
  Teuchos::RCP<Epetra_MultiVector> CalcFluxAtBoundary(
      std::vector<string>& condnames,
      const bool writetofile);

  /// compute contribution of permeable surface/interface
  void SurfacePermeability(
      Teuchos::RCP<LINALG::SparseOperator> matrix,
      Teuchos::RCP<Epetra_Vector>          rhs);

  /// returns map extractor used for convergence check either in ELCH or LOMA case
  Teuchos::RCP<LINALG::MapExtractor> Splitter() { return splitter_; }

  /// return type of physics we are solving for
  INPAR::SCATRA::ScaTraType ScaTraType() const {return scatratype_;};

  /// routine to calculate mass/heat flux vectors (in domain / at boundary)
  Teuchos::RCP<Epetra_MultiVector> CalcFlux(const bool writetofile);

  /// return flux vector
  Teuchos::RCP<Epetra_MultiVector> Flux(){return flux_;};

  /// preparations for solve
  void PrepareLinearSolve();

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// return flag for incremental solver
  bool Incremental() { return incremental_; }

  /// redistribute the scatra discretization and vectors according to nodegraph
  virtual void Redistribute(const Teuchos::RCP<Epetra_CrsGraph> nodegraph);

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

  /// calculate phidt after reinitialization for OST scheme
  virtual void CalcPhidtReinit() = 0;

  /// calculate initial field for electric potential based on ion concentrations (ELCH)
  void CalcInitialPotentialField();

  /// initialization procedure before the first time step is done
  virtual void PrepareFirstTimeStep() = 0;

  /// Set the part of the righthandside belonging to the last timestep.
  virtual void SetOldPartOfRighthandside() = 0;

  /// set time for evaluation of Neumann boundary conditions
  virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) = 0;

  /// add actual Neumann loads multipl. with time factor to the residual
  virtual void AddNeumannToResidual() = 0;

  /// add parameters specific for time-integration scheme
  virtual void AddSpecificTimeIntegrationParameters(Teuchos::ParameterList& params) = 0;

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

  /// compute approximation for fluxes and add it to a parameter list
  void AddFluxApproxToParameterList(
      Teuchos::ParameterList& p,                  ///< parameter list
      const enum INPAR::SCATRA::FluxType fluxtype ///< switch for flux calculation
  );

  /// prepare Krylov space projection
  void PrepareKrylovSpaceProjection();

  /// potential residual scaling and potential addition of Neumann terms
  void ScalingAndNeumann();

  /// call elements to calculate system matrix and rhs and assemble
  void AssembleMatAndRHS();

  /// call elements to calculate system matrix and rhs at boundaries for Taylor-Galerkin methods
  void AssembleMatAndRHS_Boundary();

  /// apply Dirichlet boundary conditions
  void ApplyDirichletToSystem();

  /// output of all reinitialization information to screen and/or file
  void ReinitializeInfo(bool printtoscreen = true, bool printtofile = true); //schott

  /// add reinitialization specific parameters
  void AddReinitializationParameters(Teuchos::ParameterList& params);

  /// solve linear system
  void LinearSolve();

  /// contains the nonlinear iteration loop
  void NonlinearSolve();

  /// check convergence (or divergence) of nonlinear iteration
  bool AbortNonlinIter(
      const int itnum,
      const int itmax,
      const double ittol,
      const double abstolres,
      double& actresidual);

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  /// note: VM3 solver still needs an explicit toggle vector for construction
  const Teuchos::RCP<const Epetra_Vector> DirichletToggle();

  /// compute values at intermediate time steps (required for generalized-alpha)
  virtual void ComputeIntermediateValues() = 0;

  /// compute values of thermodynamic pressure at intermediate time steps
  /// (required for generalized-alpha)
  virtual void ComputeThermPressureIntermediateValues() = 0;

  /// prepare AVM3-based scale separation
  void AVM3Preparation();

  /// AVM3-based scale separation
  virtual void AVM3Separation()=0;

  /// scaling of AVM3-based subgrid-diffusivity matrix
  void AVM3Scaling(Teuchos::ParameterList& eleparams);

  /// is output needed for the current time step?
  bool DoOutput()
    {return ((step_%upres_==0 )or (step_%uprestart_==0));};

  /// do we need a statistical sampling for boundary flux at the current time step?
  bool DoBoundaryFluxStatistics()
    { return ((step_>=samstart_) and (step_<=samstop_)
        and ((writeflux_==INPAR::SCATRA::flux_total_boundary)
          or (writeflux_==INPAR::SCATRA::flux_diffusive_boundary)
          or (writeflux_==INPAR::SCATRA::flux_convective_boundary)));};

  /// write state vectors (phinp and convective velocity) to BINIO
  void OutputState();

  /// write state vectors (phinp and convective velocity) to Gmsh postprocessing files
  void OutputToGmsh(
      const int step,
      const double time
  ) const;

  /// write additional data required for restart
  virtual void OutputRestart() = 0;

  /// compute contribution of electrode kinetics to eq. system
  void EvaluateElectrodeKinetics(
      Teuchos::RCP<LINALG::SparseOperator> matrix,
      Teuchos::RCP<Epetra_Vector>          rhs);

  /// evaluate Neumann inflow boundary condition
  void ComputeNeumannInflow(
      Teuchos::RCP<LINALG::SparseOperator> matrix,
      Teuchos::RCP<Epetra_Vector>          rhs);

  /// evaluate boundary condition due to convective heat transfer
  void EvaluateConvectiveHeatTransfer(
      Teuchos::RCP<LINALG::SparseOperator> matrix,
      Teuchos::RCP<Epetra_Vector>          rhs);

  /// return the right time-scaling-factor for the true residual
  virtual double ResidualScaling() const = 0;

  /// write flux vectors to BINIO
  void OutputFlux(Teuchos::RCP<Epetra_MultiVector> flux);

  /// computes different conductivity expressions for electrolyte solutions (ELCH)
  Epetra_SerialDenseVector ComputeConductivity();

  /// apply galvanostatic control (update electrode potential) (ELCH)
  bool ApplyGalvanostaticControl();

  /// check for negative values of concentrations (ELCH)
  void CheckConcentrationValues(Teuchos::RCP<Epetra_Vector> vec);

  //! compute outward pointing unit normal vectors at given bc's
  Teuchos::RCP<Epetra_MultiVector> ComputeNormalVectors(const std::vector<std::string>& condnames);

  //! update time-dependent electrode state variables at the end of an time step
  virtual void ElectrodeKineticsTimeUpdate(const bool init = false) = 0;

  //! scale lines of linear system prior to solver call (ELCH)
  void ScaleLinearSystem();

  //! set magnetic field (ELCH)
  void SetMagneticField(
      const int funcno ///< number of spatial function that defines the vector field
      );

  /// the scalar transport discretization
  Teuchos::RCP<DRT::Discretization>   discret_;
  /// the solver
  Teuchos::RCP<LINALG::Solver>        solver_;
  /// the parameter list
  Teuchos::RCP<Teuchos::ParameterList> params_;
  /// the parameter list containing extra parameters (application dependent)
  Teuchos::RCP<Teuchos::ParameterList> extraparams_;
  /// the discretization writer
  Teuchos::RCP<IO::DiscretizationWriter>  output_;
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
  bool incremental_;

  /// flag for Eulerian or ALE formulation of equation(s)
  bool isale_;

  /// flag for type of scalar transport problem
  INPAR::SCATRA::ScaTraType scatratype_;

  /// boolian if scatra object is used for level set reinitialization
  bool reinitswitch_;

  /// boolian if a reinitialization step is accepted or not
  bool reinitialization_accepted_;

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

  /// ids of scalars for which flux vectors are written (starting with 1)
  std::vector<int> writefluxids_;

  /// flag for printing out mean values of scalars
  const bool outmean_;

  /// boolean to write Gmsh postprocessing files (input parameter)
  const bool outputgmsh_;

  /// time step size (for reinitialization methods this is a pseudo time step factor)
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

  /// flag for potential Neumann inflow boundary condition
  bool neumanninflow_;

  /// flag for potential boundary condition due to convective heat transfer
  bool convheatrans_;

  /// flag for potentially skiping computation of initial time derivative
  bool skipinitder_;

  /// fine-scale subgrid-viscosity flag
  INPAR::SCATRA::FSSUGRDIFF fssgd_;

  /// ELCH-specific parameter F/RT (default value for the temperature is 298K)
  double frt_;

  /// LOMA-specific parameter: turbulence model
  enum INPAR::FLUID::TurbModelAction turbmodel_;

  /// LOMA-specific parameter: thermodynamic pressure at n+1 and n
  double thermpressnp_;
  double thermpressn_;

  /// LOMA-specific parameter: time deriv. of thermodynamic pressure at n+1 and n
  double thermpressdtnp_;
  double thermpressdtn_;

  /// LOMA-specific parameter: initial mass in domain
  double initialmass_;

  /// system matrix (either sparse matrix or block sparse matrix)
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector>    zeros_;

  /// maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  /// the vector containing body and surface forces
  Teuchos::RCP<Epetra_Vector>    neumann_loads_;

  /// The residual vector
  Teuchos::RCP<Epetra_Vector>    residual_;

  /// phi at times n+1 and n
  Teuchos::RCP<Epetra_Vector>    phinp_;
  Teuchos::RCP<Epetra_Vector>    phin_;

  /// phi at time 0 as reference for reinitialization procedure
  Teuchos::RCP<Epetra_Vector>    phistart_;

  /// time derivative of phi at times n+1 and n
  Teuchos::RCP<Epetra_Vector>  phidtnp_;
  Teuchos::RCP<Epetra_Vector>  phidtn_;

  /// true (rescaled) residual vector without zeros at Dirichlet conditions
  Teuchos::RCP<Epetra_Vector>    trueresidual_;

  /// density (nodebased!) at times n+1, n and n-1 for ELCH natural convection algorithm
  Teuchos::RCP<Epetra_Vector>    elchdensnp_;
  Teuchos::RCP<Epetra_Vector>    elchdensn_;
  Teuchos::RCP<Epetra_Vector>    elchdensnm_;

  /// histvector --- a linear combination of phinm, phin (BDF)
  ///                or phin, phidtn (One-Step-Theta)
  Teuchos::RCP<Epetra_Vector>    hist_;

  /// fluid velocity
  Teuchos::RCP<Epetra_MultiVector> vel_;

  /// convective velocity
  Teuchos::RCP<Epetra_MultiVector>    convel_;

  /// subgrid-scale velocity required for multifractal subgrid-scale modeling
  Teuchos::RCP<Epetra_MultiVector>    fsvel_;

  /// acceleration and pressure required for computation of subgrid-scale velocity
  Teuchos::RCP<Epetra_MultiVector>    accpre_;

  /// displacements at time n+1 for ALE formulation
  Teuchos::RCP<Epetra_MultiVector>    dispnp_;

  /// Nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector>    increment_;

  /// subgrid-diffusivity(-scaling) vector
  Teuchos::RCP<Epetra_Vector>    subgrdiff_;

  /// only necessary for AVM3: fine-scale subgrid-diffusivity matrix
  Teuchos::RCP<LINALG::SparseMatrix> sysmat_sd_;

  /// only necessary for AVM3: scale-separation matrix
  Teuchos::RCP<LINALG::SparseMatrix> Sep_;

  /// only necessary for AVM3: normalized fine-scale subgrid-viscosity matrix
  Teuchos::RCP<LINALG::SparseMatrix> Mnsv_;

  /// flux vector fields
  Teuchos::RCP<Epetra_MultiVector> flux_;

  /// Extractor used for convergence check either in ELCH or LOMA case
  Teuchos::RCP<LINALG::MapExtractor> splitter_;

  /// Periodic Boundary Conditions. Guess what.
  Teuchos::RCP<PeriodicBoundaryConditions> pbc_;

  /// connection between master and slave row nodes
  Teuchos::RCP<std::map<int,std::vector<int> > > pbcmapmastertoslave_;

  /// time measurement
  double dtele_;
  double dtsolve_;

  /// error file handle
  FILE* errfile_;

  /// vector for statistical evaluation of normal fluxes
  Teuchos::RCP<Epetra_SerialDenseVector> sumnormfluxintegral_;

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

  /// magnetic field (3D vector field)
  Teuchos::RCP<Epetra_MultiVector> magneticfield_;

  /// conductivity of the single species and the whole electrolyte
  Epetra_SerialDenseVector sigma_;

  /// unit outer normal vector field for flux output
  Teuchos::RCP<Epetra_MultiVector> normals_;

  /// densification coefficients for ELCH natural convection
  std::vector<double> densific_;

  /// initial concentrations of the electrolyte solution (ELCH)
  std::vector<double> c0_;

  /// flag for Krylov space projection
  bool                             project_;
  /// integrated basis vector (for each undetermined component)
  Teuchos::RCP<Epetra_MultiVector> w_;
  /// basis vector for nullspace (vector of ones for each undetermined component)
  Teuchos::RCP<Epetra_MultiVector> c_;

  /// coupling of fluid-fluid at an internal interface
  Teuchos::RCP<FLD::Meshtying>   meshtying_;

  // options for meshtying
  int msht_;

  // flag for turbulent inflow (turbulent loma specific)
  bool turbinflow_;
  // number of inflow generation time steps (turbulent loma specific)
  int numinflowsteps_;

}; // class ScaTraTimIntImpl

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_IMPLICIT_H
#endif  // #ifdef CCADISCRET
