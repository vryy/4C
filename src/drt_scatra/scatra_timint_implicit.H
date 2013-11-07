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
Maintainer: Andreas Ehrl
            ehrl@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15252
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_IMPLICIT_H
#define SCATRA_TIMINT_IMPLICIT_H

#include <Teuchos_RCP.hpp>
#include <Teuchos_ParameterList.hpp>
#include <Epetra_Vector.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_CrsGraph.h>
#include "../drt_lib/drt_dserror.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_inpar/inpar_fluid.H"

/*==========================================================================*/
// Style guide                                                    nis Mar12
/*==========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  /*--- calculate and update -----------------------------------------------*/

  /*--- query and output ---------------------------------------------------*/



/*==========================================================================*/
// forward declarations
/*==========================================================================*/

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
  class KrylovProjector;
}

namespace FLD
{
  class Meshtying;
  class DynSmagFilter;
  class Vreman;
}


/*==========================================================================*/
//! Algorithm and time integration classes for scalar transport problems
/*==========================================================================*/
namespace SCATRA
{

class HomIsoTurbScalarForcing;

/*!
 * \brief implicit time integration for scalar transport problems
 *
 * \author Georg Bauer (bauer@lnm.mw.tum.de)
 */

class ScaTraTimIntImpl
{
  friend class HomIsoTurbInitialScalarField;
  friend class HomIsoTurbScalarForcing;

public:

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  //! Standard Constructor
  ScaTraTimIntImpl(Teuchos::RCP<DRT::Discretization>       dis,
                   Teuchos::RCP<LINALG::Solver>            solver,
                   Teuchos::RCP<Teuchos::ParameterList>    params,
                   Teuchos::RCP<Teuchos::ParameterList>    extraparams,
                   Teuchos::RCP<IO::DiscretizationWriter>  output);

  //! Destructor
  virtual ~ScaTraTimIntImpl();

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! setup the variables to do a new time step
  void PrepareTimeStep();

  // TODO: Ursula: after reorganization of level-set handling, move to protected again
  //! initialization procedure before the first time step is done
  virtual void PrepareFirstTimeStep() = 0;

  //! preparations for solve
  void PrepareLinearSolve();

  //! set time and step value
  void SetTimeStep(
    double time, //!< new time value
    int step     //!< new step value
    )
  {
    time_ = time;
    step_ = step;
  };

  //! do explicit predictor step (-> better starting value for nonlinear solver)
  virtual void ExplicitPredictor() = 0;

  //! set the velocity field (zero or field by function)
  void SetVelocityField();

  //! set convective velocity field (+ pressure and acceleration field as
  //! well as fine-scale velocity field, if required)
  void SetVelocityField(
    Teuchos::RCP<const Epetra_Vector> convvel, //!< convective velocity/press. vector
    Teuchos::RCP<const Epetra_Vector> acc,     //!< acceleration vector
    Teuchos::RCP<const Epetra_Vector> vel,     //!< velocity vector
    Teuchos::RCP<const Epetra_Vector> fsvel,   //!< fine-scale velocity vector
    Teuchos::RCP<const DRT::DofSet>   dofset,  //!< dofset pointer
    Teuchos::RCP<DRT::Discretization> dis);    //!< fluid/structure discretization

  //! redistribute the scatra discretization and vectors according to nodegraph
  virtual void Redistribute(
    const Teuchos::RCP<Epetra_CrsGraph> nodegraph //!< new nodal distribution over procs
  );

  //! read restart data
  virtual void ReadRestart(int step) = 0;

  /*--- calculate and update -----------------------------------------------*/

  //! do time integration (time loop)
  void TimeLoop();

  //! general solver call for coupled algorithms (decides if linear/nonlinear internally)
  void Solve();

  //! update the solution after convergence of the nonlinear iteration.
  virtual void Update(
      const int            num=0 //!< field number
  ) = 0;

  //! apply moving mesh data
  void ApplyMeshMovement(
    Teuchos::RCP<const Epetra_Vector> dispnp, //!< displacement vector
    Teuchos::RCP<DRT::Discretization> dis     //!< fluid/structure discretization
  );

  //! calculate mass/heat flux vectors (in domain / at boundary)
  Teuchos::RCP<Epetra_MultiVector> CalcFlux(
    const bool writetofile,   //!< flag for writing flux info to file
    const int num=0           //!< field number
  );

  //! calculate mass/heat flux vector field in computational domain
  Teuchos::RCP<Epetra_MultiVector> CalcFluxInDomain(
    const INPAR::SCATRA::FluxType fluxtype //!< fluxtype enum
  );

  //! calculate mass/heat normal flux at specified boundaries
  Teuchos::RCP<Epetra_MultiVector> CalcFluxAtBoundary(
    std::vector<std::string>& condnames,   //!< ?
    const bool           writetofile, //!< ?
    const int            num=0,       //!< field number
    bool                 biogrowth=0
  );

  //! calculate error compared to analytical solution
  void EvaluateErrorComparedToAnalyticalSol();

  /*--- query and output ---------------------------------------------------*/

  //! Provide enum of time integration scheme
  enum INPAR::SCATRA::TimeIntegrationScheme MethodName(){ return timealgo_; }

  //! Provide title of time integration scheme
  std::string MethodTitle(){ return MapTimIntEnumToString(MethodName()); }

  //! return type of physics we are solving for
  INPAR::SCATRA::ScaTraType ScaTraType() const { return scatratype_; };

  //! Print information about current time step to screen
  virtual void PrintTimeStepInfo();

  //! return system matrix down-casted as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

  //! return system matrix downcasted as block sparse matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

  //! return residual
  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; };

  //! return trueresidual
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }

  //! return flag indicating if an incremental solution approach is used
  bool IsIncremental() { return incremental_; }

  //! return scalar field convective velocity at time step n+1
  Teuchos::RCP<Epetra_MultiVector> ConVel() { return convel_; }

  //! return scalar field acceleration/history of convective velocity field at time step n+1
  Teuchos::RCP<Epetra_MultiVector> ConAccPre() { return accpre_; }

  //! return scalar field fine-scale convective velocity at time step n+1
  Teuchos::RCP<Epetra_MultiVector> ConFsVel() { return fsvel_; }

  //! return flux vector
  Teuchos::RCP<Epetra_MultiVector> Flux() { return flux_; };

  //! return Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  //! return information on periodic bondary conditions (row map)
  Teuchos::RCP<std::map<int,std::vector<int> > > PBCmap() { return pbcmapmastertoslave_; };

  //! return the periodic boundary condition member
  Teuchos::RCP<PeriodicBoundaryConditions> PBC() { return pbc_; }

  //! return discretization
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

  //! return the parameter lists
  Teuchos::RCP<Teuchos::ParameterList> ScatraParameterList() { return params_;}
  Teuchos::RCP<Teuchos::ParameterList> ScatraExtraParameterList() { return extraparams_;}
  virtual Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList() = 0;

  //! Access output object: CD-Rom and DVD only - no BlueRay support!!! ;)
  IO::DiscretizationWriter& DiscWriter() { return *output_; }

  //! returns map extractor used for convergence check either in ELCH or LOMA case
  Teuchos::RCP<LINALG::MapExtractor> Splitter() { return splitter_; }

  //! output solution and restart data to file
  void Output(const int num =0);

  /*========================================================================*/
  //! @name Time, time-step and related methods
  /*========================================================================*/

  /*--- query and output ---------------------------------------------------*/

  //! return current time value
  double Time() const { return time_; }

  //! return current step number
  int Step() const { return step_; }

  //! total number of time steps ? rename StepMax?
  int NStep() const { return stepmax_;}

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! set the initial scalar field phi
  void SetInitialField(
    const INPAR::SCATRA::InitialField  init,       //!< type of initial field
    const int                          startfuncno //!< number of spatial function
  );

  /*--- calculate and update -----------------------------------------------*/

  //! calculate initial value for phidtn_ (only necessary for OnestepTheta and GenAlpha)
  void CalcInitialPhidt();

  //! calculate initial value for phidtn_ (assembly of system)
  void CalcInitialPhidtAssemble();

  //! calculate initial value for phidtn_ (solver part)
  void CalcInitialPhidtSolve();

  //! compute time derivative
  virtual void ComputeTimeDerivative() = 0;

  //! compute parameters of the Input voltage to use for the double layer current density
  virtual void ComputeTimeDerivPot0(const bool init) = 0;

  //! compute values at intermediate time steps (required for
  //generalized-alpha) ? rename?
  virtual void ComputeIntermediateValues() = 0;

  //! iterative update of phinp
  void UpdateIter(
    const Teuchos::RCP<const Epetra_Vector> inc //!< increment vector for phi
  );

  /*--- query and output ---------------------------------------------------*/

  //! return number of transported scalars
  int NumScal() const { return numscal_; }

  // return scalar field phi at time step n-1
  virtual Teuchos::RCP<Epetra_Vector> Phinm() = 0;

  //! return scalar field phi at time step n
  Teuchos::RCP<Epetra_Vector> Phin() { return phin_; }

  //! return scalar field phi at time step n+1
  Teuchos::RCP<Epetra_Vector> Phinp() { return phinp_; }

  //! return time derivative of scalar field phi at time step n
  Teuchos::RCP<Epetra_Vector> Phidtn() { return phidtn_; }

  //! return time derivative of scalar field phi at time step n+1
  Teuchos::RCP<Epetra_Vector> Phidtnp() { return phidtnp_; }

  //! return scalar field history
  Teuchos::RCP<Epetra_Vector> Hist() { return hist_; }

  //! return scalar field phi at time step n+alpha_F
  virtual Teuchos::RCP<Epetra_Vector> Phiaf() = 0;

  //! return scalar field phi at time step n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> Phiam() = 0;

  //! return time derivative of scalar field phi at time step n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> Phidtam() = 0;

  //! return fine-scale scalar field fsphi at time step n+1 or alpha_M
  virtual Teuchos::RCP<Epetra_Vector> FsPhi() = 0;

  //! output of mean values of scalar(s)
  void OutputMeanScalars(const int num=0);

  //! output of reaction(s) integral
  void OutputIntegrReac(const int num=0);

  /*========================================================================*/
  //! @name ELCH
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  /*--- calculate and update -----------------------------------------------*/

  //! compute density field from concentration(s) for ELCH natural convection
  void ComputeDensity();

  //! update density field for ELCH natural convection
  virtual void UpdateDensityElch() = 0;

  /*--- query and output ---------------------------------------------------*/

  //! return density field at time n+1 for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchNp() { return elchdensnp_; }

  //! return density field at time n for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchN() { return elchdensn_; }

  //! return density field at time n-1 for ELCH natural convection
  Teuchos::RCP<Epetra_Vector> DensElchNm() { return elchdensnm_; }

  //! output of all electrode information to screen and/or file (default: 2x true)
  Teuchos::RCP< std::vector<double> > OutputElectrodeInfo(
    bool printtoscreen = true,  //!< flag for printing information to screen
    bool printtofile = true     //!< flag for printing information to file
  );

  //! output of electrode information for a single boundary condition
  Teuchos::RCP< std::vector<double> > OutputSingleElectrodeInfo(
    DRT::Condition* condition, //!< pointer to the condition
    const int  condid,         //!< condition id
    const bool printtoscreen,  //!< flag for output to screen
    const bool printtofile,    //!< flag for output to file
    double& currentsum,        //!< net current over all boundaries (out)
    double& pottangent,        //!< tangent of current w.r.t. electrode potential (out)
    double& currresidual,      //!< neg. residual of current equation (out)
    double& electrodesurface,  //!< surface of the electrode
    double& electrodepot,      //!< electrode potential (metal side)
    double& meanoverpot        //!< mean surface overpotential
  );

  /*========================================================================*/
  //! @name low-Mach-number flow
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! set initial thermodynamic pressure for low-Mach-number flow
  void SetInitialThermPressure();

  //! predict thermodynamic pressure and time derivative for low-Mach-number flow
  virtual void PredictThermPressure() = 0;

  /*--- calculate and update -----------------------------------------------*/

  //! compute initial thermodyn. pressure time derivative for low-Mach-number flow
  void ComputeInitialThermPressureDeriv();

  //! compute initial total mass in domain for low-Mach-number flow
  void ComputeInitialMass();

  //! compute thermodynamic pressure and time derivative for low-Mach-number flow
  virtual void ComputeThermPressure() = 0;

  //! compute thermodyn. press. from mass cons. in domain for low-Mach-number flow
  void ComputeThermPressureFromMassCons();

  //! compute values of thermodynamic pressure at intermediate time steps
  //! (required for generalized-alpha)
  virtual void ComputeThermPressureIntermediateValues() = 0;

  //!  compute time derivative of thermodynamic pressure after solution
  virtual void ComputeThermPressureTimeDerivative() = 0;

  //! update thermodynamic pressure and time derivative for low-Mach-number flow
  virtual void UpdateThermPressure() = 0;

  //! check convergence (only for low-Mach-number flow, for the time being)
  bool ConvergenceCheck(
    int          itnum,   //!< ?
    int          itmax,   //!< ?
    const double ittol    //!< ?
  );

  /*--- query and output ---------------------------------------------------*/

  //! return thermo. press. at time step n for low-Mach-number flow
  double ThermPressN() const { return thermpressn_; }

  //! return thermo. press. at time step n+1 for low-Mach-number flow
  double ThermPressNp() const { return thermpressnp_; }

  //! return thermo. press. at time step n+alpha_F for low-Mach-number flow
  virtual double ThermPressAf() = 0;

  //! return thermo. press. at time step n+alpha_M for low-Mach-number flow
  virtual double ThermPressAm() = 0;

  //! return time der. of thermo. press. at time step n+1 for low-Mach-number flow
  double ThermPressDtNp() const { return thermpressdtnp_; }

  //! return time derivative of thermo. press. at time step n+alpha_F for low-Mach-number flow
  virtual double ThermPressDtAf() = 0;

  //! return time derivative of thermo. press. at time step n+alpha_M for low-Mach-number flow
  virtual double ThermPressDtAm() = 0;

  /*========================================================================*/
  //! @name functions used for reinitialization of level sets
  /*========================================================================*/

  //! update the solution after reinitialization
  virtual void UpdateReinit() = 0;

  /// routine to return scalar field phi at time step n
  Teuchos::RCP<Epetra_Vector> PhiReinitStart() { return phistart_; }

  /*========================================================================*/
  //! @name turbulence and related
  /*========================================================================*/

  ///! get access to dynamic Smagorinsky class of fluid time integration
  void AccessDynSmagFilter(Teuchos::RCP<FLD::DynSmagFilter> dynSmag);
  ///! get access to dynamic Vreman class of fluid time integration
  void AccessVreman(Teuchos::RCP<FLD::Vreman> vrem);

  ///! calculate intermediate solution to determine forcing for homogeneous isotropic turbulence
  void CalcIntermediateSolution();

  /*========================================================================*/
  //! @name Biofilm methods
  /*========================================================================*/

  //! return scatra structure growth vector
  Teuchos::RCP<Epetra_MultiVector> StrGrowth() { return scstrgrdisp_; };

  //! return scatra fluid growth vector
  Teuchos::RCP<Epetra_MultiVector> FldGrowth() { return scfldgrdisp_; };

  //! set scatra fluid displacement vector due to biofilm growth
  void SetScFldGrDisp(Teuchos::RCP<Epetra_MultiVector> scatra_fluid_growth_disp);

  //! set scatra structure displacement vector due to biofilm growth
  void SetScStrGrDisp(Teuchos::RCP<Epetra_MultiVector> scatra_struct_growth_disp);

  /*========================================================================*/
  //! @name  not classified method - to be kept clean!!!
  /*========================================================================*/

  //! compute contribution of permeable surface/interface
  void SurfacePermeability(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< ?
    Teuchos::RCP<Epetra_Vector>          rhs     //!< ?
  );

  /*========================================================================*/
  //! @name  obsolete or unused methods - to be deleted soon (at noon)!!!
  /*========================================================================*/

  // nothing here at the moment

  /*========================================================================*/

protected:

  /*========================================================================*/
  //! @name Constructors and destructors and related methods
  /*========================================================================*/

  //! don't want = operator
  //ScaTraTimIntImpl operator = (const ScaTraTimIntImpl& old);

  //! don't want copy constructor
  ScaTraTimIntImpl (const ScaTraTimIntImpl& old);

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! Set the part of the righthandside belonging to the last timestep.
  virtual void SetOldPartOfRighthandside() = 0;

  //! set time for evaluation of Neumann boundary conditions
  virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) = 0;

  //! create Krylov space projector
  void SetupKrylovSpaceProjection(DRT::Condition* kspcond);
  //! update Krylov space projector
  void UpdateKrylovSpaceProjection();

  //! export node-based multivector to column map and add it to a parameter list
  void AddMultiVectorToParameterList(
    Teuchos::ParameterList&          p,    //!< parameter list
    const std::string                name, //!< naming std::string for added vector
    Teuchos::RCP<Epetra_MultiVector> vec   //!< Epetra_MultiVector to be added
  );

  //! compute approximation for fluxes and add it to a parameter list
  void AddFluxApproxToParameterList(
    Teuchos::ParameterList&            p,       //!< parameter list
    const enum INPAR::SCATRA::FluxType fluxtype //!< switch for flux calculation
  );

  //! add parameters specific for time-integration scheme
  virtual void AddSpecificTimeIntegrationParameters(Teuchos::ParameterList& params) = 0;

  /*--- calculate and update -----------------------------------------------*/

  //! apply Dirichlet boundary conditions to linear system of equations
  void ApplyDirichletToSystem();

  //! Apply Dirichlet boundary conditions on provided state vector
  void ApplyDirichletBC(
    const double                time,  //!< evaluation time
    Teuchos::RCP<Epetra_Vector> phinp, //!< transported scalar(s) (may be = null)
    Teuchos::RCP<Epetra_Vector> phidt  //!< first time derivative (may be = null)
  );

  //! compute outward pointing unit normal vectors at given bc's
  Teuchos::RCP<Epetra_MultiVector> ComputeNormalVectors(
    const std::vector<std::string>& condnames //!< ?
  );

  //! evaluate Neumann inflow boundary condition
  void ComputeNeumannInflow(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< ?
    Teuchos::RCP<Epetra_Vector>          rhs     //!< ?
  );

  //! evaluate boundary condition due to convective heat transfer
  void EvaluateConvectiveHeatTransfer(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< ?
    Teuchos::RCP<Epetra_Vector>          rhs     //!< ?
  );

  //! potential residual scaling and potential addition of Neumann terms
  void ScalingAndNeumann();

  //! add actual Neumann loads multipl. with time factor to the residual
  virtual void AddNeumannToResidual() = 0;

  //! Apply Neumann boundary conditions
  void ApplyNeumannBC(
    const double                      time,         //!< evaluation time
    const Teuchos::RCP<Epetra_Vector> phinp,        //!< transported scalar(s) state
    Teuchos::RCP<Epetra_Vector>       neumann_loads //!< Neumann loads
  );

  //! call elements to calculate system matrix and rhs and assemble
  void AssembleMatAndRHS();

  //! return the right time-scaling-factor for the true residual
  virtual double ResidualScaling() const = 0;

  //! solve linear system
  void LinearSolve();

  //! contains the nonlinear iteration loop
  void NonlinearSolve();

  //! check convergence (or divergence) of nonlinear iteration
  bool AbortNonlinIter(
    const int itnum,        //!< current value of iteration step counter
    const int itmax,        //!< maximum number of iteration steps
    const double ittol,     //!< relative tolerance for increments
    const double abstolres, //!< absolute tolerance for the residual norm
    double& actresidual     //!< return value of the current residual
  );

  /*--- query and output ---------------------------------------------------*/

  //! helper function to get algorithm title
  std::string MapTimIntEnumToString
  (
    const enum INPAR::SCATRA::TimeIntegrationScheme term  //!< the enum
  );

  //! do we need a statistical sampling for boundary flux at the current time step?
  bool DoBoundaryFluxStatistics()
    { return( (step_>=samstart_) and (step_<=samstop_)
              and
              ( (writeflux_==INPAR::SCATRA::flux_total_boundary)
                or (writeflux_==INPAR::SCATRA::flux_diffusive_boundary)
                or (writeflux_==INPAR::SCATRA::flux_convective_boundary))); };

  //! is output needed for the current time step?
  bool DoOutput()
    { return ((step_%upres_==0 )or (step_%uprestart_==0)); };

  //! write state vectors (phinp and convective velocity) to BINIO
  void OutputState();

  //! write state vectors (phinp and convective velocity) to Gmsh postprocessing files
  void OutputToGmsh(
      const int step,
      const double time
  ) const;

  //! write additional data required for restart
  virtual void OutputRestart() = 0;

  //! write flux vectors to BINIO
  void OutputFlux(
    Teuchos::RCP<Epetra_MultiVector> flux //!< flux vector
  );

  /*========================================================================*/
  //! @name Time, time-step and related methods
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! increment time and step value
  void IncrementTimeAndStep();

  /*--- calculate and update -----------------------------------------------*/

  //! time update of time-dependent materials
  void ElementMaterialTimeUpdate();

  /*--- query and output ---------------------------------------------------*/

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  /*--- calculate and update -----------------------------------------------*/

  /*--- query and output ---------------------------------------------------*/

  /*========================================================================*/
  //! @name ELCH
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! setup of natural convection (ELCH)
  void SetupElchNatConv();

  //! set magnetic field (ELCH)
  void SetMagneticField(
    const int funcno //!< number of spatial function that defines the vector field
  );

  /*--- calculate and update -----------------------------------------------*/

  //! calculate initial field for electric potential based on ion concentrations (ELCH)
  void CalcInitialPotentialField();

  //! computes different conductivity expressions for electrolyte solutions (ELCH)
  Epetra_SerialDenseVector ComputeConductivity();

  //! apply galvanostatic control (update electrode potential) (ELCH)
  bool ApplyGalvanostaticControl();

  //! compute contribution of electrode kinetics to eq. system
  void EvaluateElectrodeKinetics(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< system matrix
    Teuchos::RCP<Epetra_Vector>          rhs     //!< rhs vector
  );

  //! Initialize Nernst-BC
  void InitNernstBC();

  //! Add Linearization for Nernst-BC
  void  LinearizationNernstCondition();

  //! update time-dependent electrode state variables at the end of an time step
  virtual void ElectrodeKineticsTimeUpdate() = 0;

  /*--- query and output ---------------------------------------------------*/

  //! check for negative values of concentrations (ELCH)
  void CheckConcentrationValues(
    Teuchos::RCP<Epetra_Vector> vec //!< current phi vector to be checked
  );

  /*========================================================================*/
  //! @name low-Mach-number flow
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  /*--- calculate and update -----------------------------------------------*/

  /*--- query and output ---------------------------------------------------*/

  /*========================================================================*/
  //! @name AVM3 and related
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! prepare AVM3-based scale separation
  void AVM3Preparation();

  //! AVM3-based scale separation
  virtual void AVM3Separation()=0;

  /*--- calculate and update -----------------------------------------------*/

  //! scaling of AVM3-based subgrid-diffusivity matrix
  void AVM3Scaling(
      Teuchos::ParameterList& eleparams //!< parameter list
  );

  /*--- query and output ---------------------------------------------------*/

  /*========================================================================*/
  //! @name turbulence and related
  /*========================================================================*/

  //! dynamic Smagorinsk model
  virtual void DynamicComputationOfCs()=0;
  
  //! dynamic Vreman model
  virtual void DynamicComputationOfCv()=0;

  //! calculate mean CsgsB to estimate CsgsD for multifractal subgrid-scale model
  void RecomputeMeanCsgsB();

  /*========================================================================*/
  //! @name elch: diffusion-conduction formulation related
  /*========================================================================*/

  void ValidParameterDiffCond();

  /*========================================================================*/
  //! @name functions used for reinitialization of level sets
  /*========================================================================*/


  /*--- calculate and update -----------------------------------------------*/

  //! add reinitialization specific parameters
  void AddReinitializationParameters(
    Teuchos::ParameterList& params //!< parameter list
  );

  //! call elements to calculate system matrix and rhs at boundaries for Taylor-Galerkin methods
  void AssembleMatAndRHS_Boundary();

  //! print information about current pseudo time step to screen
  virtual void PrintPseudoTimeStepInfoReinit()
  {
    if (myrank_==0)
      printf("TIMESTEPFACTOR: %11.4E      %s          PSEUDOSTEP = %4d/%4d \n",
          dta_,MethodTitle().c_str(),step_,stepmax_);
  }

  /*========================================================================*/
  //! @name  not classified method - to be kept clean!!!
  /*========================================================================*/

  //! scale lines of linear system prior to solver call (ELCH)
  void ScaleLinearSystem();

  /*!
   * \brief Extract the Dirichlet toggle vector based on Dirichlet BC maps
   *
   * This method provides backward compatability only. Formerly, the Dirichlet
   * conditions were handled with the Dirichlet toggle vector. Now, they are
   * stored and applied with maps, ie #dbcmaps_. Eventually, this method will
   * be removed.
   * note: VM3 solver still needs an explicit toggle vector for construction
   */
  const Teuchos::RCP<const Epetra_Vector> DirichletToggle();

  /*========================================================================*/

  /*========================================================================*/
  //! @name general framework variables
  /*========================================================================*/

  //! the solver
  Teuchos::RCP<LINALG::Solver>         solver_;

  //! the parameter list
  Teuchos::RCP<Teuchos::ParameterList> params_;

  //! the parameter list containing extra parameters (application dependent)
  Teuchos::RCP<Teuchos::ParameterList> extraparams_;

  //! processor id
  int myrank_;

  //! Extractor used for convergence check either in ELCH or LOMA case
  Teuchos::RCP<LINALG::MapExtractor> splitter_;

  //! error file handle
  FILE* errfile_;

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/

  //! flag for type of scalar transport problem
  INPAR::SCATRA::ScaTraType scatratype_;

  //! flag for Eulerian or ALE formulation of equation(s)
  bool isale_;

  //! solvertype and flags for nonlinear (always incremental) and (linear) incremental solver
  INPAR::SCATRA::SolverType solvtype_;
  //! incremental or linear full solving ? rename -> is_incremental_
  bool incremental_;

  //! flag indicating if initial velocity field has been set
  bool initialvelset_;

  //! flag for fine-scale subgrid-viscosity
  INPAR::SCATRA::FSSUGRDIFF fssgd_;

  //! LOMA-specific parameter: turbulence model
  enum INPAR::FLUID::TurbModelAction turbmodel_;

  /*--- query and output ---------------------------------------------------*/

  //! flag for writing the flux vector fields
  INPAR::SCATRA::FluxType writeflux_;

  //! ids of scalars for which flux vectors are written (starting with 1)
  std::vector<int> writefluxids_;

  //! flag for printing out mean values of scalars
  const bool outmean_;

  //! boolean to write Gmsh postprocessing files (input parameter)
  const bool outputgmsh_;

  /*========================================================================*/
  //! @name Time, time-step, and iteration variables
  /*========================================================================*/

  //! actual time
  double time_;

  //! maximum simulation time
  double maxtime_;

  //! actual step number
  int step_;

  //! maximum number of steps ? name maxtime vs. stepmax
  int stepmax_;

  //! time step size (for reinitialization methods this is a pseudo time step length)
  double dta_;

  //! time measurement element
  double dtele_;
  //! time measurement solve
  double dtsolve_;

  //! used time integration scheme
  INPAR::SCATRA::TimeIntegrationScheme timealgo_;

  /*========================================================================*/
  //! @name scalar degrees of freedom variables
  /*========================================================================*/

  //! number of transported scalars
  int numscal_;

  //! phi at times n
  Teuchos::RCP<Epetra_Vector> phin_;
  //! phi at times n+1
  Teuchos::RCP<Epetra_Vector> phinp_;

  //! time derivative of phi at times n
  Teuchos::RCP<Epetra_Vector> phidtn_;
  //! time derivative of phi at times n+1
  Teuchos::RCP<Epetra_Vector> phidtnp_;

  //! histvector --- a linear combination of phinm, phin (BDF)
  //!                or phin, phidtn (One-Step-Theta)
  Teuchos::RCP<Epetra_Vector>    hist_;

  //! initial concentrations of the electrolyte solution (ELCH)
  std::vector<double> c0_;

  //! phi at time 0 as reference for reinitialization procedure
  Teuchos::RCP<Epetra_Vector> phistart_;


  /*========================================================================*/
  //! @name electrophysiology
  /*========================================================================*/

  //! activation_time at times n+1
  Teuchos::RCP<Epetra_Vector> activation_time_np_;

  //! activation threshold for postprocessing
  double activation_threshold_;

  //! maximum expected number of material internal state variables
  int nb_max_mat_int_state_vars_;

  //! material internal state at times n+1
  Teuchos::RCP<Epetra_MultiVector> material_internal_state_np_;

  //! one component of the material internal state at times n+1 (for separated postprocessing)
  Teuchos::RCP<Epetra_Vector> material_internal_state_np_component_;

  //! maximum expected number of material ionic currents
  int nb_max_mat_ionic_currents_;

  //! material ionic currents at times n+1
  Teuchos::RCP<Epetra_MultiVector> material_ionic_currents_np_;

  //! one component of the material ionic currents at times n+1 (for separated postprocessing)
  Teuchos::RCP<Epetra_Vector> material_ionic_currents_np_component_;

  /*========================================================================*/
  //! @name velocity, pressure, and related
  /*========================================================================*/

  //! fluid velocity
  Teuchos::RCP<Epetra_MultiVector> vel_;

  //! convective velocity
  Teuchos::RCP<Epetra_MultiVector> convel_;

  //! subgrid-scale velocity required for multifractal subgrid-scale modeling
  Teuchos::RCP<Epetra_MultiVector> fsvel_;

  //! acceleration and pressure required for computation of subgrid-scale velocity
  Teuchos::RCP<Epetra_MultiVector> accpre_;

  //! displacements at time n+1 for ALE formulation
  Teuchos::RCP<Epetra_MultiVector> dispnp_;

  //! type of velocity field ? rename!!!
  const INPAR::SCATRA::VelocityField cdvel_;

  /*========================================================================*/
  //! @name coefficients and related
  /*========================================================================*/

  //! density (nodebased!) at times n-1 for ELCH natural convection algorithm
  Teuchos::RCP<Epetra_Vector> elchdensnm_;
  //! density (nodebased!) at times n for ELCH natural convection algorithm
  Teuchos::RCP<Epetra_Vector> elchdensn_;
  //! density (nodebased!) at times n+1 for ELCH natural convection algorithm
  Teuchos::RCP<Epetra_Vector> elchdensnp_;

  //! subgrid-diffusivity(-scaling) vector
  Teuchos::RCP<Epetra_Vector> subgrdiff_;

  /*========================================================================*/
  //! @name Galerkin discretization, boundary conditions, and related
  /*========================================================================*/

  //! the scalar transport discretization
  Teuchos::RCP<DRT::Discretization>    discret_;

  //! the discretization writer
  Teuchos::RCP<IO::DiscretizationWriter> output_;

  //! form of convective term
  INPAR::SCATRA::ConvForm convform_;

  //! system matrix (either sparse matrix or block sparse matrix)
  Teuchos::RCP<LINALG::SparseOperator> sysmat_;

  //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
  Teuchos::RCP<Epetra_Vector> zeros_;

  //! maps for extracting Dirichlet and free DOF sets
  Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

  //! the vector containing body and surface forces
  Teuchos::RCP<Epetra_Vector> neumann_loads_;

  //! Periodic Boundary Conditions. Guess what.
  Teuchos::RCP<PeriodicBoundaryConditions> pbc_;

  //! connection between master and slave row nodes
  Teuchos::RCP<std::map<int,std::vector<int> > > pbcmapmastertoslave_;

  //! unit outer normal vector field for flux output
  Teuchos::RCP<Epetra_MultiVector> normals_;

  //! flux vector fields
  Teuchos::RCP<Epetra_MultiVector> flux_;

  //! vector for statistical evaluation of normal fluxes
  Teuchos::RCP<Epetra_SerialDenseVector> sumnormfluxintegral_;

  //! the last step number when fluxes have been computed
  int lastfluxoutputstep_;

  //! The residual vector
  Teuchos::RCP<Epetra_Vector> residual_;

  //! true (rescaled) residual vector without zeros at Dirichlet conditions
  Teuchos::RCP<Epetra_Vector> trueresidual_;

  //! Nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector> increment_;

  //! coupling of fluid-fluid at an internal interface
  Teuchos::RCP<FLD::Meshtying> meshtying_;

  // options for meshtying ? bool?
  int msht_;

  /*========================================================================*/
  //! @name ELCH variables
  /*========================================================================*/

  //! number of iterations in galvanostatic mode (ELCH)
  int gstatnumite_;

  //! value of electric potential increment in galvanostatic mode (ELCH)
  double gstatincrement_;

  //! magnetic field (3D vector field)
  Teuchos::RCP<Epetra_MultiVector> magneticfield_;

  //! conductivity of the single species and the whole electrolyte
  Epetra_SerialDenseVector sigma_;

  //! densification coefficients for ELCH natural convection
  std::vector<double> densific_;

  //! ELCH-specific parameter F/RT (default value for the temperature is 298K)
  double frt_;

  //! true if kinetic model is butler volmer and one electrode has a double layer capacitance unequal zero
  // initialized with true and set false if no capacitance is detected in the first timestep
  bool dlcapexists_;

  //! electro-kinetics toggle
  //! Toggle defining the dof's with Nernst-BC
  Teuchos::RCP<Epetra_Vector> ektoggle_;

  /*========================================================================*/
  //! @name AVM3 variables
  /*========================================================================*/

  //! only necessary for AVM3: fine-scale subgrid-diffusivity matrix
  Teuchos::RCP<LINALG::SparseMatrix> sysmat_sd_;

  //! only necessary for AVM3: scale-separation matrix ? rename small caps
  Teuchos::RCP<LINALG::SparseMatrix> Sep_;

  //! only necessary for AVM3: normalized fine-scale subgrid-viscosity matrix ? rename small caps
  Teuchos::RCP<LINALG::SparseMatrix> Mnsv_;

  /*========================================================================*/
  //! @name low-Mach-number flow variables
  /*========================================================================*/

  //! LOMA-specific parameter: initial mass in domain
  double initialmass_;

  //! LOMA-specific parameter: thermodynamic pressure at n
  double thermpressn_;
  //! LOMA-specific parameter: thermodynamic pressure at n+1
  double thermpressnp_;

  //! LOMA-specific parameter: time deriv. of thermodynamic pressure at n
  double thermpressdtn_;
  //! LOMA-specific parameter: time deriv. of thermodynamic pressure at n+1
  double thermpressdtnp_;

  //! number of inflow generation time steps (turbulent loma specific)
  int numinflowsteps_;

  /*========================================================================*/
  //! @name turbulent flow variables
  /*========================================================================*/

  //! one instance of the filter object
  Teuchos::RCP<FLD::DynSmagFilter> DynSmag_;
  Teuchos::RCP<FLD::Vreman> Vrem_;

  //! parameters for sampling/dumping period
  int samstart_;
  int samstop_;
  int dumperiod_;

  //! flag for turbulent inflow (turbulent loma specific)
  bool turbinflow_;

  /// flag for special turbulent flow
  std::string special_flow_;

  //! the vector containing source term externally computed
  //! forcing for homogeneous isotropic turbulence
  Teuchos::RCP<Epetra_Vector> forcing_;

  //! forcing for homogeneous isotropic turbulence
  Teuchos::RCP<SCATRA::HomIsoTurbScalarForcing> homisoturb_forcing_;

  /*========================================================================*/
  //! @name variables for reinitialization of level sets
  /*========================================================================*/

  //! boolian if scatra object is used for level set reinitialization
  bool reinitswitch_;

  /*========================================================================*/
  //! @name variables for orthogonal space projection aka Krylov projection
  /*========================================================================*/

  bool updateprojection_;                             //!< bool triggering update of Krylov projection
  Teuchos::RCP<LINALG::KrylovProjector>   projector_; //!< Krylov projector himself

  /*========================================================================*/
  //! @name not classified variables - to be kept clean!!!
  /*========================================================================*/

  //! write results every upres_ steps ? writesolutionevery_
  int upres_;

  //! write restart data every uprestart_ steps ? writesolutioneveryrestart_
  int uprestart_;

  /*========================================================================*/
  //! @name obsolete or unused variables - to be deleted soon (at noon)!!!
  /*========================================================================*/

  //! flag for potential Neumann inflow boundary condition
  bool neumanninflow_;

  //! flag for potential boundary condition due to convective heat transfer
  bool convheatrans_;

  //! flag for potentially skiping computation of initial time derivative
  bool skipinitder_;

  /*========================================================================*/

  /*========================================================================*/
  //! @name Biofilm specific stuff
  /*========================================================================*/

  //! scatra fluid displacement due to growth
  Teuchos::RCP<Epetra_MultiVector> scfldgrdisp_;

  //! scatra structure displacement due to growth
  Teuchos::RCP<Epetra_MultiVector> scstrgrdisp_;

  //! flag for printing out integral values of reaction
  const bool outintegrreac_;

  //! Activate/Deactivate diffcond formulation on time integration level
  bool diffcond_;

//  std::vector<double> reacintegral_;


}; // class ScaTraTimIntImpl

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_IMPLICIT_H
