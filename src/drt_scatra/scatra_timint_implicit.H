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
#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseVector.h>

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

  //! initialize time integration
  virtual void Init();

  //! Initialization of system matrix
  virtual void InitSystemMatrix();

  //! setup meshtying system
  virtual void SetupMeshtying();

  //! Initialization of turbulence models
  void InitTurbulenceModel(
      const Epetra_Map* dofrowmap,
      const Epetra_Map* noderowmap);

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! setup the variables to do a new time step
  virtual void PrepareTimeStep();

  //! initialization procedure before the first time step is done
  virtual void PrepareFirstTimeStep() = 0;

  //! preparations for solve
  virtual void PrepareLinearSolve();

  //! set time and step value
  virtual void SetTimeStep(
    const double time, //!< new time value
    const int step     //!< new step value
    )
  {
    time_ = time;
    step_ = step;
  };

  //! do explicit predictor step (-> better starting value for nonlinear solver)
  virtual void ExplicitPredictor() = 0;

  //! set the velocity field (zero or field by function)
  virtual void SetVelocityField();

  //! set convective velocity field (+ pressure and acceleration field as
  //! well as fine-scale velocity field, if required)
  virtual void SetVelocityField(
    Teuchos::RCP<const Epetra_Vector> convvel, //!< convective velocity/press. vector
    Teuchos::RCP<const Epetra_Vector> acc,     //!< acceleration vector
    Teuchos::RCP<const Epetra_Vector> vel,     //!< velocity vector
    Teuchos::RCP<const Epetra_Vector> fsvel,   //!< fine-scale velocity vector
    Teuchos::RCP<const DRT::DofSet>   dofset,  //!< dofset pointer
    Teuchos::RCP<DRT::Discretization> dis);    //!< fluid/structure discretization

  //! read restart data
  virtual void ReadRestart(int step) = 0;

  /*--- calculate and update -----------------------------------------------*/

  //! do time integration (time loop)
  virtual void TimeLoop();

  //! general solver call for coupled algorithms (decides if linear/nonlinear internally)
  virtual void Solve();

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
    bool                 fixtotflux=0 //!< always calculate flux_total_boundary
  );

  //! calculate error compared to analytical solution
  virtual void EvaluateErrorComparedToAnalyticalSol(){return;};

  //! Calculate the reconstructed nodal gradient of phi
  Teuchos::RCP<Epetra_MultiVector> GetGradientAtNodes();

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

  //! return scalar field velocity at time step n+1
  Teuchos::RCP<Epetra_MultiVector> Vel() { return vel_; }

  //! return scalar field displacement at time step n+1
  Teuchos::RCP<Epetra_MultiVector> Disp() { return dispnp_; }

  //! return scalar field acceleration/history of convective velocity field at time step n+1
  Teuchos::RCP<Epetra_MultiVector> ConAccPre() { return accpre_; }

  //! return scalar field fine-scale convective velocity at time step n+1
  Teuchos::RCP<Epetra_MultiVector> ConFsVel() { return fsvel_; }

  //! return flux vector
  Teuchos::RCP<Epetra_MultiVector> Flux() { return flux_; };

  //! return Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

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
  virtual void Output(const int num =0);

  //! print scatra typ
  void PrintScatraType();

  //! Convergence check for two way coupled ScaTra problems.
  bool ConvergenceCheck(int          itnum,
                        int          itmax,
                        const double ittol);

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

  //! return number of newton iterations in last timestep
  double IterNum() const { return iternum_; }


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

  //! return values of scalar field phi at nodes (meshfree)
  Teuchos::RCP<Epetra_Vector> Phiatmeshfreenodes() { return phiatmeshfreenodes_; }

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
  //! @name Elch methods
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! Initialize Nernst-BC
  virtual void InitNernstBC(){return;};

  /*--- calculate and update -----------------------------------------------*/

  //! calculate initial field for electric potential based on ion concentrations (ELCH)
  virtual void CalcInitialPotentialField(){return;};

  //! compute contribution of electrode kinetics to eq. system
  virtual void EvaluateSolutionDependingBC(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< system matrix
    Teuchos::RCP<Epetra_Vector>          rhs     //!< rhs vector
  ){return;};

  //! output of all electrode information to screen and/or file (default: 2x true)
  virtual Teuchos::RCP< std::vector<double> > OutputElectrodeInfo(
    bool printtoscreen = true,  //!< flag for printing information to screen
    bool printtofile = true     //!< flag for printing information to file
  ){return Teuchos::null;};

  //! Add Linearization for Nernst-BC
  virtual void  LinearizationNernstCondition(){return;};

  //! apply galvanostatic control (update electrode potential)
  virtual bool ApplyGalvanostaticControl(){return true;};

  //! check for negative values of concentrations (ELCH)
  virtual void CheckConcentrationValues(
    Teuchos::RCP<Epetra_Vector> vec //!< current phi vector to be checked
  ){return;};

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
  //! @name  fs3i methods
  /*========================================================================*/

  //! compute contribution of permeable surface/interface
  void SurfacePermeability(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< ?
    Teuchos::RCP<Epetra_Vector>          rhs     //!< ?
  );

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

  // nothing here at the moment

  /*========================================================================*/

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

  //! create vectors for Krylov projection if necessary
  void PrepareKrylovProjection();

  /*========================================================================*/
  //! @name set element parameters
  /*========================================================================*/

  virtual void SetElementTimeParameter(bool forcedincrementalsolver=false)=0;

  //! Set backward Euler time parameter
  virtual void SetElementTimeParameterBackwardEuler(){return;};

  //! set time for evaluation of Neumann boundary conditions
  virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) = 0;

  //! Set general element parameter including stabilization
  virtual void SetElementGeneralScaTraParameter();

  //! Set general element parameter with deactivated stabilization
  virtual void SetElementGeneralScaTraParameterDeactivatedStab();

  //! Set element parameter specific for turbulence
  void SetElementTurbulenceParameter();

  //! Deactivate element parameter specific for turbulence
  void SetElementTurbulenceParameterDeactivated();

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! Set the part of the righthandside belonging to the last timestep.
  virtual void SetOldPartOfRighthandside() = 0;

  //! create Krylov space projector
  void SetupKrylovSpaceProjection(DRT::Condition* kspcond);
  //! update Krylov space projector
  void UpdateKrylovSpaceProjection();

  //! compute approximation for fluxes and add it to a parameter list
  void AddFluxApproxToParameterList(
    Teuchos::ParameterList&            p,       //!< parameter list
    const enum INPAR::SCATRA::FluxType fluxtype //!< switch for flux calculation
  );

  //! add parameters specific for time-integration scheme
  virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver=false) = 0;

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

  //! add parameters depending on the problem, i.e., loma, level-set, ...
  virtual void AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params //!< parameter list
  ){return;};

  //! add parameters depending on the problem, i.e., loma, level-set, ... with special modifications for
  //! calculation fo initial phidt
  virtual void AddProblemSpecificParametersAndVectorsForCalcInitialPhiDt(
    Teuchos::ParameterList& params //!< parameter list
  ){return;};

  //! return the right time-scaling-factor for the true residual
  virtual double ResidualScaling() const = 0;

  //! solve linear system
  void LinearSolve();

  //! contains the nonlinear iteration loop
  virtual void NonlinearSolve();

  //! check convergence (or divergence) of nonlinear iteration
  bool AbortNonlinIter(
    const int itnum,        //!< current value of iteration step counter
    const int itmax,        //!< maximum number of iteration steps
    const double ittol,     //!< relative tolerance for increments
    const double abstolres, //!< absolute tolerance for the residual norm
    double& actresidual     //!< return value of the current residual
  );

  //! Calculate problem specific norm
  virtual void CalcProblemSpecificNorm(
      double& conresnorm,
      double& incconnorm_L2,
      double& connorm_L2,
      double& incpotnorm_L2,
      double& potnorm_L2,
      double& potresnorm,
      double& conresnorminf);

  //! Calculate the reconstructed nodal gradient of phi
  void ReconstructGradientAtNodes(
      Teuchos::RCP<Epetra_MultiVector> gradPhi);

  /// return phi at evaluation time
  virtual Teuchos::RCP<Epetra_Vector> EvaluationPhi() = 0;


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
  virtual void OutputState();

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


  /*--- query and output ---------------------------------------------------*/

  /*========================================================================*/
  //! @name scalar degrees of freedom and related
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  /*--- calculate and update -----------------------------------------------*/

  /*--- query and output ---------------------------------------------------*/



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
  //! @name  not classified method - to be kept clean!!!
  /*========================================================================*/

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
  Teuchos::RCP<std::vector<int> > writefluxids_;

  //! flux vector fields
  Teuchos::RCP<Epetra_MultiVector> flux_;

  //! vector for statistical evaluation of normal fluxes
  Teuchos::RCP<Epetra_SerialDenseVector> sumnormfluxintegral_;

  //! the last step number when fluxes have been computed
  int lastfluxoutputstep_;

  //! flag for printing out mean values of scalars
  const bool outmean_;

  //! boolean to write Gmsh postprocessing files (input parameter)
  const bool outputgmsh_;

  //! boolean to write state vectore to matlab file (input parameter)
  const bool output_state_matlab_;

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

  //! time step size
  double dta_;

  //! time measurement element
  double dtele_;
  //! time measurement solve
  double dtsolve_;

  //! number of newton iterations in actual timestep
  int iternum_;

  //! used time integration scheme
  INPAR::SCATRA::TimeIntegrationScheme timealgo_;

  /*========================================================================*/
  //! @name scalar degrees of freedom variables
  /*========================================================================*/

  //! number of space dimensions
  int nsd_;

  //! number of transported scalars
  int numscal_;

  //! number of transported scalars
  int numdofpernode_;

  //! phi at times n
  Teuchos::RCP<Epetra_Vector> phin_;
  //! phi at times n+1
  Teuchos::RCP<Epetra_Vector> phinp_;
  //! phi at meshfree nodes (for non-interpolatory meshfree basis functions)
  Teuchos::RCP<Epetra_Vector> phiatmeshfreenodes_;

  //! time derivative of phi at times n
  Teuchos::RCP<Epetra_Vector> phidtn_;
  //! time derivative of phi at times n+1
  Teuchos::RCP<Epetra_Vector> phidtnp_;

  //! histvector --- a linear combination of phinm, phin (BDF)
  //!                or phin, phidtn (One-Step-Theta)
  Teuchos::RCP<Epetra_Vector>    hist_;

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

  //! unit outer normal vector field for flux output
  Teuchos::RCP<Epetra_MultiVector> normals_;

  //! The residual vector
  Teuchos::RCP<Epetra_Vector> residual_;

  //! true (rescaled) residual vector without zeros at Dirichlet conditions
  Teuchos::RCP<Epetra_Vector> trueresidual_;

  //! Nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector> increment_;

  //! coupling of fluid-fluid at an internal interface
  Teuchos::RCP<FLD::Meshtying> meshtying_;

  // options for meshtying
  enum INPAR::FLUID::MeshTying msht_;

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

  //! number of inflow generation time steps
  int numinflowsteps_;

  /// flag for special turbulent flow
  std::string special_flow_;

  //! the vector containing source term externally computed
  //! forcing for homogeneous isotropic turbulence
  Teuchos::RCP<Epetra_Vector> forcing_;

  //! forcing for homogeneous isotropic turbulence
  Teuchos::RCP<SCATRA::HomIsoTurbScalarForcing> homisoturb_forcing_;

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

  //TODO: SCATRA_ELE_CLEANING: BIOFILM
  //! scatra fluid displacement due to growth
  Teuchos::RCP<Epetra_MultiVector> scfldgrdisp_;

  //! scatra structure displacement due to growth
  Teuchos::RCP<Epetra_MultiVector> scstrgrdisp_;

  //! flag for printing out integral values of reaction
  const bool outintegrreac_;


}; // class ScaTraTimIntImpl

} // namespace SCATRA

#endif  // #ifndef SCATRA_TIMINT_IMPLICIT_H
