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
\maintainer Anh-Tu Vuong
            vuong@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089 - 289-15251
</pre>

*----------------------------------------------------------------------*/

#ifndef SCATRA_TIMINT_IMPLICIT_H
#define SCATRA_TIMINT_IMPLICIT_H

#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseVector.h>

#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_scatra.H"

/*==========================================================================*/
// Style guide                                                    nis Mar12
/*==========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  /*--- calculate and update -----------------------------------------------*/

  /*--- query and output ---------------------------------------------------*/



/*==========================================================================*/
// forward declarations
/*==========================================================================*/

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
  class DynSmagFilter;
  class Vreman;
}

/*==========================================================================*/
//! Algorithm and time integration classes for scalar transport problems
/*==========================================================================*/
namespace SCATRA
{

class HomIsoTurbScalarForcing;
class MeshtyingStrategyBase;

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

  //! Initialization of turbulence models
  void InitTurbulenceModel(
      const Epetra_Map* dofrowmap,
      const Epetra_Map* noderowmap);

  /*========================================================================*/
  //! @name general framework
  /*========================================================================*/

  /*--- set, prepare, and predict ------------------------------------------*/

  //! add global state vectors specific for time-integration scheme
  virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver=false) = 0;

  //! prepare time loop
  virtual void PrepareTimeLoop();

  //! setup the variables to do a new time step
  virtual void PrepareTimeStep();

  //! initialization procedure prior to evaluation of first time step
  virtual void PrepareFirstTimeStep();

  //! preparations for solve
  virtual void PrepareLinearSolve();

  //! set time and step value
  virtual void SetTimeStep(
      const double newtime, //!< new time value
      const int newstep     //!< new step value
      )
  {
    time_ = newtime;
    step_ = newstep;
  };

  //! set timestep value
  virtual void SetDt( const double newdt )
  {
    dta_=newdt;
    //We have to set the newdt in the ele_calc's as well:
    SetElementTimeParameter();
  };

  //! do explicit predictor step (-> better starting value for nonlinear solver)
  virtual void ExplicitPredictor() = 0;

  //! set the velocity field (zero or field by function)
  virtual void SetVelocityField(const int nds);

  //! set convective velocity field (+ pressure and acceleration field as
  //! well as fine-scale velocity field, if required)
  virtual void SetVelocityField(
    Teuchos::RCP<const Epetra_Vector>   convvel,                   //!< convective velocity/press. vector
    Teuchos::RCP<const Epetra_Vector>   acc,                       //!< acceleration vector
    Teuchos::RCP<const Epetra_Vector>   vel,                       //!< velocity vector
    Teuchos::RCP<const Epetra_Vector>   fsvel,                     //!< fine-scale velocity vector
    const int                           nds,                       //!< number of the dofset the velocity/pressure state belongs to
    const bool                          setpressure=false          //!< flag whether the fluid pressure needs to be known for the scatra
    );

  void SetWallShearStresses(Teuchos::RCP<const Epetra_Vector> wss, const int nds_wss);

  void SetPressureField(Teuchos::RCP<const Epetra_Vector> pressure, const int nds_pres);

  void SetMeanConcentration(Teuchos::RCP<const Epetra_Vector> MeanConc);

  void ClearMeanConcentration(){ meanconc_ = Teuchos::null; };

  //! read restart data
  virtual void ReadRestart(const int step) = 0;

  //! setup natural convection
  virtual void SetupNatConv();

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
    int nds                                   //!< number of the dofset the displacement state belongs to
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
      std::vector<std::string>&   condnames,     //!< ?
      const bool                  writetofile,   //!< ?
      const int                   num=0          //!< field number
      );

  //! calculate error compared to analytical solution
  virtual void EvaluateErrorComparedToAnalyticalSol();

    //! Calculate the reconstructed nodal gradient of phi
    Teuchos::RCP<Epetra_MultiVector>  ReconstructGradientAtNodes(
        const Teuchos::RCP<const Epetra_Vector> phi);

  //! finite difference check for scalar transport system matrix
  void FDCheck();

  /*--- query and output ---------------------------------------------------*/

  //! return ALE flag
  bool IsALE() const { return isale_; }

  //! Provide enum of time integration scheme
  INPAR::SCATRA::TimeIntegrationScheme MethodName() const { return timealgo_; }

  //! Provide title of time integration scheme
  std::string MethodTitle(){ return MapTimIntEnumToString(MethodName()); }

  //! print information about current time step to screen
  virtual void PrintTimeStepInfo();

  //! return system matrix as sparse operator
  Teuchos::RCP<LINALG::SparseOperator> SystemMatrixOperator() { return sysmat_; };

  //! return system matrix downcasted as sparse matrix
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

  //! return system matrix downcasted as block sparse matrix
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

  //! return residual vector
  Teuchos::RCP<Epetra_Vector> Residual() const { return residual_; };

  //! return trueresidual vector
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }

  //! return increment vector
  Teuchos::RCP<Epetra_Vector> Increment() const { return increment_; };

  //! return flag indicating if an incremental solution approach is used
  bool IsIncremental() { return incremental_; }

  //! return Krylov projector
  Teuchos::RCP<LINALG::KrylovProjector> Projector() { return projector_; }

  //! return number of dofset associated with velocity dofs
  int NdsVel() { return nds_vel_; }

  //! return number of dofset associated with displacement dofs
  int NdsDisp() { return nds_disp_; }

  //! return flux vector
  Teuchos::RCP<const Epetra_MultiVector> Flux() const { return flux_; };

  //! return Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  //! add dirichlet dofs to dbcmaps_
  void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  //! remove dirichlet dofs from dbcmaps_
  void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

  //! add vector to rhs
  void AddContributionToRHS(const Teuchos::RCP<const Epetra_Vector> contributing_vector);

  //! return pointer to const dofrowmap
  Teuchos::RCP<const Epetra_Map> DofRowMap();

  //! return pointer to const dofrowmap of specified dofset
  Teuchos::RCP<const Epetra_Map> DofRowMap(int nds);

  //! return discretization
  Teuchos::RCP<DRT::Discretization> Discretization() const { return discret_; }

  //! return the parameter lists
  Teuchos::RCP<Teuchos::ParameterList> ScatraParameterList() const { return params_;}
  Teuchos::RCP<Teuchos::ParameterList> ScatraExtraParameterList() { return extraparams_;}
  virtual Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList() = 0;

  //! Access output object: CD-Rom and DVD only - no BlueRay support!!! ;)
  IO::DiscretizationWriter& DiscWriter() { return *output_; }

  //! returns map extractor used for convergence check either in ELCH or LOMA case
  Teuchos::RCP<LINALG::MapExtractor> Splitter() const { return splitter_; }

  //! output solution and restart data to file
  virtual void Output(const int num =0);

  //! Convergence check for two way coupled ScaTra problems.
  bool ConvergenceCheck(int          itnum,
                        int          itmax,
                        const double ittol);

  //! return solver
  const Teuchos::RCP<LINALG::Solver>& Solver() const { return solver_; }

  //! return parameters for finite difference check
  INPAR::SCATRA::FDCheck FDCheckType() const {return fdcheck_;};
  double FDCheckEps() const {return fdcheckeps_;};
  double FDCheckTol() const {return fdchecktol_;};

  //! return meshtying strategy (includes standard case without meshtying)
  const Teuchos::RCP<SCATRA::MeshtyingStrategyBase>& Strategy() const {return strategy_;};

  //! return flag for scatra-scatra interface coupling
  const bool& S2ICoupling() const {return s2icoupling_;};

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
  const int& IterNum() const { return iternum_; }

  //! return time step size
  double Dt() const { return dta_; }

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

  //! compute time derivative
  virtual void ComputeTimeDerivative() = 0;

  //! compute parameters of the Input voltage to use for the double layer current density
  virtual void ComputeTimeDerivPot0(const bool init) = 0;

  //! compute values at intermediate time steps (required for generalized-alpha) ? rename?
  virtual void ComputeIntermediateValues() = 0;

  //! compute nodal density values from nodal concentration values (natural convection)
  void ComputeDensity();

  //! iterative update of phinp
  void UpdateIter(
    const Teuchos::RCP<const Epetra_Vector> inc //!< increment vector for phi
  );

  /*--- query and output ---------------------------------------------------*/

  //! return number of transported scalars
  int NumScal() const { return numscal_; }

  //! return number of dofs per node
  int NumDofPerNode() const { return numdofpernode_; }

  // return scalar field phi at time n-1
  virtual Teuchos::RCP<Epetra_Vector> Phinm() = 0;

  //! return scalar field phi at time n
  Teuchos::RCP<Epetra_Vector> Phin() { return phin_; }

  //! return scalar field phi at time n+1
  Teuchos::RCP<Epetra_Vector> Phinp() const { return phinp_; }

  //! return values of scalar field phi at nodes (meshfree)
  Teuchos::RCP<const Epetra_Vector> Phiatmeshfreenodes() const { return phiatmeshfreenodes_; }

  //! return time derivative of scalar field phi at time n
  Teuchos::RCP<Epetra_Vector> Phidtn() { return phidtn_; }

  //! return time derivative of scalar field phi at time n+1
  Teuchos::RCP<Epetra_Vector> Phidtnp() { return phidtnp_; }

  //! return scalar field history
  Teuchos::RCP<Epetra_Vector> Hist() { return hist_; }

  //! return scalar field phi at time n+alpha_F
  virtual Teuchos::RCP<Epetra_Vector> Phiaf() = 0;

  //! return scalar field phi at time n+alpha_F (gen-alpha) or n+1 (otherwise)
  virtual Teuchos::RCP<Epetra_Vector> Phiafnp() { return phinp_; }

  //! return scalar field phi at time n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> Phiam() = 0;

  //! return time derivative of scalar field phi at time n+alpha_M
  virtual Teuchos::RCP<Epetra_Vector> Phidtam() = 0;

  //! return fine-scale scalar field fsphi at time n+1 or alpha_M
  virtual Teuchos::RCP<Epetra_Vector> FsPhi() = 0;

  //! output total and mean values of transported scalars
  virtual void OutputTotalAndMeanScalars(const int num=0);

  //! output domain or boundary integrals, i.e., surface areas or volumes of specified nodesets
  void OutputDomainOrBoundaryIntegrals(const std::string condstring);

  //! output of reaction(s) integral
  void OutputIntegrReac(const int num=0);

  //! return density field at time n+alpha_F (gen-alpha) or n+1 (otherwise) for natural convection
  Teuchos::RCP<Epetra_Vector> Densafnp() { return densafnp_; }

  //! problem-specific outputs
  virtual void OutputProblemSpecific(){return;};

  virtual void RestartProblemSpecific( const int step ){return;};

  //! return time for evaluation of elements
  const double& DtEle() const {return dtele_;};

  //! return time for solution of linear system of equations
  const double& DtSolve() const {return dtsolve_;};

  //! return total values of transported scalars
  const std::map<const int,std::vector<double> >& TotalScalars() const {return totalscalars_;};

  //! return mean values of transported scalars
  const std::map<const int,std::vector<double> >& MeanScalars() const {return meanscalars_;};

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
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< system matrix
    Teuchos::RCP<Epetra_Vector>          rhs     //!< rhs vector
  );

  //! interface for fps3i problem
  void KedemKatchalsky(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< system matrix
    Teuchos::RCP<Epetra_Vector> rhs             //!< rhs vector
  );

  /*========================================================================*/
  //! @name Biofilm methods
  /*========================================================================*/

  //! return scatra structure growth vector
  Teuchos::RCP<const Epetra_MultiVector> StrGrowth() const { return scstrgrdisp_; };

  //! return scatra fluid growth vector
  Teuchos::RCP<const Epetra_MultiVector> FldGrowth() const { return scfldgrdisp_; };

  //! set scatra fluid displacement vector due to biofilm growth
  void SetScFldGrDisp(Teuchos::RCP<Epetra_MultiVector> scatra_fluid_growth_disp);

  //! set scatra structure displacement vector due to biofilm growth
  void SetScStrGrDisp(Teuchos::RCP<Epetra_MultiVector> scatra_struct_growth_disp);

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

  virtual void SetElementTimeParameter(bool forcedincrementalsolver=false)  const=0;

  //! Set backward Euler time parameter
  virtual void SetElementTimeParameterBackwardEuler() const {return;};

  //! set time for evaluation of Neumann boundary conditions
  virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) = 0;

  //! Set general element parameters
  void SetElementGeneralParameters(bool calcinitialtimederivative = false) const;

  //! Set additional problem-specific parameters for non-standard scalar transport problems (electrochemistry etc.)
  virtual void SetElementSpecificScaTraParameters(Teuchos::ParameterList& eleparams) const {return;};

  //! Set element parameter specific for turbulence
  void SetElementTurbulenceParameters(bool calcinitialtimederivative = false) const;

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

  //! calculate consistent initial scalar time derivatives in compliance with initial scalar field
  //! this function is never called directly, but only via overloading
  virtual void CalcInitialTimeDerivative();

  //! initialize meshtying strategy (including standard case without meshtying)
  virtual void CreateMeshtyingStrategy();

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
    const Teuchos::RCP<Epetra_Vector>&   neumann_loads //!< Neumann loads
  );

  //! call elements to calculate system matrix and rhs and assemble
  virtual void AssembleMatAndRHS();

  //! add parameters depending on the problem, i.e., loma, level-set, ...
  virtual void AddProblemSpecificParametersAndVectors(
    Teuchos::ParameterList& params //!< parameter list
  ){return;};

  //! return the right time-scaling-factor for the true residual
  virtual double ResidualScaling() const = 0;

  //! solve linear system
  void LinearSolve();

  //! contains the nonlinear iteration loop
  virtual void NonlinearSolve();

  //! Calculate the reconstructed nodal gradient of phi
  Teuchos::RCP<Epetra_MultiVector> ComputeNodalL2Projection(
    Teuchos::RCP<const Epetra_Vector> state,
    const std::string statename,
    const int numvec,
    Teuchos::ParameterList& params,
    const int solvernumber
    );

  //! compute contributions of solution-depending boundary and interface conditions to global system of equations
  virtual void EvaluateSolutionDependingConditions(
    Teuchos::RCP<LINALG::SparseOperator> systemmatrix,      //!< system matrix
    Teuchos::RCP<Epetra_Vector>          rhs                //!< rhs vector
  );

  //! compute contribution of Robin boundary condition to eq. system
  void EvaluateRobinBoundaryConditions(
    Teuchos::RCP<LINALG::SparseOperator> matrix, //!< system matrix
    Teuchos::RCP<Epetra_Vector>          rhs     //!< rhs vector
  );

  //! compute contribution of mechanical state to eq. system
  void EvaluateStructSolutionDependingConditions(
    Teuchos::RCP<Epetra_Vector>          rhs     //!< rhs vector
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

  //! adapt number of transported scalars if necessary
  virtual void AdaptNumScal(){return;};

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
   * This method provides backward compatibility only. Formerly, the Dirichlet
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

  //! solver
  Teuchos::RCP<LINALG::Solver>         solver_;

  //! parameter list
  const Teuchos::RCP<Teuchos::ParameterList> params_;

  //! parameter list containing extra parameters (application dependent)
  const Teuchos::RCP<Teuchos::ParameterList> extraparams_;

  //! processor id
  int myrank_;

  //! Extractor used for convergence check either in ELCH or LOMA case
  Teuchos::RCP<LINALG::MapExtractor> splitter_;

  //! error file handle
  FILE* errfile_;

  //! meshtying strategy (includes standard case without meshtying)
  Teuchos::RCP<SCATRA::MeshtyingStrategyBase> strategy_;

  /*========================================================================*/
  //! @name flags and enums
  /*========================================================================*/

  //! flag for Eulerian or ALE formulation of equation(s)
  bool isale_;

  //! solvertype and flags for nonlinear (always incremental) and (linear) incremental solver
  INPAR::SCATRA::SolverType solvtype_;

  //! incremental or linear full solving? rename -> is_incremental_
  bool incremental_;

  //! flag indicating if initial velocity field has been set
  bool initialvelset_;

  //! flag for fine-scale subgrid-viscosity
  INPAR::SCATRA::FSSUGRDIFF fssgd_;

  //! LOMA-specific parameter: turbulence model
  INPAR::FLUID::TurbModelAction turbmodel_;

  //! flag for scatra-scatra interface coupling
  bool s2icoupling_;

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

  //! flag for printing out total and mean values of transported scalars
  const bool outputscalars_;

  //! boolean to write Gmsh postprocessing files (input parameter)
  const bool outputgmsh_;

  //! boolean to write state vectore to matlab file (input parameter)
  const bool output_state_matlab_;

  //! flag for finite difference check
  INPAR::SCATRA::FDCheck fdcheck_;

  //! perturbation magnitude for finite difference check
  const double fdcheckeps_;

  //! relative tolerance for finite difference check
  const double fdchecktol_;

  //! flag for computation of domain and boundary integrals, i.e., of surface areas and volumes associated with specified nodesets
  INPAR::SCATRA::ComputeIntegrals computeintegrals_;

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

  //! phi at time n
  Teuchos::RCP<Epetra_Vector> phin_;
  //! phi at time n+1
  Teuchos::RCP<Epetra_Vector> phinp_;
  //! phi at meshfree nodes (for non-interpolatory meshfree basis functions)
  Teuchos::RCP<Epetra_Vector> phiatmeshfreenodes_;

  //! time derivative of phi at time n
  Teuchos::RCP<Epetra_Vector> phidtn_;
  //! time derivative of phi at time n+1
  Teuchos::RCP<Epetra_Vector> phidtnp_;

  //! histvector --- a linear combination of phinm, phin (BDF)
  //!                or phin, phidtn (One-Step-Theta)
  Teuchos::RCP<Epetra_Vector> hist_;

  //! density at time n+alpha_F (gen-alpha) or n+1 (otherwise) for natural convection algorithm
  Teuchos::RCP<Epetra_Vector> densafnp_;

  //! total values of transported scalars
  std::map<const int,std::vector<double> > totalscalars_;

  //! mean values of transported scalars
  std::map<const int,std::vector<double> > meanscalars_;

  /*========================================================================*/
  //! @name velocity, pressure, and related
  /*========================================================================*/

  //! subgrid-scale velocity required for multifractal subgrid-scale modeling
  Teuchos::RCP<Epetra_MultiVector> fsvel_;

  //! type of velocity field ? rename!!!
  const INPAR::SCATRA::VelocityField cdvel_;

  //! mean concentration in interface bewteen a scatracoupling (needed for instance for type fps3i)
  Teuchos::RCP<const Epetra_Vector> meanconc_;

  //! number of dofset associated with velocity related dofs
  int nds_vel_;

  //! number of dofset associated with displacement dofs
  int nds_disp_;

  //! number of dofset associated with pressure dofs
  int nds_pres_;

  //! number of dofset associated with wall shear stress dofs
  int nds_wss_;

  /*========================================================================*/
  //! @name coefficients and related
  /*========================================================================*/

  //! subgrid-diffusivity(-scaling) vector
  Teuchos::RCP<Epetra_Vector> subgrdiff_;

  //! densification coefficients for natural convection
  std::vector<double> densific_;

  //! initial concentrations for natural convection
  std::vector<double> c0_;

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

  //! residual vector
  Teuchos::RCP<Epetra_Vector> residual_;

  //! true (rescaled) residual vector without zeros at Dirichlet conditions
  Teuchos::RCP<Epetra_Vector> trueresidual_;

  //! nonlinear iteration increment vector
  Teuchos::RCP<Epetra_Vector> increment_;

  //! options for meshtying
  INPAR::FLUID::MeshTying msht_;

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

  //! flag for potentially skipping computation of initial time derivative
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
