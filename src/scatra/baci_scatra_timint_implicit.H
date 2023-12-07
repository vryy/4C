/*----------------------------------------------------------------------*/
/*! \file
\brief Associated with control routine for con.-diff.(in)stat. solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o implicit characteristic Galerkin (ICG) time-integration scheme (level-set transport)

     o explicit taylor galerkin (TG) time-integration schemes (level-set transport)

     and stationary solver.

\level 1



*----------------------------------------------------------------------*/

#ifndef BACI_SCATRA_TIMINT_IMPLICIT_H
#define BACI_SCATRA_TIMINT_IMPLICIT_H

#include "baci_adapter_scatra_wrapper.H"
#include "baci_inpar_fluid.H"
#include "baci_inpar_scatra.H"
#include "baci_io_runtime_csv_writer.H"
#include "baci_linalg_serialdensevector.H"

#include <Epetra_MultiVector.h>

#include <memory>
#include <optional>
#include <set>


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
  class Problem;
  class ResultTest;
}  // namespace DRT

namespace IO
{
  class DiscretizationReader;
  class DiscretizationWriter;
  class InputControl;
  class RuntimeCsvWriter;
}  // namespace IO

namespace CORE::LINALG
{
  class Solver;
  class SparseMatrix;
  class MapExtractor;
  class MultiMapExtractor;
  class BlockSparseMatrixBase;
  class SparseOperator;
  class KrylovProjector;
  enum class EquilibrationMethod;
  enum class MatrixType;
}  // namespace CORE::LINALG

namespace FLD
{
  class DynSmagFilter;
  class Vreman;
}  // namespace FLD

namespace SCATRA
{
  class HomIsoTurbScalarForcing;
  class MeshtyingStrategyBase;
  class ScalarHandler;
  class OutputScalarsStrategyBase;
  class OutputScalarsStrategyDomain;
  class OutputScalarsStrategyCondition;
  class OutputDomainIntegralStrategy;

  /*!
   * \brief implicit time integration for scalar transport problems
   */

  class ScaTraTimIntImpl : public ::ADAPTER::ScatraInterface
  {
    friend class HomIsoTurbInitialScalarField;
    friend class HomIsoTurbScalarForcing;
    friend class OutputScalarsStrategyBase;
    friend class OutputScalarsStrategyDomain;
    friend class OutputScalarsStrategyCondition;

   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    ScaTraTimIntImpl(Teuchos::RCP<DRT::Discretization> actdis,  //!< discretization
        Teuchos::RCP<CORE::LINALG::Solver> solver,              //!< linear solver
        Teuchos::RCP<Teuchos::ParameterList> params,            //!< parameter list
        Teuchos::RCP<Teuchos::ParameterList> extraparams,       //!< supplementary parameter list
        Teuchos::RCP<IO::DiscretizationWriter> output,          //!< output writer
        const int probnum = 0                                   //!< global problem number
    );

    //! don't want copy constructor
    ScaTraTimIntImpl(const ScaTraTimIntImpl& old) = delete;

    /*! \brief Initialize this object

    Hand in all objects/parameters/etc. from outside.
    Construct and manipulate internal objects.

    \note Try to only perform actions in Init(), which are still valid
          after parallel redistribution of discretizations.
          If you have to perform an action depending on the parallel
          distribution, make sure you adapt the affected objects after
          parallel redistribution.
          Example: cloning a discretization from another discretization is
          OK in Init(...). However, after redistribution of the source
          discretization do not forget to also redistribute the cloned
          discretization.
          All objects relying on the parallel distribution are supposed to
          the constructed in \ref Setup().

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Init();

    /*! \brief Setup all class internal objects and members

     Setup() is not supposed to have any input arguments !

     Must only be called after Init().

     Construct all objects depending on the parallel distribution and
     relying on valid maps like, e.g. the state vectors, system matrices, etc.

     Call all Setup() routines on previously initialized internal objects and members.

    \note Must only be called after parallel (re-)distribution of discretizations is finished !
          Otherwise, e.g. vectors may have wrong maps.

    \warning none
    \return void
    \date 08/16
    \author rauch  */
    virtual void Setup();

    //! Initialization of turbulence models
    void InitTurbulenceModel(const Epetra_Map* dofrowmap, const Epetra_Map* noderowmap);

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! add global state vectors specific for time-integration scheme
    void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) override;

    //! initialize system matrix
    Teuchos::RCP<CORE::LINALG::SparseOperator> InitSystemMatrix() const;

    //! prepare time loop
    virtual void PrepareTimeLoop();

    //! setup the variables to do a new time step
    virtual void PrepareTimeStep();

    //! initialization procedure prior to evaluation of first time step
    virtual void PrepareFirstTimeStep();

    //! preparations for solve
    virtual void PrepareLinearSolve();

    //! set time and step value
    virtual void SetTimeStep(const double newtime,  //!< new time value
        const int newstep                           //!< new step value
    )
    {
      time_ = newtime;
      step_ = newstep;
    }

    //! set time stepping information from this time integration to micro scales
    void SetTimeSteppingToMicroScale();

    //! set timestep value
    virtual void SetDt(const double newdt)
    {
      dta_ = newdt;
      // We have to set the newdt in the ele_calc's as well:
      SetElementTimeParameter();
    }

    //! do explicit predictor step to obtain better starting value for Newton-Raphson iteration
    virtual void ExplicitPredictor() const;

    //! set the velocity field (zero or field by function)
    virtual void SetVelocityField();

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    virtual void SetVelocityField(
        Teuchos::RCP<const Epetra_Vector> convvel,  //!< convective velocity/press. vector
        Teuchos::RCP<const Epetra_Vector> acc,      //!< acceleration vector
        Teuchos::RCP<const Epetra_Vector> vel,      //!< velocity vector
        Teuchos::RCP<const Epetra_Vector> fsvel,    //!< fine-scale velocity vector
        const bool setpressure =
            false  //!< flag whether the fluid pressure needs to be known for the scatra
    );

    void SetWallShearStresses(Teuchos::RCP<const Epetra_Vector> wss);

    void SetPressureField(Teuchos::RCP<const Epetra_Vector> pressure);

    void SetMembraneConcentration(Teuchos::RCP<const Epetra_Vector> MembraneConc);

    void SetMeanConcentration(Teuchos::RCP<const Epetra_Vector> MeanConc);

    void ClearExternalConcentrations()
    {
      mean_conc_ = Teuchos::null;
      membrane_conc_ = Teuchos::null;
    };

    //! read restart data
    virtual void ReadRestart(const int step, Teuchos::RCP<IO::InputControl> input = Teuchos::null);

    //! setup natural convection
    virtual void SetupNatConv();

    //! set number of dofset to write displacement values on
    void SetNumberOfDofSetDisplacement(int nds_disp)
    {
      dsassert(nds_disp_ == -1, "Don't set 'nds_disp_' twice!");
      nds_disp_ = nds_disp;
    }

    //! set number of dofset to write interface growth values on
    void SetNumberOfDofSetGrowth(int nds_growth)
    {
      dsassert(nds_growth_ == -1, "Don't set 'nds_growth_' twice!");
      nds_growth_ = nds_growth;
    }

    //! set number of dofset to write micro scale values on
    void SetNumberOfDofSetMicroScale(int nds_micro)
    {
      dsassert(nds_micro_ == -1, "Don't set 'nds_micro_' twice!");
      nds_micro_ = nds_micro;
    }

    //! set number of dofset to write pressure values on
    void SetNumberOfDofSetPressure(int nds_pressure)
    {
      dsassert(nds_pres_ == -1, "Don't set 'nds_pres_' twice!");
      nds_pres_ = nds_pressure;
    }

    //! set number of dofset to write scalar transport values on
    void SetNumberOfDofSetScaTra(int nds_scatra)
    {
      dsassert(nds_scatra_ == -1, "Don't set 'nds_scatra_' twice!");
      nds_scatra_ = nds_scatra;
    }

    //! set number of dofset to write thermo values on
    void SetNumberOfDofSetThermo(int nds_thermo)
    {
      dsassert(nds_thermo_ == -1, "Don't set 'nds_thermo_' twice!");
      nds_thermo_ = nds_thermo;
    }

    //! set number of dofset to write two-tensor quantities on, e.g. stresses, strains
    void SetNumberOfDofSetTwoTensorQuantity(int nds_two_tensor_quantitiy)
    {
      dsassert(nds_two_tensor_quantitiy_ == -1, "Don't set 'nds_two_tensor_quantitiy_' twice!");
      nds_two_tensor_quantitiy_ = nds_two_tensor_quantitiy;
    }

    //! set number of dofset to write velocity values on
    void SetNumberOfDofSetVelocity(int nds_velocity)
    {
      dsassert(nds_vel_ == -1, "Don't set 'nds_vel_' twice!");
      nds_vel_ = nds_velocity;
    }

    //! set number of dofset to write wall shear stress values on
    void SetNumberOfDofSetWallShearStress(int nds_wall_shear_stress)
    {
      dsassert(nds_wss_ == -1, "Don't set 'nds_wss_' twice!");
      nds_wss_ = nds_wall_shear_stress;
    }

    //! returns the maximum dofset number that is set
    [[nodiscard]] int GetMaxDofSetNumber() const;

    //! store reaction coefficient for macro-micro coupling with deforming macro dis
    void SetMacroMicroReaCoeff(const double macro_micro_rea_coeff)
    {
      macro_micro_rea_coeff_ = macro_micro_rea_coeff;
    }

    //! create result test for scalar transport field
    virtual Teuchos::RCP<DRT::ResultTest> CreateScaTraFieldTest();

    //! Add tests to global problem and start tests
    virtual void TestResults();

    /*--- calculate and update -----------------------------------------------*/

    //! do time integration (time loop)
    virtual void TimeLoop();

    //! operator for manipulations before call to \ref Solve() ; May be overridden by subclass.
    virtual void PreSolve(){};

    //! general solver call for coupled algorithms (decides if linear/nonlinear internally)
    virtual void Solve();

    //! operator for manipulations after call to \ref Solve() ; May be overridden by subclass.
    virtual void PostSolve(){};

    //! update solution after convergence of the nonlinear Newton-Raphson iteration
    virtual void Update();

    /*!
     * @brief apply moving mesh data
     *
     * @param[in] dispnp  displacement vector
     */
    void ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> dispnp);

    //! calculate fluxes inside domain and/or on boundary
    void CalcFlux(const bool writetofile  //!< flag for writing flux info to file
    );

    //! calculate flux vector field inside computational domain
    Teuchos::RCP<Epetra_MultiVector> CalcFluxInDomain();

    //! calculate mass/heat normal flux at specified boundaries and write result to file if @p
    //! writetofile is true
    Teuchos::RCP<Epetra_MultiVector> CalcFluxAtBoundary(const bool writetofile);

    //! calculation of relative error with reference to analytical solution
    virtual void EvaluateErrorComparedToAnalyticalSol();

    //! Calculate the reconstructed nodal gradient of phi from L2-projection
    Teuchos::RCP<Epetra_MultiVector> ReconstructGradientAtNodesL2Projection(
        const Teuchos::RCP<const Epetra_Vector> phi,
        bool scalenormal = false,  ///< Scale the smoothed normal field to 1
        bool returnnodal = false   ///< Return nodal based vector
    );

    //! Calculate the reconstructed nodal gradient of phi from super convergent patch recovery
    Teuchos::RCP<Epetra_MultiVector> ReconstructGradientAtNodesPatchRecon(
        const Teuchos::RCP<const Epetra_Vector> phi, const int dimension = 3,
        bool scalenormal = false,  ///< Scale the smoothed normal field to 1
        bool returnnodal = false   ///< Return nodal based vector
    );

    //! Calculate the reconstructed nodal gradient of phi from node mean averaging
    Teuchos::RCP<Epetra_MultiVector> ReconstructGradientAtNodesMeanAverage(
        Teuchos::RCP<const Epetra_Vector> phi,
        bool scalenormal = false,  ///< Scale the smoothed normal field to 1
        bool returnnodal = false   ///< Return nodal based vector
    );

    //! Calculate the reconstructed nodal gradient of phi
    void ScaleGradientsToOne(Teuchos::RCP<Epetra_MultiVector> state);

    //! finite difference check for scalar transport system matrix
    virtual void FDCheck();

    //! apply Neumann and Dirichlet BC to system
    void ApplyBCToSystem();

    void EvaluateInitialTimeDerivative(
        Teuchos::RCP<CORE::LINALG::SparseOperator> matrix, Teuchos::RCP<Epetra_Vector> rhs);

    //! prepare time integrator specific things before calculation of initial time derivative
    virtual void PreCalcInitialTimeDerivative(){};

    //! clean up settings from PreCalcInitialTimeDerivative() after initial time derivative is
    //! calculated
    virtual void PostCalcInitialTimeDerivative(){};

    //! calculate mean concentrations of micro discretization at nodes
    void CalcMeanMicroConcentration();

    /*--- query and output ---------------------------------------------------*/

    //! return ALE flag
    bool IsALE() const { return isale_; }

    //! return flag for macro scale in multi-scale simulations
    bool MacroScale() const { return macro_scale_; };

    //! return type of equilibration of global system of scalar transport equations
    CORE::LINALG::EquilibrationMethod EquilibrationMethod() const { return equilibrationmethod_; }

    //! return type of global system matrix in global system of equations
    CORE::LINALG::MatrixType MatrixType() const { return matrixtype_; }

    //! Provide enum of time integration scheme
    INPAR::SCATRA::TimeIntegrationScheme MethodName() const { return timealgo_; }

    //! Provide title of time integration scheme
    std::string MethodTitle() { return MapTimIntEnumToString(MethodName()); }

    //! return flag for micro scale in multi-scale simulations
    bool MicroScale() const { return micro_scale_; };

    //! return flag for electromagnetic diffusion simulations
    bool IsEMD() const { return isemd_; };

    //! print information about current time step to screen
    virtual void PrintTimeStepInfo();

    //! return system matrix as sparse operator
    Teuchos::RCP<CORE::LINALG::SparseOperator> SystemMatrixOperator() { return sysmat_; };

    //! return system matrix downcasted as sparse matrix
    Teuchos::RCP<CORE::LINALG::SparseMatrix> SystemMatrix();

    //! return system matrix downcasted as block sparse matrix
    Teuchos::RCP<CORE::LINALG::BlockSparseMatrixBase> BlockSystemMatrix();

    //! return map extractor associated with blocks of global system matrix
    Teuchos::RCP<const CORE::LINALG::MultiMapExtractor> BlockMaps() const { return blockmaps_; }

    //! return residual vector
    Teuchos::RCP<Epetra_Vector> Residual() const { return residual_; };

    //! return trueresidual vector
    Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }

    //! return increment vector
    Teuchos::RCP<Epetra_Vector> Increment() const { return increment_; };

    //! return flag indicating if an incremental solution approach is used
    bool IsIncremental() { return incremental_; }

    //! return Krylov projector
    Teuchos::RCP<CORE::LINALG::KrylovProjector> Projector() { return projector_; }

    //! return number of dofset associated with displacement dofs
    int NdsDisp() const override { return nds_disp_; }

    //! return number of dofset associated with interface growth dofs
    int NdsGrowth() const { return nds_growth_; }

    //! return number of dofset to store nodal micro quantities on macro dis
    int NdsMicro() const { return nds_micro_; }

    //! return number of dofset associated with pressure dofs
    int NdsPressure() const { return nds_pres_; }

    //! return number of dofset associated with scalar transport dofs
    int NdsScaTra() const { return nds_scatra_; }

    //! return number of dofset associated with thermo dofs
    int NdsThermo() const { return nds_thermo_; }

    //! return number of dofset associated with two-tensor quantity dofs, e.g. stresses, strains
    int NdsTwoTensorQuantity() const { return nds_two_tensor_quantitiy_; }

    //! return number of dofset associated with velocity dofs
    int NdsVel() const { return nds_vel_; }

    //! return number of dofset associated with wall shear stress dofs
    int NdsWallShearStress() const { return nds_wss_; }

    //! return domain flux vector
    Teuchos::RCP<const Epetra_MultiVector> FluxDomain() const { return flux_domain_; };

    //! return boundary flux vector
    Teuchos::RCP<const Epetra_MultiVector> FluxBoundary() const { return flux_boundary_; };

    //! return Dirichlet map
    Teuchos::RCP<const CORE::LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

    //! add dirichlet dofs to dbcmaps_
    void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

    //! remove dirichlet dofs from dbcmaps_
    void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

    //! return pointer to const dofrowmap
    Teuchos::RCP<const Epetra_Map> DofRowMap();

    //! return pointer to const dofrowmap of specified dofset
    Teuchos::RCP<const Epetra_Map> DofRowMap(int nds);

    //! return discretization
    Teuchos::RCP<DRT::Discretization> Discretization() const override { return discret_; }

    //! return the parameter lists
    Teuchos::RCP<Teuchos::ParameterList> ScatraParameterList() const { return params_; }
    Teuchos::RCP<Teuchos::ParameterList> ScatraExtraParameterList() { return extraparams_; }
    virtual Teuchos::RCP<Teuchos::ParameterList> ScatraTimeParameterList() = 0;

    //! Access output object: CD-Rom and DVD only - no BlueRay support!!! ;)
    const Teuchos::RCP<IO::DiscretizationWriter>& DiscWriter() const { return output_; }

    //! returns map extractor used for convergence check either in ELCH or LOMA case
    Teuchos::RCP<CORE::LINALG::MapExtractor> Splitter() const { return splitter_; }

    //! Checks if output of results or restart information is required and writes data to disk
    virtual void CheckAndWriteOutputAndRestart();

    //! write restart data to disk
    virtual void WriteRestart() const;

    //! write results to disk
    void WriteResult();

    //! Convergence check for two way coupled ScaTra problems.
    bool ConvergenceCheck(int itnum, int itmax, const double ittol);

    //! return solver
    const Teuchos::RCP<CORE::LINALG::Solver>& Solver() const { return solver_; }

    //! return parameters for finite difference check
    INPAR::SCATRA::FDCheck FDCheckType() const { return fdcheck_; };
    double FDCheckEps() const { return fdcheckeps_; };
    double FDCheckTol() const { return fdchecktol_; };

    //! return meshtying strategy (includes standard case without meshtying)
    const Teuchos::RCP<SCATRA::MeshtyingStrategyBase>& Strategy() const override
    {
      return strategy_;
    };

    //! return flag indicating availability of scatra-scatra interface kinetics condition(s)
    bool S2IKinetics() const { return s2ikinetics_; };

    //! return flag for scatra-scatra interface mesh tying
    bool S2IMeshtying() const { return s2imeshtying_; };

    //! return relative errors of scalar fields in L2 and H1 norms
    const Teuchos::RCP<std::vector<double>>& RelErrors() const { return relerrors_; };

    //! output performance statistics associated with linear solver into *.csv file
    static void OutputLinSolverStats(const CORE::LINALG::Solver& solver,  //!< linear solver
        const double& time,    //!< solver time maximized over all processors
        const int& step,       //!< time step
        const int& iteration,  //!< Newton-Raphson iteration number
        const int& size        //!< size of linear system
    );

    //! output performance statistics associated with nonlinear solver into *.csv file
    static void OutputNonlinSolverStats(
        const int& iterations,   //!< iteration count of nonlinear solver
        const double& time,      //!< solver time maximized over all processors
        const int& step,         //!< time step
        const Epetra_Comm& comm  //!< communicator
    );

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- calculate and update -----------------------------------------------*/

    //! determine whether there are still time steps to be evaluated
    [[nodiscard]] virtual bool NotFinished() const
    {
      return step_ < stepmax_ and time_ + 1.e-12 < maxtime_;
    }

    /*--- query and output ---------------------------------------------------*/

    //! return current time value
    double Time() const { return time_; }

    //! return current step number
    int Step() const { return step_; }

    //! total number of time steps ? rename StepMax?
    int NStep() const { return stepmax_; }

    //! return number of newton iterations in last timestep
    const int& IterNum() const { return iternum_; }

    //! return number of outer iterations in partitioned simulations
    const unsigned& IterNumOuter() const { return iternum_outer_; };

    //! return time step size
    double Dt() const { return dta_; }

    //! return if time step was changed during AdaptTimeStepSize()
    bool TimeStepAdapted() const { return timestepadapted_; };

    /*========================================================================*/
    //! @name scalar degrees of freedom and related
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! set the initial scalar field phi
    virtual void SetInitialField(const INPAR::SCATRA::InitialField init,  //!< type of initial field
        const int startfuncno  //!< number of spatial function
    );

    /*========================================================================*/
    //! @name Preconditioning
    /*========================================================================*/

    virtual void SetupSplitter(){};

    //! set up the (block) maps of the scatra system matrix and the meshtying object
    void SetupMatrixBlockMapsAndMeshtying();

    //! set up the (block) maps of the scatra system matrix
    void SetupMatrixBlockMaps();

    //! some of the set up of the (block) maps of the scatra system matrix has to be done after
    //! SetupMeshtying() has been called
    void PostSetupMatrixBlockMaps();

    /*!
     * @brief build maps associated with blocks of global system matrix
     *
     * @param[in]  partitioningconditions  domain partitioning conditions
     * @param[out] blockmaps               empty vector for maps to be built
     */
    virtual void BuildBlockMaps(
        const std::vector<Teuchos::RCP<DRT::Condition>>& partitioningconditions,
        std::vector<Teuchos::RCP<const Epetra_Map>>& blockmaps) const;

    //! build null spaces associated with blocks of global system matrix. Hand in solver to access
    //! parameter list and initial number of block (e.g. for coupled problems)
    virtual void BuildBlockNullSpaces(
        Teuchos::RCP<CORE::LINALG::Solver> solver, int init_block_number) const;

    /*--- calculate and update -----------------------------------------------*/

    //! call elements to calculate system matrix and rhs and assemble
    virtual void AssembleMatAndRHS();

    //! compute time derivatives of discrete state variables
    virtual void ComputeTimeDerivative();

    //! compute parameters of the Input voltage to use for the double layer current density
    virtual void ComputeTimeDerivPot0(const bool init) = 0;

    //! compute values at intermediate time steps (required for generalized-alpha) ? rename?
    virtual void ComputeIntermediateValues() = 0;

    //! compute values at the interior of the elements (required for hdg)
    virtual void ComputeInteriorValues() = 0;

    //! compute nodal density values from nodal concentration values (natural convection)
    void ComputeDensity();

    //! evaluate macro-micro coupling on micro scale in multi-scale scalar transport problems
    void EvaluateMacroMicroCoupling();

    //! iterative update of phinp
    void UpdateIter(const Teuchos::RCP<const Epetra_Vector> inc  //!< increment vector for phi
    );

    /*--- query and output ---------------------------------------------------*/

    //! return number of transported scalars
    int NumScal() const;

    //! return number of dofs per node
    int NumDofPerNode() const;

    //! return number of dofs per node in condition
    int NumDofPerNodeInCondition(const DRT::Condition& condition) const;

    //! return number of transported scalars per node in condition
    virtual int NumScalInCondition(const DRT::Condition& condition) const
    {
      return NumDofPerNodeInCondition(condition);
    };

    //! return relaxation parameters
    std::vector<double>& Omega() { return omega_; };

    //! return relaxation parameters
    const std::vector<double>& Omega() const { return omega_; };

    //! return scalar field phi at time n
    Teuchos::RCP<Epetra_Vector> Phin() override { return phin_; }

    //! return scalar field phi at time n+1
    Teuchos::RCP<Epetra_Vector> Phinp() const { return phinp_; }

    //! get mean concentration of micro discretization
    Teuchos::RCP<Epetra_Vector> PhinpMicro() const { return phinp_micro_; }

    //! return increment of scalar field phi at time n+1 for partitioned simulations
    Teuchos::RCP<Epetra_Vector>& PhinpInc() { return phinp_inc_; };

    //! return increment of scalar field phi at time n+1 for partitioned simulations
    const Teuchos::RCP<Epetra_Vector>& PhinpInc() const { return phinp_inc_; };

    //! return increment of scalar field phi at time n+1 from previous outer iteration step for
    //! partitioned simulations
    Teuchos::RCP<Epetra_Vector>& PhinpIncOld() { return phinp_inc_old_; };

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
    virtual void OutputTotalAndMeanScalars(const int num = 0);

    //! output domain or boundary integrals, i.e., surface areas or volumes of specified nodesets
    void OutputDomainOrBoundaryIntegrals(const std::string& condstring);

    //! output of reaction(s) integral
    void OutputIntegrReac(const int num = 0);

    //! return density field at time n+alpha_F (gen-alpha) or n+1 (otherwise) for natural convection
    Teuchos::RCP<Epetra_Vector> Densafnp() { return densafnp_; }

    //! problem-specific outputs
    virtual void OutputProblemSpecific(){};

    //! problem-specific restart
    virtual void ReadRestartProblemSpecific(const int step, IO::DiscretizationReader& reader){};

    //! return time for evaluation of elements
    const double& DtEle() const { return dtele_; };

    //! return time for solution of linear system of equations
    const double& DtSolve() const { return dtsolve_; };

    //! return total values of transported scalars
    const std::map<const int, std::vector<double>>& TotalScalars() const;

    //! return mean values of transported scalars
    const std::map<const int, std::vector<double>>& MeanScalars() const;

    //! return values of domain integrals
    const std::vector<double>& DomainIntegrals() const;

    //! return values of boundary integrals
    const std::vector<double>& BoundaryIntegrals() const;

    //! return micro-scale coupling flux for macro-micro coupling in multi-scale simulations
    const double& Q() const { return q_; };

    //! derivative of micro-scale coupling flux w.r.t. macro-scale state variable for macro-micro
    //! coupling in multi-scale simulations
    const std::vector<double>& DqDphi() const { return dq_dphi_; };

    //! return rcp ptr to neumann loads vector
    Teuchos::RCP<Epetra_Vector> GetNeumannLoadsPtr() override { return neumann_loads_; };

    //! returns if restart information is needed for the current time step
    [[nodiscard]] bool IsRestartStep() const
    {
      // write restart info if the simulation ends
      const bool is_finished = not NotFinished();

      return (step_ % uprestart_ == 0 and step_ != 0) or is_finished;
    }

    //! returns if output of results is needed for the current time step
    [[nodiscard]] bool IsResultStep() const { return ((step_ % upres_ == 0) or IsRestartStep()); }

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
    void SurfacePermeability(Teuchos::RCP<CORE::LINALG::SparseOperator> matrix,  //!< system matrix
        Teuchos::RCP<Epetra_Vector> rhs                                          //!< rhs vector
    );

    //! interface for fps3i problem
    void KedemKatchalsky(Teuchos::RCP<CORE::LINALG::SparseOperator> matrix,  //!< system matrix
        Teuchos::RCP<Epetra_Vector> rhs                                      //!< rhs vector
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

    //! set ptr to wrapper of this time integrator
    void SetModelEvaluatroPtr(ADAPTER::AdapterScatraWrapper* adapter_scatra_wrapper)
    {
      additional_model_evaluator_ = adapter_scatra_wrapper;
    };

   protected:
    //! create vectors for Krylov projection if necessary
    void PrepareKrylovProjection();

    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    virtual void SetElementTimeParameter(bool forcedincrementalsolver = false) const = 0;

    //! Set backward Euler time parameter
    virtual void SetElementTimeParameterBackwardEuler() const {};

    //! set time for evaluation of Neumann boundary conditions
    virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) = 0;

    //! Set general element parameters
    void SetElementGeneralParameters(bool calcinitialtimederivative = false) const;

    //! Set node set parameters
    void SetElementNodesetParameters() const;

    //! Set additional problem-specific parameters for non-standard scalar transport problems
    //! (electrochemistry etc.)
    virtual void SetElementSpecificScaTraParameters(Teuchos::ParameterList& eleparams) const {};

    //! Set element parameter specific for turbulence
    void SetElementTurbulenceParameters(bool calcinitialtimederivative = false) const;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! compute history vector, i.e., the history part of the right-hand side vector with all
    //! contributions from the previous time step
    virtual void SetOldPartOfRighthandside();

    //! create Krylov space projector
    void SetupKrylovSpaceProjection(DRT::Condition* kspcond);
    //! update Krylov space projector
    void UpdateKrylovSpaceProjection();

    //! compute approximation for fluxes and add it to a parameter list
    void AddFluxApproxToParameterList(Teuchos::ParameterList& p);

    //! calculate consistent initial scalar time derivatives in compliance with initial scalar field
    //! this function is never called directly, but only via overloading
    virtual void CalcInitialTimeDerivative();

    //! initialize meshtying strategy (including standard case without meshtying)
    virtual void CreateMeshtyingStrategy();

    /*--- calculate and update -----------------------------------------------*/

    //! apply Dirichlet boundary conditions to linear system of equations
    void ApplyDirichletToSystem();

    //! Apply Dirichlet boundary conditions on provided state vector
    virtual void ApplyDirichletBC(const double time,  //!< evaluation time
        Teuchos::RCP<Epetra_Vector> phinp,            //!< transported scalar(s) (may be = null)
        Teuchos::RCP<Epetra_Vector> phidt             //!< first time derivative (may be = null)
    );

    //! compute outward pointing unit normal vectors at given bc's
    Teuchos::RCP<Epetra_MultiVector> ComputeNormalVectors(
        const std::vector<std::string>& condnames  //!< ?
    );

    //! evaluate Neumann inflow boundary condition
    void ComputeNeumannInflow(Teuchos::RCP<CORE::LINALG::SparseOperator> matrix,  //!< ?
        Teuchos::RCP<Epetra_Vector> rhs                                           //!< ?
    );

    //! evaluate boundary condition due to convective heat transfer
    void EvaluateConvectiveHeatTransfer(Teuchos::RCP<CORE::LINALG::SparseOperator> matrix,  //!< ?
        Teuchos::RCP<Epetra_Vector> rhs                                                     //!< ?
    );

    //! potential residual scaling and potential addition of Neumann terms
    void ScalingAndNeumann();

    //! add actual Neumann loads multipl. with time factor to the residual
    virtual void AddNeumannToResidual() = 0;

    //! evaluate Neumann boundary conditions
    virtual void ApplyNeumannBC(const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
    );

    //! add parameters depending on the problem, i.e., loma, level-set, ...
    virtual void AddProblemSpecificParametersAndVectors(
        Teuchos::ParameterList& params  //!< parameter list
    );

    //! return the right time-scaling-factor for the true residual
    virtual double ResidualScaling() const = 0;

    //! solve linear system
    void LinearSolve();

    //! contains the nonlinear iteration loop
    virtual void NonlinearSolve();

    //! contains the nonlinear iteration loop for truly partitioned multi-scale simulations
    void NonlinearMultiScaleSolve();

    //! solve micro scale in truly partitioned multi-scale simulations
    void NonlinearMicroScaleSolve();

    //! Calculate the reconstructed nodal gradient of phi by means of L2-projection
    Teuchos::RCP<Epetra_MultiVector> ComputeNodalL2Projection(
        const Teuchos::RCP<const Epetra_Vector>& state, const std::string& statename,
        const int numvec, Teuchos::ParameterList& params, const int solvernumber);

    //! Calculate the reconstructed nodal gradient of phi by means of SPR
    Teuchos::RCP<Epetra_MultiVector> ComputeSuperconvergentPatchRecovery(
        Teuchos::RCP<const Epetra_Vector> state, const std::string& statename, const int numvec,
        Teuchos::ParameterList& params, const int dim);

    //! compute contributions of solution-depending boundary and interface conditions to global
    //! system of equations
    virtual void EvaluateSolutionDependingConditions(
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,  //!< system matrix
        Teuchos::RCP<Epetra_Vector> rhs                           //!< rhs vector
    );

    //! compute contribution of Robin boundary condition to eq. system
    void EvaluateRobinBoundaryConditions(
        Teuchos::RCP<CORE::LINALG::SparseOperator> matrix,  //!< system matrix
        Teuchos::RCP<Epetra_Vector> rhs                     //!< rhs vector
    );

    //! compute contributions of additional solution-depending models to global system of equations
    virtual void EvaluateAdditionalSolutionDependingModels(
        Teuchos::RCP<CORE::LINALG::SparseOperator> systemmatrix,  //!< system matrix
        Teuchos::RCP<Epetra_Vector> rhs                           //!< rhs vector
    );

    //! perform Aitken relaxation
    virtual void PerformAitkenRelaxation(Epetra_Vector& phinp,  //!< state vector to be relaxed
        const Epetra_Vector&
            phinp_inc_diff  //!< difference between current and previous state vector increments
    );

    /*--- query and output ---------------------------------------------------*/

    //! returns true if \ref Setup() was called and is still valid
    bool IsSetup() const { return issetup_; };

    //! returns true if \ ref Init() was called and is still valid
    bool IsInit() const { return isinit_; };

    //! check if \ref Setup() was called
    void CheckIsSetup() const;

    //! check if \ref Init() was called
    void CheckIsInit() const;

    //! helper function to get algorithm title
    std::string MapTimIntEnumToString(
        const enum INPAR::SCATRA::TimeIntegrationScheme term  //!< the enum
    );

    //! do we need a statistical sampling for boundary flux at the current time step?
    bool DoBoundaryFluxStatistics()
    {
      return ((step_ >= samstart_) and (step_ <= samstop_) and
              ((calcflux_boundary_ == INPAR::SCATRA::flux_total) or
                  (calcflux_boundary_ == INPAR::SCATRA::flux_diffusive) or
                  (calcflux_boundary_ == INPAR::SCATRA::flux_convective)));
    };

    //! write state vectors (phinp and convective velocity) to BINIO
    virtual void OutputState();

    //! write state vectors (phinp and convective velocity) to Gmsh postprocessing files
    void OutputToGmsh(const int step, const double time) const;

    //! write flux vectors to BINIO
    virtual void OutputFlux(Teuchos::RCP<Epetra_MultiVector> flux,  //!< flux vector
        const std::string& fluxtype  //!< flux type ("domain" or "boundary")
    );

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! adapt time step size if desired
    void AdaptTimeStepSize();

    //! compute time step size
    virtual void ComputeTimeStepSize(double& dt);

    //! increment time and step value
    void IncrementTimeAndStep();

    /*--- calculate and update -----------------------------------------------*/


    /*--- query and output ---------------------------------------------------*/

    /*========================================================================*/
    //! @name scalar degrees of freedom and related
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/
    //! compute null space information associated with global system matrix if applicable
    virtual void ComputeNullSpaceIfNecessary() const;

    //! create scalar handler
    virtual void CreateScalarHandler();

    /*--- calculate and update -----------------------------------------------*/

    /*--- query and output ---------------------------------------------------*/

    /*========================================================================*/
    //! @name AVM3 and related
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! prepare AVM3-based scale separation
    void AVM3Preparation();

    //! AVM3-based scale separation
    virtual void AVM3Separation() = 0;

    /*--- calculate and update -----------------------------------------------*/

    //! scaling of AVM3-based subgrid-diffusivity matrix
    void AVM3Scaling(Teuchos::ParameterList& eleparams  //!< parameter list
    );

    /*--- query and output ---------------------------------------------------*/

    /*========================================================================*/
    //! @name turbulence and related
    /*========================================================================*/

    //! dynamic Smagorinsk model
    virtual void DynamicComputationOfCs() = 0;

    //! dynamic Vreman model
    virtual void DynamicComputationOfCv() = 0;

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
    Teuchos::RCP<const Epetra_Vector> DirichletToggle();

    /*========================================================================*/

    /*========================================================================*/
    //! @name general framework variables
    /*========================================================================*/

    //! problem
    DRT::Problem* const problem_;

    //! problem number
    const int probnum_;

    //! solver
    Teuchos::RCP<CORE::LINALG::Solver> solver_;

    //! parameter list
    const Teuchos::RCP<Teuchos::ParameterList> params_;

    //! parameter list containing extra parameters (application dependent)
    const Teuchos::RCP<Teuchos::ParameterList> extraparams_;

    //! processor id
    int myrank_;

    //! Extractor used for convergence check either in ELCH or LOMA case
    Teuchos::RCP<CORE::LINALG::MapExtractor> splitter_;

    //! meshtying strategy (includes standard case without meshtying)
    Teuchos::RCP<SCATRA::MeshtyingStrategyBase> strategy_;

    //! Ptr to time integration wrapper.
    //! That wrapper holds a ptr to this time integrator in turn.
    //! This Ptr is uneqal nullptr only if a scatra adapter was constructed.
    //! Class AdapterScatraWrapper sets this pointer during construction
    //! by calling \ref SetModelEvaluatroPtr .
    ADAPTER::AdapterScatraWrapper* additional_model_evaluator_;

    /*========================================================================*/
    //! @name flags and enums
    /*========================================================================*/

    //! flag for Eulerian or ALE formulation of equation(s)
    bool isale_;

    //! solvertype and flags for nonlinear (always incremental) and (linear) incremental solver
    INPAR::SCATRA::SolverType solvtype_;

    //! type of equilibration of global system of scalar transport equations
    const CORE::LINALG::EquilibrationMethod equilibrationmethod_;

    //! type of global system matrix in global system of equations
    const CORE::LINALG::MatrixType matrixtype_;

    //! incremental or linear full solving? rename -> is_incremental_
    bool incremental_;

    //! flag for fine-scale subgrid-viscosity
    INPAR::SCATRA::FSSUGRDIFF fssgd_;

    //! LOMA-specific parameter: turbulence model
    INPAR::FLUID::TurbModelAction turbmodel_;

    //! flag indicating availability of scatra-scatra interface kinetics condition(s)
    bool s2ikinetics_;

    //! flag for scatra-scatra interface mesh tying
    bool s2imeshtying_;

    //! flag for artery-scatra interface coupling
    const bool arterycoupling_;

    //! flag for scatra-scatra heterogeneous reaction coupling
    const bool heteroreaccoupling_;

    //! flag for macro scale in multi-scale simulations
    const bool macro_scale_;

    //! flag for micro scale in multi-scale simulations
    const bool micro_scale_;

    //! flag for electromagnetic diffusion simulations
    bool isemd_;

    //! electromagnetic diffusion current source function
    int emd_source_;

    /*--- query and output ---------------------------------------------------*/

    //! flag for calculating flux vector field inside domain
    INPAR::SCATRA::FluxType calcflux_domain_;

    //! flag for approximate domain flux calculation involving matrix lumping
    const bool calcflux_domain_lumped_;

    //! flag for calculating flux vector field on boundary
    INPAR::SCATRA::FluxType calcflux_boundary_;

    //! flag for approximate boundary flux calculation involving matrix lumping
    const bool calcflux_boundary_lumped_;

    //! ids of scalars for which flux vectors are written (starting with 1)
    Teuchos::RCP<std::vector<int>> writefluxids_;

    //! flux vector field inside domain
    Teuchos::RCP<Epetra_MultiVector> flux_domain_;

    //! flux vector field on boundary
    Teuchos::RCP<Epetra_MultiVector> flux_boundary_;

    //! map extractor associated with boundary segments for flux calculation
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> flux_boundary_maps_;

    //! vector for statistical evaluation of normal fluxes
    Teuchos::RCP<CORE::LINALG::SerialDenseVector> sumnormfluxintegral_;

    //! the last step number when fluxes have been computed
    int lastfluxoutputstep_;

    //! flag for printing out total and mean values of transported scalars
    const INPAR::SCATRA::OutputScalarType outputscalars_;

    //! boolean to write Gmsh postprocessing files (input parameter)
    const bool outputgmsh_;

    //! boolean to write state vectore to matlab file (input parameter)
    const bool output_state_matlab_;

    //! flag for finite difference check
    const INPAR::SCATRA::FDCheck fdcheck_;

    //! perturbation magnitude for finite difference check
    const double fdcheckeps_;

    //! relative tolerance for finite difference check
    const double fdchecktol_;

    //! flag for computation of domain and boundary integrals, i.e., of surface areas and volumes
    //! associated with specified nodesets
    const INPAR::SCATRA::ComputeIntegrals computeintegrals_;

    //! flag for calculation of relative error with reference to analytical solution
    const INPAR::SCATRA::CalcError calcerror_;

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

    //! number of outer iterations in partitioned simulations
    unsigned iternum_outer_;

    //! used time integration scheme
    INPAR::SCATRA::TimeIntegrationScheme timealgo_;

    /*========================================================================*/
    //! @name scalar degrees of freedom variables
    /*========================================================================*/

    //! number of space dimensions
    int nsd_;

    //! scalar mangager
    Teuchos::RCP<ScalarHandler> scalarhandler_;

    //! scalar mangager
    Teuchos::RCP<OutputScalarsStrategyBase> outputscalarstrategy_;

    //! domain integral manager
    Teuchos::RCP<OutputDomainIntegralStrategy> outputdomainintegralstrategy_;

    //! phi at time n
    Teuchos::RCP<Epetra_Vector> phin_;
    //! phi at time n+1
    Teuchos::RCP<Epetra_Vector> phinp_;
    //! increment of phi at time n+1 for partitioned simulations
    Teuchos::RCP<Epetra_Vector> phinp_inc_;
    //! increment of phi at time n+1 from previous outer iteration step for partitioned simulations
    Teuchos::RCP<Epetra_Vector> phinp_inc_old_;
    //! relaxation parameters
    std::vector<double> omega_;

    //! time derivative of phi at time n
    Teuchos::RCP<Epetra_Vector> phidtn_;
    //! time derivative of phi at time n+1
    Teuchos::RCP<Epetra_Vector> phidtnp_;

    //! histvector --- a linear combination of phinm, phin (BDF)
    //!                or phin, phidtn (One-Step-Theta)
    Teuchos::RCP<Epetra_Vector> hist_;

    //! density at time n+alpha_F (gen-alpha) or n+1 (otherwise) for natural convection algorithm
    Teuchos::RCP<Epetra_Vector> densafnp_;

    //! relative errors of scalar fields in L2 and H1 norms
    Teuchos::RCP<std::vector<double>> relerrors_;

    /*========================================================================*/
    //! @name velocity, pressure, and related
    /*========================================================================*/

    //! subgrid-scale velocity required for multifractal subgrid-scale modeling
    Teuchos::RCP<Epetra_MultiVector> fsvel_;

    //! type of velocity field
    const INPAR::SCATRA::VelocityField velocity_field_type_;

    //! mean in time at the interface concentration
    Teuchos::RCP<const Epetra_Vector> mean_conc_;

    //! Membrane concentration in interface bewteen a scatracoupling (needed for instance for type
    //! fps3i)
    Teuchos::RCP<const Epetra_Vector> membrane_conc_;

    //! mean concentration of micro discretization  on macro dis
    Teuchos::RCP<Epetra_Vector> phinp_micro_;

   private:
    //! number of dofset associated with displacement dofs
    int nds_disp_;

    //! number of dofset associated with interface growth dofs
    int nds_growth_;

    //! number of dofset to write micro scale values on
    int nds_micro_;

    //! number of dofset associated with pressure dofs
    int nds_pres_;

    //! number of dofset associated with scatra dofs
    int nds_scatra_;

    //! number of dofset associated with thermo dofs
    int nds_thermo_;

    //! number of dofset associated with two-tensor quantity dofs, e.g. stresses, strains
    int nds_two_tensor_quantitiy_;

    //! number of dofset associated with velocity related dofs
    int nds_vel_;

    //! number of dofset associated with wall shear stress dofs
    int nds_wss_;

    /*========================================================================*/
    //! @name coefficients and related
    /*========================================================================*/
   protected:
    //! subgrid-diffusivity(-scaling) vector
    Teuchos::RCP<Epetra_Vector> subgrdiff_;

    //! densification coefficients for natural convection
    std::vector<double> densific_;

    //! initial concentrations for natural convection
    std::vector<double> c0_;

    //! reaction coefficient
    double macro_micro_rea_coeff_;

    /*========================================================================*/
    //! @name Galerkin discretization, boundary conditions, and related
    /*========================================================================*/

    //! the scalar transport discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! the discretization writer
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //! form of convective term
    INPAR::SCATRA::ConvForm convform_;

    //! system matrix (either sparse matrix or block sparse matrix)
    Teuchos::RCP<CORE::LINALG::SparseOperator> sysmat_;

    //! map extractor associated with blocks of global system matrix
    Teuchos::RCP<CORE::LINALG::MultiMapExtractor> blockmaps_;

    //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;

    //! maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<CORE::LINALG::MapExtractor> dbcmaps_;

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
    Teuchos::RCP<CORE::LINALG::SparseMatrix> sysmat_sd_;

    //! only necessary for AVM3: scale-separation matrix ? rename small caps
    Teuchos::RCP<CORE::LINALG::SparseMatrix> Sep_;

    //! only necessary for AVM3: normalized fine-scale subgrid-viscosity matrix ? rename small caps
    Teuchos::RCP<CORE::LINALG::SparseMatrix> Mnsv_;

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

    bool updateprojection_;  //!< bool triggering update of Krylov projection
    Teuchos::RCP<CORE::LINALG::KrylovProjector> projector_;  //!< Krylov projector himself

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

    //! macro-scale state variables for macro-micro coupling in multi-scale simulations
    std::vector<double> phinp_macro_;

    //! micro-scale coupling flux for macro-micro coupling in multi-scale simulations
    double q_;

    //! derivatives of micro-scale coupling flux w.r.t. macro-scale state variables for macro-micro
    //! coupling in multi-scale simulations
    std::vector<double> dq_dphi_;

    /*========================================================================*/

    /*========================================================================*/
    //! @name Biofilm specific stuff
    /*========================================================================*/

    // TODO: SCATRA_ELE_CLEANING: BIOFILM
    //! scatra fluid displacement due to growth
    Teuchos::RCP<Epetra_MultiVector> scfldgrdisp_;

    //! scatra structure displacement due to growth
    Teuchos::RCP<Epetra_MultiVector> scstrgrdisp_;

    //! flag for printing out integral values of reaction
    const bool outintegrreac_;

   private:
    /*========================================================================*/
    //! @name flags and enums
    /*========================================================================*/

    //! flag for potentially skipping computation of initial time derivative
    bool skipinitder_;

    //! flag indicating if time step was changed
    bool timestepadapted_;

    //! flag indicating if class is setup
    bool issetup_;

    //! flag indicating if class is initialized
    bool isinit_;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    //! set flag true after setup or false if setup became invalid
    void SetIsSetup(bool trueorfalse) { issetup_ = trueorfalse; };

    //! set flag true after init or false if init became invalid
    void SetIsInit(bool trueorfalse) { isinit_ = trueorfalse; };

  };  // class ScaTraTimIntImpl

  /*========================================================================*/
  /*========================================================================*/
  /*!
   * \brief Helper class for managing different number of degrees of freedom per node
   */
  class ScalarHandler
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    ScalarHandler() : numdofpernode_(), equalnumdof_(true), issetup_(false){};

    /**
     * Virtual destructor.
     */
    virtual ~ScalarHandler() = default;

    //! set up scalar handler
    virtual void Setup(const ScaTraTimIntImpl* const scatratimint);

    /*========================================================================*/
    //! @name Access and Query methods
    /*========================================================================*/

    //! return maximum number of dofs per node
    int NumDofPerNodeInCondition(const DRT::Condition& condition,
        const Teuchos::RCP<const DRT::Discretization>& discret) const;

    //! return maximum number of transported scalars per node
    virtual int NumScalInCondition(const DRT::Condition& condition,
        const Teuchos::RCP<const DRT::Discretization>& discret) const
    {
      return NumDofPerNodeInCondition(condition, discret);
    };

    //! return maximum number of dofs per node
    virtual int NumDofPerNode() const;

    //! return maximum number of transported scalars per node
    virtual int NumScal() const { return NumDofPerNode(); }

    //! return flag indicating equal number of DOFs per node in whole discretization
    bool EqualNumDof() { return equalnumdof_; };

   protected:
    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/
    //! check if \ref Setup() was called
    void CheckIsSetup() const;

    /*========================================================================*/
    //! @name Internal variables
    /*========================================================================*/
    //! number of transported scalars
    std::set<int> numdofpernode_;

    //! flag indicating equal number of DOFs per node in whole discretization
    bool equalnumdof_;

   private:
    /*========================================================================*/
    //! @name Internal variables
    /*========================================================================*/
    //! flag indicating \ref Setup() call
    bool issetup_;
  };

  /*========================================================================*/
  /*========================================================================*/
  /*!
   * \brief Base class for output of mean and total scalar values
   */
  class OutputScalarsStrategyBase
  {
   public:
    /**
     * Virtual destructor.
     */
    virtual ~OutputScalarsStrategyBase() = default;

    //! initialize time integration
    void Init(const ScaTraTimIntImpl* const scatratimint);

    //! do the output
    void OutputTotalAndMeanScalars(const ScaTraTimIntImpl* const scatratimint, const int num);

    /*========================================================================*/
    //! @name Access methods
    /*========================================================================*/

    //! return total values of transported scalars
    const std::map<const int, std::vector<double>>& TotalScalars() const { return totalscalars_; };

    //! return mean values of transported scalars
    const std::map<const int, std::vector<double>>& MeanScalars() const { return meanscalars_; };

   protected:
    /*========================================================================*/
    //! @name Helper methods
    /*========================================================================*/

    //! evaluate mean and total scalars and print them to file and screen
    virtual void EvaluateIntegrals(const ScaTraTimIntImpl* const scatratimint) = 0;

    //! print bar to screen as bottom of table
    void FinalizeScreenOutput();

    //! init objects that are specific for output strategy
    virtual void InitStrategySpecific(const ScaTraTimIntImpl* const scatratimint) = 0;

    //! Pass evaluated data to csv writer
    virtual void PassToCSVWriter() = 0;

    //! fill parameter list and set variables in discretization for evaluation of mean scalars
    void PrepareEvaluate(
        const ScaTraTimIntImpl* const scatratimint, Teuchos::ParameterList& eleparams);

    //! print header of table for summary of mean values to screen
    void PrintHeaderToScreen(const std::string& dis_name);

    //! Print evaluated data to screen
    virtual void PrintToScreen() = 0;

    /*========================================================================*/
    //! @name Internal variables
    /*========================================================================*/
    //! size of domain
    std::map<const int, double> domainintegral_;

    //! mean values of transported scalars
    std::map<const int, std::vector<double>> meanscalars_;

    //! mean values of gradient of transported scalars
    std::map<const int, std::vector<double>> meangradients_;

    //! total of micro scalars
    std::map<const int, std::vector<double>> micrototalscalars_;

    //! process number
    int myrank_;

    //! do output of mean of gradient
    bool output_mean_grad_;

    //! do output of micro dis
    bool output_micro_dis_;

    //! writes evaluated data to output
    std::optional<IO::RuntimeCsvWriter> runtime_csvwriter_;

    //! total values of transported scalars
    std::map<const int, std::vector<double>> totalscalars_;
  };

  /*!
   * \brief Strategy evaluating total and mean scalars on entire domain
   */
  class OutputScalarsStrategyDomain : virtual public OutputScalarsStrategyBase
  {
   public:
    OutputScalarsStrategyDomain() : dummy_domain_id_(-1), numdofpernode_(0), numscal_(0){};

   protected:
    void EvaluateIntegrals(const ScaTraTimIntImpl* const scatratimint) override;

    void InitStrategySpecific(const ScaTraTimIntImpl* const scatratimint) override;

    void PassToCSVWriter() override;

    void PrintToScreen() override;

   private:
    const int dummy_domain_id_;
    //! number of degrees of freedom per node
    int numdofpernode_;

    //! number of transported scalars
    int numscal_;
  };

  /*========================================================================*/
  /*========================================================================*/
  /*!
   * \brief Strategy evaluating total and mean scalars on given condition
   */
  class OutputScalarsStrategyCondition : virtual public OutputScalarsStrategyBase
  {
   public:
    OutputScalarsStrategyCondition() : conditions_(), numdofpernodepercondition_(){};

   protected:
    void EvaluateIntegrals(const ScaTraTimIntImpl* const scatratimint) override;

    void InitStrategySpecific(const ScaTraTimIntImpl* const scatratimint) override;

    void PassToCSVWriter() override;

    void PrintToScreen() override;

   private:
    //! vector of 'TotalAndMeanScalar'-conditions
    std::vector<DRT::Condition*> conditions_;

    //! number of degrees of freedom per node per 'TotalAndMeanScalar'-conditions
    std::map<int, int> numdofpernodepercondition_;

    //! number of scalars per 'TotalAndMeanScalar'-conditions
    std::map<int, int> numscalpercondition_;
  };

  /*========================================================================*/
  /*========================================================================*/
  /*!
   * \brief Strategy evaluating total and mean scalars on entire domain and on given condition
   */
  class OutputScalarsStrategyDomainAndCondition : public OutputScalarsStrategyDomain,
                                                  public OutputScalarsStrategyCondition
  {
   protected:
    void EvaluateIntegrals(const ScaTraTimIntImpl* const scatratimint) override;

    void InitStrategySpecific(const ScaTraTimIntImpl* const scatratimint) override;

    void PassToCSVWriter() override;

    void PrintToScreen() override;
  };

  /*========================================================================*/
  /*========================================================================*/
  /*!
   * \brief Strategy evaluating domain integrals on given condition
   */
  class OutputDomainIntegralStrategy
  {
   public:
    //! Standard Constructor
    OutputDomainIntegralStrategy()
        : conditionsdomain_(),
          conditionsboundary_(),
          domainintegralvalues_(),
          boundaryintegralvalues_(){};

    //! initialize time integration
    void Init(const ScaTraTimIntImpl* const scatratimint);

    //! evaluate domain integrals and print to screen
    void EvaluateIntegralsAndPrintResults(
        const ScaTraTimIntImpl* const scatratimint, const std::string& condstring);

    /*========================================================================*/
    //! @name Access methods
    /*========================================================================*/

    //! return values of domain integrals
    const std::vector<double>& DomainIntegrals() const { return domainintegralvalues_; };

    //! return values of boundary
    const std::vector<double>& BoundaryIntegrals() const { return boundaryintegralvalues_; };

   private:
    //! vector of 'DomainIntegral'-conditions
    std::vector<DRT::Condition*> conditionsdomain_;
    //! vector of 'BoundaryIntegral'-conditions
    std::vector<DRT::Condition*> conditionsboundary_;
    //! vector of 'DomainIntegral'-values
    std::vector<double> domainintegralvalues_;
    //! vector of 'BoundaryIntegral'-values
    std::vector<double> boundaryintegralvalues_;
  };

}  // namespace SCATRA
#endif  // BACI_SCATRA_TIMINT_IMPLICIT_H
