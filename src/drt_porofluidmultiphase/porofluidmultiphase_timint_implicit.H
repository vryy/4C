/*----------------------------------------------------------------------*/
/*!
 \brief base class for implicit time integration schemes for
        multiphas porous fluid problems

   \level 3

   \maintainer  Johannes Kremheller
 *----------------------------------------------------------------------*/

#ifndef SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_TIMINT_IMPLICIT_H_
#define SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_TIMINT_IMPLICIT_H_



#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_Map.h>

#include <Teuchos_ParameterList.hpp>

#include "../drt_lib/drt_discret.H"

#include "../drt_adapter/ad_porofluidmultiphase.H"
#include "../drt_adapter/ad_art_net.H"

#include "../drt_inpar/inpar_porofluidmultiphase.H"


/*==========================================================================*/
// forward declarations
/*==========================================================================*/

namespace DRT
{
  class DofSet;
  class Condition;
  class ResultTest;
}  // namespace DRT

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
}  // namespace LINALG

namespace ADAPTER
{
  class ArtNet;
}


namespace POROFLUIDMULTIPHASE
{
  // forward declaration
  class MeshtyingStrategyBase;

  /*!
   * \brief implicit time integration for porous multiphase flow problems
   */

  class TimIntImpl : public ADAPTER::PoroFluidMultiphase
  {
   public:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimIntImpl(Teuchos::RCP<DRT::Discretization> dis, const int linsolvernumber,
        const Teuchos::ParameterList& probparams, const Teuchos::ParameterList& poroparams,
        FILE* errfile, Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~TimIntImpl();

    //! initialize time integration
    virtual void Init(bool isale, int nds_disp, int nds_vel, int nds_solidpressure, int nds_scalar,
        const std::map<int, std::set<int>>* nearbyelepairs);

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! add global state vectors specific for time-integration scheme
    virtual void AddTimeIntegrationSpecificVectors() = 0;

    //! prepare time loop
    virtual void PrepareTimeLoop();

    //! setup the variables to do a new time step
    virtual void PrepareTimeStep();

    //! initialization procedure prior to evaluation of first time step
    virtual void PrepareFirstTimeStep();

    //! initialization procedure prior to evaluation of first time step
    virtual void CalcInitialTimeDerivative();

    //! read restart data
    virtual void ReadRestart(int step);

    /// create result test for porous fluid field
    Teuchos::RCP<DRT::ResultTest> CreateFieldTest();

    //! finite difference check for system matrix
    void FDCheck();

    /*--- calculate and update -----------------------------------------------*/

    //! do time integration (time loop)
    virtual void TimeLoop();

    //! general solver call for coupled algorithms
    virtual void Solve();

    //! update the solution after convergence of the nonlinear iteration.
    virtual void Update();

    ///  compute time derivative
    virtual void ComputeTimeDerivative() = 0;

    ///  compute intermediate values if necessary
    virtual void ComputeIntermediateValues() = 0;

    //! apply moving mesh data
    void ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> dispnp  //!< displacement vector
    );

    //! set convective velocity field (+ pressure and acceleration field as
    //! well as fine-scale velocity field, if required)
    void SetVelocityField(Teuchos::RCP<const Epetra_Vector> vel  //!< velocity vector
    );

    //! set state on discretization
    void SetState(unsigned nds, const std::string& name, Teuchos::RCP<const Epetra_Vector> state);

    //! calculate error compared to analytical solution
    virtual void EvaluateErrorComparedToAnalyticalSol();

    /*--- query and output ---------------------------------------------------*/

    //! print information about current time step to screen
    virtual void PrintTimeStepInfo();

    //! iterative update of phinp
    virtual void UpdateIter(
        const Teuchos::RCP<const Epetra_Vector> inc  //!< increment vector for phi
    );

    //! build linear system tangent matrix, rhs/force residual
    void Evaluate();

    //! apply Dirichlet Boundary Condition
    void PrepareSystemForNewtonSolve();

    //! direct access to system matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix()
    {
      return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_);
    };

    //! Return MapExtractor for Dirichlet boundary conditions
    const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const
    {
      return dbcmaps_with_volfracpress_;
    }

    //! right-hand side alias the dynamic force residual
    Teuchos::RCP<const Epetra_Vector> RHS() const { return residual_; }

    //! right-hand side alias the dynamic force residual for coupled system
    Teuchos::RCP<const Epetra_Vector> ArteryPorofluidRHS() const;

    //! return discretization
    Teuchos::RCP<DRT::Discretization> Discretization() const { return discret_; }

    //! access dof row map
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds) const;

    //! access dof row map
    Teuchos::RCP<const Epetra_Map> ArteryDofRowMap() const;

    //! direct access to block system matrix of artery poro problem
    virtual Teuchos::RCP<LINALG::BlockSparseMatrixBase> ArteryPorofluidSysmat() const;

    //! output solution and restart data to file
    virtual void Output();

    //! output solution and restart data to file
    virtual void PrintHeader();

    //! write additional data required for restart
    virtual void OutputRestart() = 0;

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- query and output ---------------------------------------------------*/

    //! return current time value
    double Time() const { return time_; }

    //! return current step number
    int Step() const { return step_; }

    //! return number of newton iterations in last timestep
    double IterNum() const { return iternum_; }

    //! return time step size
    double Dt() const { return dt_; }

    /*========================================================================*/
    //! @name degrees of freedom and related
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! set the initial scalar field phi
    virtual void SetInitialField(
        const INPAR::POROFLUIDMULTIPHASE::InitialField init,  //!< type of initial field
        const int startfuncno                                 //!< number of spatial function
    );

    /*--- query and output ---------------------------------------------------*/

    //! return pressure field at time n+1
    Teuchos::RCP<const Epetra_Vector> Phinp() const { return phinp_; }

    //! return scalar field phi at time n
    Teuchos::RCP<const Epetra_Vector> Phin() const { return phin_; }

    //! return time derivative of scalar field phi at time n
    Teuchos::RCP<const Epetra_Vector> Phidtn() const { return phidtn_; }

    //! return time derivative of scalar field phi at time n+1
    Teuchos::RCP<const Epetra_Vector> Phidtnp() const { return phidtnp_; }

    //! return scalar field history
    Teuchos::RCP<const Epetra_Vector> Hist() const { return hist_; }

    //! return solid pressure field
    Teuchos::RCP<const Epetra_Vector> SolidPressure() const { return solidpressure_; }

    //! return pressure field
    Teuchos::RCP<const Epetra_Vector> Pressure() const { return pressure_; }

    //! return saturation field
    Teuchos::RCP<const Epetra_Vector> Saturation() const { return saturation_; }

    //! return phase flux field at time n+1
    Teuchos::RCP<const Epetra_MultiVector> Flux() const { return flux_; }

    //! return number of dof set associated with solid pressure
    int GetDofSetNumberOfSolidPressure() const { return nds_solidpressure_; };

    //! return valid volume fraction species
    Teuchos::RCP<const Epetra_Vector> ValidVolFracSpecDofs() const
    {
      return valid_volfracspec_dofs_;
    }

   protected:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! don't want copy constructor
    TimIntImpl(const TimIntImpl& old);

    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    //! Set element time step parameters (varying every time step)
    virtual void SetElementTimeStepParameter() const = 0;

    //! set time for evaluation of Neumann boundary conditions
    virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) = 0;

    //! Set general element parameters
    void SetElementGeneralParameters() const;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! Set the part of the righthandside belonging to the last timestep.
    virtual void SetOldPartOfRighthandside() = 0;

    /*--- calculate and update -----------------------------------------------*/

    //! Apply Dirichlet boundary conditions on provided state vector
    void ApplyDirichletBC(const double time,  //!< evaluation time
        Teuchos::RCP<Epetra_Vector> prenp,    //!< pressure (may be = null)
        Teuchos::RCP<Epetra_Vector> predt     //!< first time derivative (may be = null)
    );

    //! potential residual scaling and potential addition of Neumann terms
    void ScalingAndNeumann();

    //! add actual Neumann loads multipl. with time factor to the residual
    virtual void AddNeumannToResidual() = 0;

    //! Apply Neumann boundary conditions
    void ApplyNeumannBC(const Teuchos::RCP<Epetra_Vector>& neumann_loads  //!< Neumann loads
    );

    //! call elements to calculate system matrix and rhs and assemble
    virtual void AssembleMatAndRHS();

    //! call elements to find the valid volume frac pressures and species
    virtual void EvaluateValidVolumeFracPressAndSpec();

    //! apply the additional volume fraction pressures as DBC
    virtual void ApplyAdditionalDBCForVolFracPress();

    //! call elements to calculate fluid coupling matrix with structure and assemble
    virtual void AssembleFluidStructCouplingMat(Teuchos::RCP<LINALG::SparseOperator> k_fs);

    //! call elements to calculate fluid coupling matrix with scatra and assemble
    virtual void AssembleFluidScatraCouplingMat(Teuchos::RCP<LINALG::SparseOperator> k_pfs);

    //! return the right time-scaling-factor for the true residual
    virtual double ResidualScaling() const = 0;

    //! contains the nonlinear iteration loop
    virtual void NonlinearSolve();

    //! check convergence (or divergence) of nonlinear iteration
    bool AbortNonlinIter(const int itnum,  //!< current value of iteration step counter
        const int itmax,                   //!< maximum number of iteration steps
        const double abstolres,            //!< absolute tolerance for the residual norm
        double& actresidual                //!< return value of the current residual
    );

    //! linear solve
    virtual void LinearSolve(bool isadapttol, double actresidual, double adaptolbetter);

    //! reconstruct pressures and saturation from current solution
    void ReconstructPressuresAndSaturations();

    //! reconstruct solid pressures from current solution
    void ReconstructSolidPressures();

    //! reconstruct fluxes from current solution
    void ReconstructFlux();

    //! reconstruct porosity from current solution
    void ReconstructPorosity();

    //! evaluate domain integrals
    void EvaluateDomainIntegrals();

    /*--- query and output ---------------------------------------------------*/

    //! is output needed for the current time step?
    bool DoOutput() { return ((step_ % upres_ == 0) or (step_ % uprestart_ == 0)); };

    //! write state vectors prenp to BINIO
    virtual void OutputState();

    //! print header of convergence table to screen
    virtual void PrintConvergenceHeader();

    //! print first line of convergence table to screen
    virtual void PrintConvergenceValuesFirstIter(
        const int& itnum,         //!< current Newton-Raphson iteration step
        const int& itemax,        //!< maximum number of Newton-Raphson iteration steps
        const double& ittol,      //!< relative tolerance for Newton-Raphson scheme
        const double& preresnorm  //!< norm of pressure residual
    );

    //! print current line of convergence table to screen
    virtual void PrintConvergenceValues(
        const int& itnum,          //!< current Newton-Raphson iteration step
        const int& itemax,         //!< maximum number of Newton-Raphson iteration steps
        const double& ittol,       //!< relative tolerance for Newton-Raphson scheme
        const double& preresnorm,  //!< norm of pressure residual
        const double& incprenorm,  //!< norm of pressure increment
        const double& prenorm      //!< norm of pressure state vector
    );

    //! print finish line of convergence table to screen
    virtual void PrintConvergenceFinishLine();

    // return arterial network time integrator
    virtual Teuchos::RCP<ADAPTER::ArtNet> ArtNetTimInt();

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! increment time and step value
    void IncrementTimeAndStep();

    /*========================================================================*/
    //! @name general framework variables
    /*========================================================================*/

    //! linear solver
    Teuchos::RCP<LINALG::Solver> solver_;

    //! solver number in input file
    const int linsolvernumber_;

    //! parameter list of global control problem
    const Teuchos::ParameterList& params_;

    //! parameter list of poro fluid multiphase problem
    const Teuchos::ParameterList& poroparams_;

    //! processor id
    int myrank_;

    //! error file handle
    FILE* errfile_;

    //! number of space dimensions
    int nsd_;

    /*========================================================================*/
    //! @name flags and enums
    /*========================================================================*/

    //! flag for Eulerian or ALE formulation of equation(s)
    bool isale_;

    //! flag if initial time derivative should be skipped
    bool skipinitder_;

    //! flag if porosity should be output
    bool output_porosity_;

    //! flag for biot stabilization
    bool stab_biot_;

    /*--- query and output ---------------------------------------------------*/

    //! parameters for domain integration
    std::vector<int> domainint_funct_;
    int num_domainint_funct_;

    //! flag for error calculation
    const INPAR::POROFLUIDMULTIPHASE::CalcError calcerr_;

    //! flag for flux reconstruction
    const INPAR::POROFLUIDMULTIPHASE::FluxReconstructionMethod fluxrecon_;

    //! solver number for flux reconstruction
    const int fluxreconsolvernum_;

    //! what to do when nonlinear solution fails
    enum INPAR::POROFLUIDMULTIPHASE::DivContAct divcontype_;

    //! flag for finite difference check
    const INPAR::POROFLUIDMULTIPHASE::FDCheck fdcheck_;

    //! perturbation magnitude for finite difference check
    const double fdcheckeps_;

    //! relative tolerance for finite difference check
    const double fdchecktol_;

    //! scaling factor for biot stabilization
    double stab_biot_scaling_;

    /*========================================================================*/
    //! @name Time, time-step, and iteration variables
    /*========================================================================*/

    //! actual time
    double time_;

    //! maximum simulation time
    double maxtime_;

    //! actual step number
    int step_;

    //! maximum number of steps
    const int stepmax_;

    //! time step size
    double dt_;

    //! time measurement element
    double dtele_;

    //! time measurement solve
    double dtsolve_;

    //! number of newton iterations in actual timestep
    int iternum_;

    //! maximum number of newton iterations
    const int itemax_;

    //! write results every upres_ steps ? writesolutionevery_
    const int upres_;

    //! write restart data every uprestart_ steps ? writesolutioneveryrestart_
    const int uprestart_;

    // vector norm for residuals
    enum INPAR::POROFLUIDMULTIPHASE::VectorNorm vectornormfres_;
    // vector norm for increments
    enum INPAR::POROFLUIDMULTIPHASE::VectorNorm vectornorminc_;

    //! convergence tolerance for increments
    double ittolres_;
    //! convergence tolerance for residuals
    double ittolinc_;

    /*========================================================================*/
    //! @name degrees of freedom variables
    /*========================================================================*/

    //! phi at time n
    Teuchos::RCP<Epetra_Vector> phin_;
    //! phi at time n+1
    Teuchos::RCP<Epetra_Vector> phinp_;

    //! time derivative of phi at time n
    Teuchos::RCP<Epetra_Vector> phidtn_;
    //! time derivative of phi at time n+1
    Teuchos::RCP<Epetra_Vector> phidtnp_;

    //! histvector --- a linear combination of phinm, phin (BDF)
    //!                or phin, phidtn (One-Step-Theta)
    Teuchos::RCP<Epetra_Vector> hist_;

    /*========================================================================*/
    //! @name degrees of freedom and related
    /*========================================================================*/

    //! pressure at time n+1
    Teuchos::RCP<Epetra_Vector> pressure_;

    //! saturation at time n+1
    Teuchos::RCP<Epetra_Vector> saturation_;

    //! solid pressure at time n+1
    Teuchos::RCP<Epetra_Vector> solidpressure_;

    //! porosity at time n+1
    Teuchos::RCP<Epetra_Vector> porosity_;

    //! vector with valid volume fraction pressure dofs, this vector identifies volume fraction
    //! pressure DOFs,
    //  which actually have to be evaluated with a double >= 1.0, see also
    //  EvaluatorValidVolFracPressures: if at least one nodal volume fraction value of an element is
    //  bigger than a threshold (min volfrac), the volume fraction pressure is a valid (physically
    //  meaningful) quantity in this element and the respective Darcy equation has to be solved
    //  for volume fraction species we only evaluate if all nodal volume fraction values of the
    //  element are bigger than the threshold (min volfrac), this turned out to be the most stable
    //  approach
    Teuchos::RCP<Epetra_Vector> valid_volfracpress_dofs_;
    Teuchos::RCP<Epetra_Vector> valid_volfracspec_dofs_;

    //! velocity of each phase at time n+1 (post-processed from pressure solution)
    Teuchos::RCP<Epetra_MultiVector> flux_;

    //! number of dofset associated with displacement dofs
    int nds_disp_;

    //! number of dofset associated with velocity related dofs
    int nds_vel_;

    //! number of dofset associated with solid pressure dofs
    int nds_solidpressure_;

    //! number of dofset associated with scatra dofs
    int nds_scatra_;

    /*========================================================================*/
    //! @name Galerkin discretization, boundary conditions, and related
    /*========================================================================*/

    //! the porous multiphase flow discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! the discretization writer
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //! system matrix (either sparse matrix or block sparse matrix)
    Teuchos::RCP<LINALG::SparseOperator> sysmat_;

    //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;

    //! maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

    //! maps for extracting Dirichlet and free DOF sets, here the additional dofs which have to be
    //! zeroed out for
    //  the volume fraction pressure since it is not defined if the corresponding volume fraction is
    //  equal to zero ( or smaller than minvolfrac) have been added
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_with_volfracpress_;

    //! the vector containing body and surface forces
    Teuchos::RCP<Epetra_Vector> neumann_loads_;

    //! residual vector
    Teuchos::RCP<Epetra_Vector> residual_;

    //! true (rescaled) residual vector without zeros at Dirichlet conditions
    Teuchos::RCP<Epetra_Vector> trueresidual_;

    //! nonlinear iteration increment vector
    Teuchos::RCP<Epetra_Vector> increment_;

    //! meshtying strategy (includes standard case without meshtying)
    Teuchos::RCP<POROFLUIDMULTIPHASE::MeshtyingStrategyBase> strategy_;

    /*========================================================================*/

  };  // class TimIntImpl
}  // namespace POROFLUIDMULTIPHASE


#endif /* SRC_DRT_POROFLUIDMULTIPHASE_POROFLUIDMULTIPHASE_TIMINT_IMPLICIT_H_ */
