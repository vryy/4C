/*--------------------------------------------------------------------------*/
/*!
\file lubrication_timint_implicit.H

\brief Associated with control routine for Lubrication solvers,

     including stationary solver.

\maintainer Mostafa Faraji

\level 3

*/
/*--------------------------------------------------------------------------*/
#ifndef SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_IMPLICIT_H_
#define SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_IMPLICIT_H_

#include <Epetra_MultiVector.h>
#include <Epetra_SerialDenseVector.h>
#include <Epetra_Map.h>

#include <Teuchos_ParameterList.hpp>

#include "../drt_inpar/inpar_lubrication.H"

#include "../drt_lib/drt_discret.H"

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
  class DofSet;
  class Condition;
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

namespace FLD
{
  class DynSmagFilter;
  class Vreman;
}  // namespace FLD


/*==========================================================================*/
//! Algorithm and time integration classes for lubrication problems
/*==========================================================================*/
namespace LUBRICATION
{
  /*!
   * \brief implicit time integration for lubrication problems
   */

  class TimIntImpl
  {
   public:
    virtual Teuchos::RCP<IO::DiscretizationWriter> DiscWriter() { return output_; }

    Teuchos::RCP<Epetra_Vector>& InfGapToggle() { return inf_gap_toggle_lub_; }

    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! Standard Constructor
    TimIntImpl(Teuchos::RCP<DRT::Discretization> dis, Teuchos::RCP<LINALG::Solver> solver,
        Teuchos::RCP<Teuchos::ParameterList> params,
        Teuchos::RCP<Teuchos::ParameterList> extraparams,
        Teuchos::RCP<IO::DiscretizationWriter> output);

    //! Destructor
    virtual ~TimIntImpl();

    //! initialize time integration
    virtual void Init();

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! set the nodal film height
    void SetHeightField(const int nds, Teuchos::RCP<const Epetra_Vector> height);

    //! set the time derivative of the height (film thickness) by OST
    void SetHeightDotField(const int nds, Teuchos::RCP<const Epetra_Vector> heightdot);

    //! set average tangential interface velocity for Reynolds equation
    void SetAverageVelocityField(const int nds, Teuchos::RCP<const Epetra_Vector> av_vel);

    //! add global state vectors specific for time-integration scheme
    virtual void AddTimeIntegrationSpecificVectors(bool forcedincrementalsolver = false) = 0;

    //! prepare time loop
    virtual void PrepareTimeLoop();

    //! setup the variables to do a new time step
    virtual void PrepareTimeStep();

    //! initialization procedure prior to evaluation of first time step
    virtual void PrepareFirstTimeStep();

    //! read restart data
    virtual void ReadRestart(int step) = 0;

    /*--- calculate and update -----------------------------------------------*/

    //! do time integration (time loop)
    virtual void TimeLoop();

    //! general solver call for coupled algorithms (decides if linear/nonlinear internally)
    virtual void Solve();

    //! update the solution after convergence of the nonlinear iteration.
    virtual void Update(const int num = 0  //!< field number
        ) = 0;

    //! apply moving mesh data
    void ApplyMeshMovement(Teuchos::RCP<const Epetra_Vector> dispnp,  //!< displacement vector
        int nds  //!< number of the dofset the displacement state belongs to
    );

    //! calculate error compared to analytical solution
    virtual void EvaluateErrorComparedToAnalyticalSol();

    /*--- query and output ---------------------------------------------------*/

    //! print information about current time step to screen
    virtual void PrintTimeStepInfo();

    //! return system matrix downcasted as sparse matrix
    Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix();

    //! update Newton step
    virtual void UpdateNewton(Teuchos::RCP<const Epetra_Vector> prei);

    //! Update iteration incrementally
    //!
    //! This update is carried out by computing the new #raten_
    //! from scratch by using the newly updated #prenp_. The method
    //! respects the Dirichlet DOFs which are not touched.
    //! This method is necessary for certain predictors
    //! (like #PredictConstTempConsistRate)
    virtual void UpdateIterIncrementally() = 0;

    //! Update iteration incrementally with prescribed residual
    //! pressures
    void UpdateIterIncrementally(
        const Teuchos::RCP<const Epetra_Vector> prei  //!< input residual pressures
    );

    //! build linear system tangent matrix, rhs/force residual
    //! Monolithic EHL accesses the linearised lubrication problem
    void Evaluate();

    //! non-overlapping DOF map for multiple dofsets
    Teuchos::RCP<const Epetra_Map> DofRowMap(unsigned nds = 0)
    {
      const Epetra_Map* dofrowmap = discret_->DofRowMap(nds);
      return Teuchos::rcp(new Epetra_Map(*dofrowmap));
    }

    //! Return MapExtractor for Dirichlet boundary conditions
    const Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor() const { return dbcmaps_; }

    //! right-hand side alias the dynamic force residual
    Teuchos::RCP<const Epetra_Vector> RHS() { return residual_; }

    //! return flag indicating if an incremental solution approach is used
    bool IsIncremental() { return incremental_; }

    //! return discretization
    Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

    //! output solution and restart data to file
    virtual void Output(const int num = 0);

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
    double Dt() const { return dta_; }

    /*========================================================================*/
    //! @name pressure degrees of freedom and related
    /*========================================================================*/

    /*--- query and output ---------------------------------------------------*/

    //! return pressure field pre at time n+1
    Teuchos::RCP<Epetra_Vector> Prenp() { return prenp_; }

    //! output mean values of pressure(s)
    virtual void OutputMeanPressures(const int num = 0);

    //! output domain or boundary integrals, i.e., surface areas or volumes of specified nodesets
    void OutputDomainOrBoundaryIntegrals(const std::string condstring);

   protected:
    /*========================================================================*/
    //! @name Constructors and destructors and related methods
    /*========================================================================*/

    //! don't want = operator
    // TimIntImpl operator = (const TimIntImpl& old);

    //! don't want copy constructor
    TimIntImpl(const TimIntImpl& old);

    /*========================================================================*/
    //! @name set element parameters
    /*========================================================================*/

    virtual void SetElementTimeParameter() const = 0;

    //! set time for evaluation of Neumann boundary conditions
    virtual void SetTimeForNeumannEvaluation(Teuchos::ParameterList& params) = 0;

    //! Set general element parameters
    void SetElementGeneralParameters() const;

    /*========================================================================*/
    //! @name general framework
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/


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

    //! return the right time-scaling-factor for the true residual
    virtual double ResidualScaling() const = 0;

    //! penalty term to ensure positive pressures (cavitation)
    virtual void AddCavitationPenalty();

    //! contains the nonlinear iteration loop
    virtual void NonlinearSolve();

    //! check convergence (or divergence) of nonlinear iteration
    bool AbortNonlinIter(const int itnum,  //!< current value of iteration step counter
        const int itmax,                   //!< maximum number of iteration steps
        const double ittol,                //!< relative tolerance for increments
        const double abstolres,            //!< absolute tolerance for the residual norm
        double& actresidual                //!< return value of the current residual
    );

    //! Calculate problem specific norm
    virtual void CalcProblemSpecificNorm(
        double& preresnorm, double& incprenorm_L2, double& prenorm_L2, double& preresnorminf);

    /*--- query and output ---------------------------------------------------*/

    //! is output needed for the current time step?
    bool DoOutput() { return ((step_ % upres_ == 0) or (step_ % uprestart_ == 0)); };

    //! write state vectors prenp to BINIO
    virtual void OutputState();

    //! write state vectors prenp to Gmsh postprocessing files
    void OutputToGmsh(const int step, const double time) const;

    //! print header of convergence table to screen
    virtual void PrintConvergenceHeader();

    //! print first line of convergence table to screen
    virtual void PrintConvergenceValuesFirstIter(
        const int& itnum,            //!< current Newton-Raphson iteration step
        const int& itemax,           //!< maximum number of Newton-Raphson iteration steps
        const double& ittol,         //!< relative tolerance for Newton-Raphson scheme
        const double& preresnorm,    //!< L2 norm of pressure residual
        const double& preresnorminf  //!< infinity norm of pressure residual
    );

    //! print current line of convergence table to screen
    virtual void PrintConvergenceValues(
        const int& itnum,             //!< current Newton-Raphson iteration step
        const int& itemax,            //!< maximum number of Newton-Raphson iteration steps
        const double& ittol,          //!< relative tolerance for Newton-Raphson scheme
        const double& preresnorm,     //!< L2 norm of pressure residual
        const double& incprenorm_L2,  //!< L2 norm of pressure increment
        const double& prenorm_L2,     //!< L2 norm of pressure state vector
        const double& preresnorminf   //!< infinity norm of pressure residual
    );

    //! print finish line of convergence table to screen
    virtual void PrintConvergenceFinishLine();

    /*========================================================================*/
    //! @name Time, time-step and related methods
    /*========================================================================*/

    /*--- set, prepare, and predict ------------------------------------------*/

    //! increment time and step value
    void IncrementTimeAndStep();

    /*========================================================================*/
    //! @name general framework variables
    /*========================================================================*/

    //! solver
    Teuchos::RCP<LINALG::Solver> solver_;

    //! parameter list
    const Teuchos::RCP<Teuchos::ParameterList> params_;

    //! processor id
    int myrank_;

    //! error file handle
    FILE* errfile_;

    /*========================================================================*/
    //! @name flags and enums
    /*========================================================================*/

    //! flag for Eulerian or ALE formulation of equation(s)
    bool isale_;

    //! incremental or linear full solving? rename -> is_incremental_
    bool incremental_;

    /*--- query and output ---------------------------------------------------*/

    //! flag for printing out mean values of pressures
    const bool outmean_;

    //! boolean to write Gmsh postprocessing files (input parameter)
    const bool outputgmsh_;

    //! boolean to write state vector to matlab file (input parameter)
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

    /*========================================================================*/
    //! @name pressure degrees of freedom variables
    /*========================================================================*/

    //! number of space dimensions
    int nsd_;

    //! pressure at time n+1
    Teuchos::RCP<Epetra_Vector> prenp_;

    /*========================================================================*/
    //! @name velocity, pressure, and related
    /*========================================================================*/

    //! number of dofset associated with displacement dofs
    int nds_disp_;

    /*========================================================================*/
    //! @name Galerkin discretization, boundary conditions, and related
    /*========================================================================*/

    //! the lubrication discretization
    Teuchos::RCP<DRT::Discretization> discret_;

    //! the discretization writer
    Teuchos::RCP<IO::DiscretizationWriter> output_;

    //! system matrix (either sparse matrix or block sparse matrix)
    Teuchos::RCP<LINALG::SparseOperator> sysmat_;

    //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
    Teuchos::RCP<Epetra_Vector> zeros_;

    //! maps for extracting Dirichlet and free DOF sets
    Teuchos::RCP<LINALG::MapExtractor> dbcmaps_;

    //! the vector containing body and surface forces
    Teuchos::RCP<Epetra_Vector> neumann_loads_;

    //! residual vector
    Teuchos::RCP<Epetra_Vector> residual_;

    //! true (rescaled) residual vector without zeros at Dirichlet conditions
    Teuchos::RCP<Epetra_Vector> trueresidual_;

    //! nonlinear iteration increment vector
    Teuchos::RCP<Epetra_Vector> increment_;

    Teuchos::RCP<Epetra_Vector> prei_;  //!< residual pressures
                                        //!< \f$\Delta{p}^{<k>}_{n+1}\f$

    //! Dirchlet toggle vector for unprojectable nodes (i.e. infinite gap)
    Teuchos::RCP<Epetra_Vector> inf_gap_toggle_lub_;

    /*========================================================================*/
    //! @name not classified variables - to be kept clean!!!
    /*========================================================================*/

    //! write results every upres_ steps ? writesolutionevery_
    int upres_;

    //! write restart data every uprestart_ steps ? writesolutioneveryrestart_
    int uprestart_;

    /*========================================================================*/

  };  // class TimIntImpl
}  // namespace LUBRICATION


#endif /* SRC_DRT_LUBRICATION_LUBRICATION_TIMINT_IMPLICIT_H_ */
