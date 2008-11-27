/*!----------------------------------------------------------------------
\file fluidimplicitintegration.H
\brief Associated with control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef FLUIDIMPLICITINTEGRATION_H
#define FLUIDIMPLICITINTEGRATION_H

#include "turbulence_statistic_manager.H"

#include "../drt_fluid/dyn_smag.H"

#include "fluidimpedancecondition.H"

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
#include "../drt_lib/linalg_mapextractor.H"
#include "../drt_lib/linalg_sparsematrix.H"
#include "../drt_io/io.H"
#include "drt_periodicbc.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_function.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>


using namespace std;
using namespace Teuchos;

/// fluid algorithms
namespace FLD
{

// forwarded declarations
class TurbulenceStatisticManager;

/*!
\brief time integration for fluid problems

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class FluidImplicitTimeInt
{
  friend class TurbulenceStatisticManager;
  friend class FluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  FluidImplicitTimeInt(RCP<DRT::Discretization>  dis,
                       LINALG::Solver&           solver,
                       ParameterList&            params,
                       IO::DiscretizationWriter& output,
                       bool                      alefluid=false);


  /*!
  \brief Destructor

  */
  virtual ~FluidImplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate();

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop();

  /*!
  \brief Solve stationary problem

  */
  void SolveStationaryProblem();

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /*!
  \brief do nonlinear iteration, e.g. full Newton, Newton like or
         Fixpoint iteration

  */
  void NonlinearSolve();

  /*!
  \brief solve linearised fluid

  */
  void LinearSolve();

  /*!
    \brief build linear system matrix and rhs

    Monolithic FSI needs to access the linear fluid problem.

    \param vel new guess at velocity and pressure
  */
  void Evaluate(Teuchos::RCP<const Epetra_Vector> vel);

  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  void TimeUpdate();

  /*!
  \brief lift'n'drag forces, statistics time sample and
         output of solution and statistics

  */
  void StatisticsAndOutput();

  /*!
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief set initial flow field for analytical test problems

  */
  void SetInitialFlowField(
    int whichinitialfield,
    int startfuncno
    );

  /*!
  \brief calculate error between a analytical solution and the
         numerical solution of a test problems

  */
  void EvaluateErrorComparedToAnalyticalSol();

  /*!
  \brief calculate stresses

  */
  Teuchos::RCP<Epetra_Vector> CalcStresses();

  /*!
  \brief read restart data

  */
  void ReadRestart(int step);

  //! @name access methods for composite algorithms
  /// monolithic FSI needs to access the linear fluid problem

  Teuchos::RCP<Epetra_Vector> InitialGuess() { return incvel_; }
  Teuchos::RCP<Epetra_Vector> Residual() { return residual_; }
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return trueresidual_; }
  Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }
  Teuchos::RCP<Epetra_Vector> Velaf() { return velaf_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }
  Teuchos::RCP<Epetra_Vector> GridVel() { return gridv_; }

  /// provide access to the Dirichlet map
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /// Expand the Dirichlet DOF set
  ///
  /// The method expands the DOF set (map) which contains the DOFs
  /// subjected to Dirichlet boundary conditions. For instance, the method is 
  /// called by the staggered FSI in which the velocities on the FSI
  /// interface are prescribed by the other fields.
  void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  /// Extract the Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> Dirichlet();

  /// Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps
  ///
  /// This method provides backward compatability only. Formerly, the Dirichlet conditions
  /// were handled with the Dirichlet toggle vector. Now, they are stored and applied
  /// with maps, ie #dbcmaps_. Eventually, this method will be removed.
  const Teuchos::RCP<const Epetra_Vector> InvDirichlet();

  Teuchos::RCP<Epetra_Vector> Dispnp() { return dispnp_; }
  Teuchos::RCP<Epetra_Vector> Dispn() { return dispn_; }
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> MeshMoveMatrix() { return meshmovematrix_; }
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }
  Teuchos::RCP<Epetra_Vector> RobinRHS() { return robinrhs_; }

  LINALG::MapExtractor& VelPresSplitter() { return velpressplitter_; }
  Teuchos::RCP<const Epetra_Map> VelocityRowMap() { return velpressplitter_.OtherMap(); }
  Teuchos::RCP<const Epetra_Map> PressureRowMap() { return velpressplitter_.CondMap(); }

  double ResidualScaling() const { return density_/dta_/theta_; }
  double Dt() const { return dta_; }

  double Time() const { return time_; }
  int Step() const { return step_; }

  int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

  /*!
  \brief return type of time integration scheme

  */
  const FLUID_TIMEINTTYPE TimIntScheme() const {return timealgo_;}

  /*!
  \brief set time-step-related fields for low-Mach-number flow

  */
  void SetTimeLomaFields(RCP<const Epetra_Vector> densnp,
                         RCP<const Epetra_Vector> densn,
                         RCP<const Epetra_Vector> densnm);

  /*!
  \brief set outer-iteration-related fields for low-Mach-number flow

  */
  void SetIterLomaFields(RCP<const Epetra_Vector> densnp,
                         RCP<const Epetra_Vector> densdtnp);

  /// introduce free surface extractor object
  /*!
    This method must (and will) be called during setup with a properly
    initialized extractor object iff we are on an ale mesh.
   */
  void SetFreeSurface(const LINALG::MapExtractor* freesurface) { freesurface_ = freesurface; }

  /// introduce (robin) fsi surface extractor object
  /*!
    Introduce the FSI surface to the fluid algorithm. This way special case
    algorithms (robin coupling) can work on the surface.
   */
  void SetFSISurface(const LINALG::MapExtractor* fsisurface) { fsisurface_ = fsisurface; }

  /// determine grid velocity
  void UpdateGridv();

  /// prepare AVM3-based scale separation
  void AVM3Preparation();

  /// AVM3-based scale separation
  void AVM3Separation();

  /// calculate lift&drag forces and angular momenta
  void LiftDrag() const;

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// switch fluid field to block matrix
  void UseBlockMatrix(Teuchos::RCP<std::set<int> > condelements,
                      const LINALG::MultiMapExtractor& domainmaps,
                      const LINALG::MultiMapExtractor& rangemaps,
                      bool splitmatrix=true);

  /// linear solve with prescribed dirichlet conditions and without history
  /*!
    This is the linear solve as needed for steepest descent FSI.
   */
  void LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  //@}

  //! @name methods for turbulence models (dynamic Smagorinsky)
  /*!

  \brief filter quantities for dynamic Smagorinsky model.
         Compute averaged values for LijMij and MijMij.

  */
  void ApplyFilterForDynamicComputationOfCs();

protected:

  // don't want = operator and cctor
  FluidImplicitTimeInt operator = (const FluidImplicitTimeInt& old);
  FluidImplicitTimeInt (const FluidImplicitTimeInt& old);

  //! @name general algorithm parameters
  //! fluid discretization
  RCP<DRT::Discretization>      discret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  IO::DiscretizationWriter&     output_;
  //! the processor ID from the communicator
  int                           myrank_;
  //! do we move the fluid mesh and calculate the fluid on this moving mesh?
  bool                          alefluid_;
  bool                          dynamic_smagorinsky_;
  //@}

  //! @name time stepping variables
  double  time_;    ///< physical time
  int     step_;    ///< timestep
  int     stepmax_; ///< maximal number of timesteps
  double  maxtime_; ///< maximal physical computation time
  //@}

  /// time algorithm flag
  FLUID_TIMEINTTYPE timealgo_;

  /// constant density extracted from element material for incompressible flow
  double density_;

  //! use (or not) linearisation of reactive terms on the element
  string newton_;

  /// whether we do some extrapolation at the begin of each time step
  bool extrapolationpredictor_;

  //! predictor used in nonlinear iteration
  string  predictor_;

  //! @name restart variables
  int uprestart_;
  int upres_;
  int writestresses_;
  //@}

  //! @name time step sizes
  double dta_;
  double dtp_;

  //! @name time-integration-scheme factors
  double theta_;
  double alphaM_;
  double alphaF_;
  double gamma_;
  //@}

  //! @name parameters for sampling/dumping period
  int samstart_;
  int samstop_;
  int dumperiod_;
  //@}

  /// flag for low-Mach-number solver (no: standard incompressible solver)
  string loma_;

  /// flag for special flow: currently channel flow or flow in a lid-driven cavity
  string special_flow_;

  /// flag for potential Neumann-type outflow stabilization
  string outflow_stab_;

  /// fine-scale subgrid-viscosity flag
  string fssgv_;

  /// (standard) system matrix
  RCP<LINALG::SparseOperator> sysmat_;

  /// linearization with respect to mesh motion
  RCP<LINALG::BlockSparseMatrixBase> meshmovematrix_;

  /// maps for extracting Dirichlet and free DOF sets
  RCP<LINALG::MapExtractor> dbcmaps_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// the vector containing potential Neumann-type outflow stabilization
  RCP<Epetra_Vector>    outflow_stabil_;

  /// (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>    residual_;

  /// rhs: right hand side vector for linearised flow solution
  RefCountPtr<Epetra_Vector>    rhs_;

  /// robinrhs_: additional rhs when using robin-BC
  RCP<Epetra_Vector>    robinrhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  RCP<Epetra_Vector>    incvel_;

  //! @name accelerations and density derivatives at time n+1, n and n+alpha_M
  RCP<Epetra_Vector>    accnp_;
  RCP<Epetra_Vector>    accn_;
  RCP<Epetra_Vector>    accam_;
  //@}

  //! @name velocities and pressures at time n+1, n, n-1 and n+alpha_F
  RCP<Epetra_Vector>    velnp_;
  RCP<Epetra_Vector>    veln_;
  RCP<Epetra_Vector>    velnm_;
  RCP<Epetra_Vector>    velaf_;
  //@}

  //! @name velocities and densities at time n+1, n, n-1, n+alpha_M and n+alpha_F
  RCP<Epetra_Vector>    vedenp_;
  RCP<Epetra_Vector>    veden_;
  RCP<Epetra_Vector>    vedenm_;
  RCP<Epetra_Vector>    vedeam_;
  RCP<Epetra_Vector>    vedeaf_;
  //@}

  //! @name displacemets at time n+1, n and n-1
  RCP<Epetra_Vector>    dispnp_;
  RCP<Epetra_Vector>    dispn_;
  RCP<Epetra_Vector>    dispnm_;
  //@}

  /// grid velocity (set from the adapter!)
  RCP<Epetra_Vector>    gridv_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  /// vm3 solver
  RCP<VM3_Solver>       vm3_solver_;

  //! @name only necessary for vm3 solver: fine-scale solution vector
  RCP<Epetra_Vector>    fsvelnp_;
  //@}


  //! manager for turbulence statistics 
  RCP<TurbulenceStatisticManager> statisticsmanager_;

  //! @name Dynamic Smagorinsky model: methods and variables 
  //        -------------------------

  //! one instance of the filter object
  RCP<DynSmagFilter> DynSmag_;

  //@}

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

  /// free surface row dof map extractor
  const LINALG::MapExtractor* freesurface_;

  /// (robin) fsi coupling row dof map extractor
  const LINALG::MapExtractor* fsisurface_;

  /// Periodic Boundary Conditions. Guess what.
  RCP<PeriodicBoundaryConditions> pbc_;

  /// connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > >     pbcmapmastertoslave_;

  /// bio related special (outflow) boudaries
  RCP<UTILS::FluidImpedanceWrapper> impedancebc_;

  /// @name special relaxation state

  bool inrelaxation_;

  Teuchos::RCP<LINALG::SparseMatrix> dirichletlines_;

  Teuchos::RCP<LINALG::SparseMatrix> meshmatrix_;

  //@}

}; // class FluidImplicitTimeInt

}


#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
