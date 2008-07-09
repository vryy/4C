/*!----------------------------------------------------------------------
\file xfluidimplicitintegration.H
\brief Associated with control routine for fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Axel Gerstenberger
            gerstenberger@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef XFLUIDIMPLICITINTEGRATION_H
#define XFLUIDIMPLICITINTEGRATION_H

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
#include "../drt_lib/linalg_systemmatrix.H"
#include "../drt_io/io.H"
#include "../drt_lib/drt_function.H"
#include "../drt_xfem/dof_management.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>


using namespace Teuchos;

/*!
\brief time integration for fluid problems

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class XFluidImplicitTimeInt
{
  friend class XFluidResultTest;

public:


  /*!
  \brief Standard Constructor

  */
  XFluidImplicitTimeInt(RCP<DRT::Discretization>  dis,
                       LINALG::Solver&           solver,
                       ParameterList&            params,
                       IO::DiscretizationWriter& output,
                       bool                      alefluid=false);


  /*!
  \brief Destructor

  */
  virtual ~XFluidImplicitTimeInt();

  /*!
  \brief start time loop for startingalgo, normal problems and restarts

  */
  void Integrate(
          RCP<DRT::Discretization> cutterdiscret
          );

  /*!
  \brief Do time integration (time loop)

  */
  void TimeLoop(
          RCP<DRT::Discretization> cutterdiscret     ///< surface discretization used to intersect fluid elements
          );

  /*!
  \brief Solve stationary problem

  */
  void SolveStationaryProblem(
          RCP<DRT::Discretization> cutterdiscret     ///< surface discretization used to intersect fluid elements
          );

  /// setup the variables to do a new time step
  void PrepareTimeStep();

  /// setup the variables to do a new nonlinear iteration
  void PrepareNonlinearSolve();

  /*!
  \brief do nonlinear iteration, e.g. full Newton, Newton like or
         Fixpoint iteration

  */
  void NonlinearSolve(
          RCP<DRT::Discretization> cutterdiscret,     ///< surface discretization used to intersect fluid elements
          RCP<Epetra_Vector>       idispcol,          ///< interface displacement with coloumn distribution
          RCP<Epetra_Vector>       ivelcol,           ///< interface velocity with coloumn distribution
          RCP<Epetra_Vector>       itruerescol        ///< interface residual with coloumn distribution
          );

  /*!
  \brief compute Interface
  */
  void ComputeInterfaceAndSetDOFs(
          RCP<DRT::Discretization>  cutterdiscret,
          const Epetra_Vector&      idispcol
          );

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
  \brief update configuration and output to file/screen

  */
  void Output();

  /*!
  \brief output solution fields to Gmsh for debugging

  */
  void OutputToGmsh();

  /*!
  \brief set initial flow field for analytical test problems

  */
  void SetInitialFlowField(
    RCP<DRT::Discretization> cutterdiscret,     ///< surface discretization used to intersect fluid elements
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
  Teuchos::RCP<Epetra_Vector> Velnp() { return state_.velnp_; }
  Teuchos::RCP<Epetra_Vector> Veln() { return state_.veln_; }
  Teuchos::RCP<Epetra_Vector> Dirichlet() { return dirichtoggle_; }
  Teuchos::RCP<Epetra_Vector> InvDirichlet() { return invtoggle_; }
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }
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

protected:

  // don't want = operator and cctor
  XFluidImplicitTimeInt operator = (const XFluidImplicitTimeInt& old);
  XFluidImplicitTimeInt (const XFluidImplicitTimeInt& old);

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
  //@}

  //! @name time stepping variables
  double  time_;    ///< physical time
  int     step_;    ///< timestep
  int     stepmax_; ///< maximal number of timesteps
  double  maxtime_; ///< maximal physical computation time
  //@}

  /// time algorithm flag
  FLUID_TIMEINTTYPE timealgo_;

  /// density extracted from element material
  double density_;

  /// whether we do some extrapolation at the begin of each time step
  bool extrapolationpredictor_;

  //! @name restart variables
  int uprestart_;
  int upres_;
  int writestresses_;
  //@}

  //! @name time step sizes
  double dta_;
  double dtp_;
  double theta_;
  //@}

  /*! contains connection between elemental DofKey and a position in the global state arrays
   *  should always fit to the current Epetra_Maps in this class
   */
  XFEM::NodalDofPosMap       nodalDofDistributionMap_;
  /*! contains connection between nodal DofKey and a position in the global state arrays
   *  should always fit to the current Epetra_Maps in this class
   */
  XFEM::ElementalDofPosMap   elementalDofDistributionMap_;

  RCP<XFEM::InterfaceHandle> ihForOutput_;
  RCP<XFEM::DofManager>      dofmanagerForOutput_;
  DRT::DofSet                dofset_out_;

  /// (standard) system matrix
  RCP<LINALG::SparseOperator> sysmat_;

  /// toggle vector indicating which dofs have Dirichlet BCs
  RCP<Epetra_Vector>    dirichtoggle_;

  RCP<Epetra_Vector>    invtoggle_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>    residual_;

  /// robinrhs_: additional rhs when using robin-BC
  RCP<Epetra_Vector>    robinrhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  RCP<Epetra_Vector>    incvel_;

  struct State
  {
  //! @name accelerations at time n and n-1
  RCP<Epetra_Vector>    accn_;
  RCP<Epetra_Vector>    accnm_;
  //@}

  //! @name velocities and pressures at time n+1, n and n-1
  RCP<Epetra_Vector>    velnp_;
  RCP<Epetra_Vector>    veln_;
  RCP<Epetra_Vector>    velnm_;
  //@}

  } state_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  //! Extractor used for convergence check
  LINALG::MapExtractor      velpressplitter_;

  //! Extractor used for output
  LINALG::MapExtractor      velpressplitterForOutput_;

  /// free surface row dof map extractor
  const LINALG::MapExtractor* freesurface_;

  /// (robin) fsi coupling row dof map extractor
  const LINALG::MapExtractor* fsisurface_;

}; // class FluidImplicitTimeInt




#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
