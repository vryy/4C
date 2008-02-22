/*!----------------------------------------------------------------------
\file fluidimplicitintegration.H
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

#include "turbulence_statistics_ldc.H"

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
#include "../io/io_drt.H"
#include "drt_periodicbc.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_function.H"
#include "../drt_xfem/dof_management.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>


using namespace std;
using namespace Teuchos;

/*!
\class XFluidImplicitTimeInt

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
  XFluidImplicitTimeInt(
		  RCP<DRT::Discretization> dis,
		  RCP<DRT::Discretization> cutterdis,
		  LINALG::Solver&       solver,
		  ParameterList&        params,
		  IO::DiscretizationWriter& fluidoutput,
		  IO::DiscretizationWriter& solidoutput,
		  bool alefluid);


  /*!
  \brief Destructor

  */
  ~XFluidImplicitTimeInt();

  void SetupXFluidSplit(
          RCP<XFEM::DofManager> dofman,
          LINALG::MapExtractor& extractor);

  /*!
  \brief do explicit predictor step to start nonlinear iteration from
         a better initial value

  */
  void ExplicitPredictor();

    /*!
  \brief Set the part of the righthandside belonging to the last
         timestep.
  */
  void SetOldPartOfRighthandside();

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
  RCP<Epetra_Vector> CalcStresses();

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
  Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }
  Teuchos::RCP<Epetra_Vector> GridVel() { return gridv_; }
  Teuchos::RCP<Epetra_Vector> Dirichlet() { return dirichtoggle_; }
  Teuchos::RCP<Epetra_Vector> InvDirichlet() { return invtoggle_; }
  Teuchos::RCP<Epetra_Vector> Dispnp() { return dispnp_; }
  Teuchos::RCP<Epetra_Vector> Dispn() { return dispn_; }
  Teuchos::RCP<Epetra_Operator> SysMat() { return sysmat_->EpetraOperator(); }
  Teuchos::RCP<DRT::Discretization> Discretization() { return discret_; }

  Teuchos::RCP<const Epetra_Map> VelocityRowMap() { return velpressplitter_.CondMap(); }
  Teuchos::RCP<const Epetra_Map> PressureRowMap() { return velpressplitter_.OtherMap(); }

  double ResidualScaling() const { return density_/dta_/theta_; }

  double Time() const { return time_; }
  int Step() const { return step_; }

  int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

  /*!
  \brief determine grid velocity

  */
  void UpdateGridv();

  /*!
  \brief calculate lift&drag forces and angular momenta

  */
  void LiftDrag() const;

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// linear solve with prescribed dirichlet conditions and without history
  /*!
    This is the linear solve as needed for steepest descent FSI.
   */
  void LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);

  //@}

protected:

  // don't want = operator and cctor
  XFluidImplicitTimeInt operator = (const XFluidImplicitTimeInt& old);
  XFluidImplicitTimeInt (const XFluidImplicitTimeInt& old);

  RCP<DRT::Discretization> discret_;
  RCP<DRT::Discretization> cutterdiscret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  IO::DiscretizationWriter&         output_;
  IO::DiscretizationWriter&         solidoutput_;
  int                           myrank_;
  bool                          alefluid_;

  double time_;
  int step_;
  int stepmax_;
  double maxtime_;

  FLUID_TIMEINTTYPE timealgo_;

  /// density extracted from element material
  double density_;

  //! use (or not) linearisation of reactive terms on the element
  bool newton_;

  int restartstep_;
  int uprestart_;
  int writestep_;
  int upres_;
  int writestresses_;

  double dta_;
  double dtp_;
  double theta_;

  /// parameters for sampling/dumping period
  int samstart_;
  int samstop_;
  int dumperiod_;

  /// fine-scale subgrid-viscosity flag and Smagorinsky model parameter
  int fssgv_;
  double Cs_fs_;

  /// (standard) system matrix
  RCP<LINALG::SparseMatrix> sysmat_;

  /// toggle vector indicating which dofs have Dirichlet BCs
  RCP<Epetra_Vector>    dirichtoggle_;

  RCP<Epetra_Vector>    invtoggle_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RCP<Epetra_Vector>    neumann_loads_;

  /// (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>    residual_;

  /// rhs: right hand side vector for linearised flow solution
  RefCountPtr<Epetra_Vector>    rhs_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  RCP<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  RCP<Epetra_Vector>    incvel_;

  /// accelerations at time n and n-1
  RCP<Epetra_Vector>    accn_;
  RCP<Epetra_Vector>    accnm_;

  /// velocities and pressures at time n+1, n and n-1
  RCP<Epetra_Vector>    velnp_;
  RCP<Epetra_Vector>    veln_;
  RCP<Epetra_Vector>    velnm_;

  /// displacemets at time n+1, n and n-1
  RCP<Epetra_Vector>    dispnp_;
  RCP<Epetra_Vector>    dispn_;
  RCP<Epetra_Vector>    dispnm_;

  RCP<Epetra_Vector>    gridv_;

  /// histvector --- a linear combination of velnm, veln (BDF)
  ///                or veln, accn (One-Step-Theta)
  RCP<Epetra_Vector>    hist_;

  // solid displacement (might be moved somewhere else later)
  RCP<Epetra_Vector>    soliddispnp_;

  /// vm3 solver
  RCP<VM3_Solver>       vm3_solver_;

  /// only necessary for vm3 solver: fine-scale solution vector
  RCP<Epetra_Vector>    fsvelnp_;

  /// turbulence statistics
  RCP<TurbulenceStatisticsLdc> turbulencestatistics_ldc_;


  // Extractor used for convergence check
  // ------------------------------------

  LINALG::MapExtractor      velpressplitter_;


  // time measurement
  // ----------------
  RCP<Teuchos::Time> timetotal_   ;
  RCP<Teuchos::Time> timeinit_    ;
  RCP<Teuchos::Time> timetimeloop_;
  RCP<Teuchos::Time> timenlnitlin_;
  RCP<Teuchos::Time> timeelement_ ;
  RCP<Teuchos::Time> timeavm3_    ;
  RCP<Teuchos::Time> timeapplydbc_;
  RCP<Teuchos::Time> timesolver_  ;
  RCP<Teuchos::Time> timeout_     ;

  RCP<TimeMonitor> tm0_ref_ ;
  RCP<TimeMonitor> tm1_ref_ ;
  RCP<TimeMonitor> tm2_ref_ ;
  RCP<TimeMonitor> tm3_ref_ ;
  RCP<TimeMonitor> tm4_ref_ ;
  RCP<TimeMonitor> tm5_ref_ ;
  RCP<TimeMonitor> tm6_ref_ ;
  RCP<TimeMonitor> tm7_ref_ ;
  RCP<TimeMonitor> tm8_ref_ ;

  double dtele_;
  double dtsolve_;

}; // class FluidImplicitTimeInt




#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
