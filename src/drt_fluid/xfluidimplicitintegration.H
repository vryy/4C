/*!----------------------------------------------------------------------
\file xfluidimplicitintegration.H
\brief Control routine for fluid time integration. Associated with
       fluidimplicitintegration.cpp


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
#include "../io/io_drt.H"
#include "../drt_xfem/dof_management.H"
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
#ifdef RESULTTEST
  friend class XFluidResultTest;
#endif

public:


  /*!
  \brief Standard Constructor

  */
  XFluidImplicitTimeInt(
		  RefCountPtr<DRT::Discretization> dis,
		  RefCountPtr<DRT::Discretization> cutterdis,
		  LINALG::Solver&       solver,
		  ParameterList&        params,
		  IO::DiscretizationWriter& fluidoutput,
		  IO::DiscretizationWriter& solidoutput,
		  bool alefluid);


  /*!
  \brief Destructor

  */
  virtual ~XFluidImplicitTimeInt();


  /*!
  \brief Create default parameter list

  */
  static void SetDefaults(ParameterList& params);

  /*!
  \brief Create single field rowmaps used for error control of iteration 

  */
  void ComputeSingleFieldRowMaps(RCP<XFEM::DofManager> dofman);
  
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
  \brief compute convergence check after a nonlinear iteration step

  */
  void NonlinearConvCheck(
    bool&  stopnonliniter,
    int    itnum         ,
    double dtele         ,
    double dtsolve
    );

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
  RefCountPtr<Epetra_Vector> CalcStresses();
  
  /*!
  \brief read restart data

  */
  void ReadRestart(int step);


protected:

  // don't want = operator and cctor
  XFluidImplicitTimeInt operator = (const XFluidImplicitTimeInt& old);
  XFluidImplicitTimeInt (const XFluidImplicitTimeInt& old);

  RefCountPtr<DRT::Discretization> discret_;
  RefCountPtr<DRT::Discretization> cutterdiscret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  IO::DiscretizationWriter&         output_;
  IO::DiscretizationWriter&         solidoutput_;
  int                           maxentriesperrow_;
  bool                          alefluid_;

  double time_;
  int step_;
  int stepmax_;
  double maxtime_;

  FLUID_TIMEINTTYPE timealgo_;

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

  /// system matrix --- stiffness and mass are assembled in one matrix!
  RefCountPtr<Epetra_CrsMatrix> sysmat_;

  /// toggle vector indicating which dofs have Dirichlet BCs
  RefCountPtr<Epetra_Vector>    dirichtoggle_;

  RefCountPtr<Epetra_Vector>    invtoggle_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RefCountPtr<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RefCountPtr<Epetra_Vector>    neumann_loads_;

  /// The residual vector --- more or less the rhs for the incremental
  /// formulation!!!
  RefCountPtr<Epetra_Vector>    residual_;

  /// true (rescaled) residual vector without zeros at dirichlet positions
  RefCountPtr<Epetra_Vector>    trueresidual_;

  /// Nonlinear iteration increment vector
  RefCountPtr<Epetra_Vector>    incvel_;

  /// accelerations at time n and n-1
  RefCountPtr<Epetra_Vector>    accn_;
  RefCountPtr<Epetra_Vector>    accnm_;

  /// velocities and pressures at time n+1, n and n-1
  RefCountPtr<Epetra_Vector>    velnp_;
  RefCountPtr<Epetra_Vector>    veln_;
  RefCountPtr<Epetra_Vector>    velnm_;

  RefCountPtr<Epetra_Vector>    dispnp_;
  RefCountPtr<Epetra_Vector>    dispn_;

  RefCountPtr<Epetra_Vector>    gridv_;

  // histvector --- a linear combination of velnm, veln (BDF)
  //                or veln, accn (One-Step-Theta)
  RefCountPtr<Epetra_Vector>    hist_;

  // solid displacement (might be moved somewhere else later)
  RefCountPtr<Epetra_Vector>    soliddispnp_;

  // Vectors used for convergence check
  // ----------------------------------
  RefCountPtr<Epetra_Map>       velrowmap_;
  RefCountPtr<Epetra_Map>       prerowmap_;

  RefCountPtr<Epetra_Vector>    onlyvel;
  RefCountPtr<Epetra_Vector>    onlypre;


  // time measurement
  // ----------------
  RefCountPtr<Time> timedyntot_     ;
  RefCountPtr<Time> timedyninit_    ;
  RefCountPtr<Time> timedynloop_    ;
  RefCountPtr<Time> timenlnloop_    ;
  RefCountPtr<Time> timeeleloop_    ;
  RefCountPtr<Time> timeapplydirich_;
  RefCountPtr<Time> timesolver_     ;

  RefCountPtr<TimeMonitor> tm0_ref_ ;
  RefCountPtr<TimeMonitor> tm1_ref_ ;
  RefCountPtr<TimeMonitor> tm2_ref_ ;
  RefCountPtr<TimeMonitor> tm3_ref_ ;
  RefCountPtr<TimeMonitor> tm4_ref_ ;
  RefCountPtr<TimeMonitor> tm5_ref_ ;
  RefCountPtr<TimeMonitor> tm6_ref_ ;

private:

/*!
\brief calculate lift&drag forces and angular momenta

*/
  void LiftDrag();


}; // class FluidImplicitTimeInt




#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
