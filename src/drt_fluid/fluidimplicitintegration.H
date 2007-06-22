/*!----------------------------------------------------------------------
\file
\brief Control routine for fluid time integration. Associated with
       fluidimplicitintegration.cpp


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef FLUIDIMPLICITINTEGRATION_H
#define FLUIDIMPLICITINTEGRATION_H

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

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>


extern "C" /* stuff which is c and is accessed from c++ */
{
#include "../headers/standardtypes.h"
#include "../fluid_full/fluid_prototypes.h"
}


using namespace std;
using namespace Teuchos;

/*!
\class FluidImplicitTimeInt

\brief time integration for fluid problems

\author gammi (gamnitzer@lnm.mw.tum.de)
*/
class FluidImplicitTimeInt
{
#ifdef RESULTTEST
  friend class FluidResultTest;
#endif

public:


  /*!
  \brief Standard Constructor

  */
  FluidImplicitTimeInt(RefCountPtr<DRT::Discretization> dis,
                       LINALG::Solver&       solver,
                       ParameterList&        params,
                       IO::DiscretizationWriter& output,
                       bool alefluid=false);


  /*!
  \brief Destructor

  */
  virtual ~FluidImplicitTimeInt();




  /*!
  \brief Create default parameter list

  */
  static void SetDefaults(ParameterList& params);

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
  void NonlinearSolve(bool is_stat=false);

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


  void ReadRestart(int step);

protected:

  // don't want = operator and cctor
  FluidImplicitTimeInt operator = (const FluidImplicitTimeInt& old);
  FluidImplicitTimeInt (const FluidImplicitTimeInt& old);

  RefCountPtr<DRT::Discretization> discret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  IO::DiscretizationWriter&         output_;
  int                           myrank_;
  int                           maxentriesperrow_;
  bool                          alefluid_;

  double time_;
  int step_;
  int stepmax_;
  double maxtime_;

  FLUID_TIMEINTTYPE timealgo_;

  int restartstep_;
  int uprestart_;
  int writestep_;
  int upres_;

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

}; // class FluidImplicitTimeInt




#endif  // #ifndef FLUIDIMPLICITINTEGRATION_H
#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
