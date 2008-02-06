/*!----------------------------------------------------------------------
\file
\brief Associated with control routine for con.-diff.(in)stat. solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme 
       (with potential one-step-theta start algorithm)

     and stationary solver.

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET
#ifndef CONDIFIMPLICITINTEGRATION_H
#define CONDIFIMPLICITINTEGRATION_H

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
\class CondifImplicitTimeInt

\brief time integration for convection-diffusion problems

\author vg (vgravem@lnm.mw.tum.de)
*/
class CondifImplicitTimeInt
{

public:


  /*!
  \brief Standard Constructor

  */
  CondifImplicitTimeInt(RefCountPtr<DRT::Discretization> dis,
                        LINALG::Solver&                  solver,
                        ParameterList&                   params,
                        IO::DiscretizationWriter&        output);


  /*!
  \brief Destructor

  */
  virtual ~CondifImplicitTimeInt();


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

  void Solve(bool is_stat=false);

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
  \brief read restart data

  */
  void ReadRestart(int step);

protected:

  // don't want = operator and cctor
  CondifImplicitTimeInt operator = (const CondifImplicitTimeInt& old);
  CondifImplicitTimeInt (const CondifImplicitTimeInt& old);

  RefCountPtr<DRT::Discretization> discret_;
  LINALG::Solver&               solver_;
  ParameterList&                params_;
  IO::DiscretizationWriter&     output_;
  int                           myrank_;
  int                           maxentriesperrow_;

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

  int cdvel_;
  int fssgd_;

  /// system matrices: standard (stabilized) matrix and (fine-scale)
  /// subgrid-diffusivity matrix
  RefCountPtr<Epetra_CrsMatrix> sysmat_;
  RefCountPtr<Epetra_CrsMatrix> sysmat_sd_;

  /// scale-separation matrix for VM3-based solution approach
  RefCountPtr<Epetra_CrsMatrix> scalesep_;

  /// toggle vector indicating which dofs have Dirichlet BCs
  RefCountPtr<Epetra_Vector>    dirichtoggle_;

  RefCountPtr<Epetra_Vector>    invtoggle_;

  /// a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RefCountPtr<Epetra_Vector>    zeros_;

  /// the vector containing body and surface forces
  RefCountPtr<Epetra_Vector>    neumann_loads_;

  /// The residual vector
  RefCountPtr<Epetra_Vector>    residual_;
  RefCountPtr<Epetra_Vector>    sugrvisc_;

  /// time derivatives of phi at time n and n-1
  RefCountPtr<Epetra_Vector>    phidtn_;
  RefCountPtr<Epetra_Vector>    phidtnm_;
  /// phi at time n+1, n and n-1
  RefCountPtr<Epetra_Vector>    phinp_;
  RefCountPtr<Epetra_Vector>    phin_;
  RefCountPtr<Epetra_Vector>    phinm_;
  // histvector --- a linear combination of phinm, phin (BDF)
  //                or phin, phidtn (One-Step-Theta)
  RefCountPtr<Epetra_Vector>    hist_;


  // time measurement
  // ----------------
  RefCountPtr<Time> timedyntot_     ;
  RefCountPtr<Time> timedyninit_    ;
  RefCountPtr<Time> timedynloop_    ;
  RefCountPtr<Time> timeeleloop_    ;
  RefCountPtr<Time> timeapplydirich_;
  RefCountPtr<Time> timesolver_     ;

  RefCountPtr<TimeMonitor> tm0_ref_ ;
  RefCountPtr<TimeMonitor> tm1_ref_ ;
  RefCountPtr<TimeMonitor> tm2_ref_ ;
  RefCountPtr<TimeMonitor> tm3_ref_ ;
  RefCountPtr<TimeMonitor> tm4_ref_ ;
  RefCountPtr<TimeMonitor> tm5_ref_ ;

}; // class CondifImplicitTimeInt




#endif  // #ifndef CONDIFIMPLICITINTEGRATION_H
#endif  // #ifdef CCADISCRET
