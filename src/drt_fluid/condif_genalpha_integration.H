/*!----------------------------------------------------------------------
\file condif_genalpha_integration.H

\class CondifGenAlphaIntegration

\brief Associated with control routine for con.-diff. instationary solver
       based on generalized-alpha time-integration scheme

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef CONDIF_GENALPHA_INTEGRATION_H
#define CONDIF_GENALPHA_INTEGRATION_H

#include "../drt_lib/drt_nodematchingoctree.H"
#include "drt_periodicbc.H"

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
#include "../drt_lib/drt_function.H"
#include "../io/io_drt.H"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>

using namespace std;
using namespace Teuchos;


class CondifGenAlphaIntegration
{

public:

  /*!

  \brief Standard Constructor (public)

  o trigger numbering of dofs

  o create empty system matrix and right hand side

  o allocate distributed solution vectors

  */
  CondifGenAlphaIntegration(RefCountPtr<DRT::Discretization> actdis,
                            LINALG::Solver&                  solver,
                            ParameterList&                   params,
                            IO::DiscretizationWriter&        output);

   /*!
  \brief Destructor

  */
  virtual ~CondifGenAlphaIntegration();

  //! @name loops and iterations

  /*!
  \brief Do time integration from this step/time to endstep/endtime

  o calculate and set gamma parameter for generalized-alpha method
    according to the condition necessary for achieving second order

  o Do the timeloop

  */
  void GenAlphaIntegrateTo(
    int                endstep,
    double             endtime
    );

  /*!
  \brief do predictor-(1)corrector procedure

  In each step, the linear system based on the (partially intermediate)
  phi and dphi/dt is assembled and increments for the dphi/dt at the new
  timestep are calculated.
  The solutions (dphi/dt and phi) at the new timestep "*np" are updated
  using those increments.

  */
  void DoGenAlphaPredictorCorrector();


  //@}


  //! @name service functions for time loop

  /*!
  \brief Apply Dirichlet boundary conditions to phi vector.
  Apply surface Neumann conditions.

  */
  void GenAlphaApplyDirichletAndNeumann();

  /*!
  \brief calculate initial dphi/dt according to predicted phi and boundary
  values

  */
  void GenAlphaCalcInitialTimeDeriv();

  /*!
  \brief Update: current solution becomes old solution of next timestep.
  */
  void GenAlphaTimeUpdate();

  /*!
  \brief Write solution to file for IO and restart if we have a
  restart step.

  */
  void GenAlphaOutput();


  //@}


  //! @name service methods for predictor-(1)corrector procedure

  /*!
  \brief Evaluate dphi/dt and phi at the intermediate time levels
  n+alpha_M and n+alpha_F

  */
  void GenAlphaComputeIntermediateSol();

  /*!
  \brief Assemble systemmatrix and right hand side (residual).

  o Zeros out systemmatrix

  o Calls elements for assembly

  o Zeros out increment_ (solution vector)

  o Applies Dirichlet conditions to incremental linear system

  o Zeros out residual

  o Adds Neumann loads to the residual

  o the true residual is store in vector force_

  */
  void GenAlphaAssembleResidualAndMatrix();

  /*!
  \brief update the current dphi/dt and phi (*np)

  */
  void GenAlphaSolUpdate(
    );

  //@}

  //! @name initialization

  /*!
  \brief do a restart

  */
  void ReadRestart(int step);

  //@}


protected:

private:
  // input arguments
  RefCountPtr<DRT::Discretization> discret_;
  LINALG::Solver&                  solver_;
  ParameterList&                   params_;
  IO::DiscretizationWriter&            output_;

  // parallelism
  int                              myrank_;

  // matrix layout
  int                              maxentriesperrow_;

  // current time and step
  double                           time_;
  int                              step_;

  // current time step size
  double                           dt_;

  // parameters for generalised alpha
  double                           gamma_;
  double                           alphaM_;
  double                           alphaF_;

  // write restart every uprestart_ steps --- restartstep_ is
  // the counter used for this purpose
  int                              restartstep_;
  int                              uprestart_;
  int                              writestep_;
  int                              upres_;

  // system matrices: standard (stabilized) matrix and complete system matrix
  // including discontinuity-capturing (artificial diffusivity) matrix
  RefCountPtr<Epetra_CrsMatrix> sysmat_;
  RefCountPtr<Epetra_CrsMatrix> sysmat_dc_;

  // The residual vector --- more or less the rhs for the incremental
  // formulation!!!
  RefCountPtr<Epetra_Vector>       residual_;

  // the force vector is a copy of the residual without dirichlet
  // forces applied
  RefCountPtr<Epetra_Vector>       force_;

  // increment vector contains dphi/dt increment
  RefCountPtr<Epetra_Vector>       increment_;


  // toggle vector indicating which dofs have Dirichlet BCs
  RefCountPtr<Epetra_Vector>       dirichtoggle_;

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RefCountPtr<Epetra_Vector>       zeros_;

  // the vector containing body and surface forces
  RefCountPtr<Epetra_Vector>       neumann_loads_;

  // acceleration and pressure at time n+1, n and intermediate time step n+alpha_M
  // in all acc vectors, every component associated to a pressure dof is
  // completely meaningless --- the Navier Stokes equations are a differential
  // algebraic equation, there's no evolution equation for the pressure!
  RefCountPtr<Epetra_Vector>       phidtnp_;
  RefCountPtr<Epetra_Vector>       phidtn_;
  RefCountPtr<Epetra_Vector>       phidtam_;

  // velocities and pressures at time n+1, n and intermediate time step n+alpha_F
  // in velaf_, every component associated to a pressure dof is completely
  // meaningless since the pressure is treated implicitly.
  RefCountPtr<Epetra_Vector>       phinp_;
  RefCountPtr<Epetra_Vector>       phin_;
  RefCountPtr<Epetra_Vector>       phiaf_;


  // time measurement
  // ----------------
  RefCountPtr<Time> timedyntot_     ;
  RefCountPtr<Time> timedyninit_    ;
  RefCountPtr<Time> timedynloop_    ;
  RefCountPtr<Time> timepcloop_     ;
  RefCountPtr<Time> timeeleloop_    ;
  RefCountPtr<Time> timeevaldirich_ ;
  RefCountPtr<Time> timeapplydirich_;
  RefCountPtr<Time> timesolver_     ;
  RefCountPtr<Time> timeout_        ;
  RefCountPtr<Time> timesolupdate_  ;

  RefCountPtr<TimeMonitor> tm0_ref_ ;
  RefCountPtr<TimeMonitor> tm1_ref_ ;
  RefCountPtr<TimeMonitor> tm2_ref_ ;
  RefCountPtr<TimeMonitor> tm3_ref_ ;
  RefCountPtr<TimeMonitor> tm4_ref_ ;
  RefCountPtr<TimeMonitor> tm5_ref_ ;
  RefCountPtr<TimeMonitor> tm6_ref_ ;
  RefCountPtr<TimeMonitor> tm7_ref_ ;
  RefCountPtr<TimeMonitor> tm8_ref_ ;
  RefCountPtr<TimeMonitor> tm9_ref_ ;

  double                   dtele_;
  double                   dtsolve_;
};




#endif  // not CONDIF_GENALPHA_INTEGRATION_H

#endif  // #ifdef CCADISCRET
