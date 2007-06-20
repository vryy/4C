/*!----------------------------------------------------------------------
\file fluid_genalpha_integration.H

\brief Time integration according to dis. C. Whiting


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifdef TRILINOS_PACKAGE
#ifndef FLUID_GENALPHA_INTEGRATION_H
#define FLUID_GENALPHA_INTEGRATION_H

#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_periodicbc.H"

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

using namespace std;
using namespace Teuchos;


class FluidGenAlphaIntegration
{
public:

  /*!

  \brief Standard Constructor (public)

  o trigger numbering of dofs

  o create empty system matrix and right hand side

  o allocate distributed solution vectors

  o generate maps to split velocity and pressure dofs

  */
  FluidGenAlphaIntegration(RefCountPtr<DRT::Discretization> actdis,
                           LINALG::Solver&                  solver,
                           ParameterList&                   params,
                           IO::DiscretizationWriter&            output);

   /*!
  \brief Destructor

  */
  virtual ~FluidGenAlphaIntegration();

  //! @name loops and iterations

  /*!
  \brief Do time integration (time loop)

  */
  void GenAlphaIntegrateTo(
    int                endstep,
    double             endtime
    );

  /*!
  \brief do nonlinear iteration, e.g. Newton like iteration

  */
  void DoGenAlphaPredictorCorrectorIteration();



  //@}


  //! @name service functions for time loop

  /*!
  \brief Predict velocity and pressure of the new timestep. Up to now, we
  use a constant predictor for the velocity and the pressure.

  Remark: For Dirichlet nodes, no matter what was set here, velnp
          will be overwritten by the prescribed value. The
          accelerations are calculated after these Dirichlet values
          have been set.
  */
  void GenAlphaPredictNewSolutionValues();

  /*!
  \brief Apply Dirichlet boundary conditions to velocity vector.
  Apply surface Neumann conditions.

  */
  void GenAlphaApplyDirichletAndNeumann();

  /*!
  \brief calculate initial acceleration according to predicted
  velocities and boundary values

  */
  void GenAlphaCalcInitialAccelerations();



  /*!
  \brief Update the solution after convergence of the nonlinear
         iteration. Current solution becomes old solution of next
         timestep.
  */
  void GenAlphaTimeUpdate();

  /*!
  \brief Write solution to file

  */
  void GenAlphaOutput();


  //@}


  //! @name service methods for nonlinear iteration


  /*!
  \brief Evaluate acceleration and velocity at the intermediate
  time level n+alpha_M and n+alpha_F

  */
  void GenAlphaComputeIntermediateSol();

  /*!
  \brief Assemble linearised system and solve linear problem

  */
  void GenAlphaCalcIncrement(
    int          itnum
    );

  /*!
  \brief update the current acceleration, velocity and pressure

  */
  void GenAlphaNonlinearUpdate(
    );

  /*!
  \brief check for convergence of nonlinear iteration

  */
  bool GenAlphaNonlinearConvergenceCheck(
    int          itnum
    );

  //@}

  //! @name initialisation

  /*!
  \brief set an initial flow field

  */

  void SetInitialFlowField(
    int whichinitialfield
    );

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

  // the dimension of the problem --- determines which dofs are pressure dofs
  int                              numdim_;

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

  // system matrix --- stiffness and mass are assembled in one matrix!
  RefCountPtr<Epetra_CrsMatrix>    sysmat_;

  // The residual vector --- more or less the rhs for the incremental
  // formulation!!!
  RefCountPtr<Epetra_Vector>       residual_;

  // Nonlinear iteration increment vector contains accleeration
  // and pressure increment
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
  RefCountPtr<Epetra_Vector>       accnp_;
  RefCountPtr<Epetra_Vector>       accn_;
  RefCountPtr<Epetra_Vector>       accam_;

  // velocities and pressures at time n+1, n and intermediate time step n+alpha_F
  // in velaf_, every component associated to a pressure dof is completely
  // meaningless since the pressure is treated implicitly.
  RefCountPtr<Epetra_Vector>       velnp_;
  RefCountPtr<Epetra_Vector>       veln_;
  RefCountPtr<Epetra_Vector>       velaf_;


  // Vectors used for convergence check
  // ----------------------------------
  RefCountPtr<Epetra_Map>          velrowmap_;
  RefCountPtr<Epetra_Map>          prerowmap_;

  RefCountPtr<Epetra_Vector>       onlyvel;
  RefCountPtr<Epetra_Vector>       onlypre;


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
};




#endif  // not FLUID_GENALPHA_INTEGRATION_H

#endif  // #ifdef TRILINOS_PACKAGE
#endif  // #ifdef CCADISCRET
