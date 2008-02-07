/*!----------------------------------------------------------------------
\file fluid_genalpha_integration.H

\class FluidGenAlphaIntegration

\brief Generalised alpha time integration scheme for incompressible flow
problems.

Based on the method described in:

<pre>
  |    "A generalized-alpha method for integrating the filtered
  |     Navier-Stokes equations with a stabilized finite element method"
  |
  |    K.E. Jansen, C. Whiting and G.M. Hulbert
  |
  |    Computer Methods in Applied Mechanics and Engineering
</pre>

The idea for the implicit treatment of the pressure is taken from:

<pre>
  |    "Stabilized finite element methods for fluid dynamics using a
  |     hirarchical basis"
  |
  |    Christian H. Whiting
  |
  |    Dissertation (1999)
</pre>

In addition to what is described there, we do not evaluate the continuity
equation time n+alpha_F but at time n+1.
If the continuity equation is fulfilled exactly at time n, both variants
are equivalent --- but in the case of a randomly disturbed flow profile,
which might be used for turbulence calculations, the only fully implicit
method would provide consistency with the implicit treatment of the
pressure.


<pre>
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

*----------------------------------------------------------------------*/
#ifdef CCADISCRET
#ifndef FLUID_GENALPHA_INTEGRATION_H
#define FLUID_GENALPHA_INTEGRATION_H

#include "../drt_lib/drt_nodematchingoctree.H"
#include "../drt_lib/drt_periodicbc.H"
#include "turbulence_statistics.H"

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
#include "../drt_lib/linalg_mapextractor.H"
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
#ifdef RESULTTEST
  friend class FluidResultTest;
#endif

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
                           IO::DiscretizationWriter&        output);

   /*!
  \brief Destructor

  */
  virtual ~FluidGenAlphaIntegration();

  //! @name loops and iterations

  /*!
  \brief Do time integration from this step/time to endstep/endtime

  o calculate and set gamma parameter for generalised alpha method
    according to the second order necessary condition

  o Do the timeloop

  */
  void GenAlphaIntegrateTo(
    int                endstep,
    double             endtime
    );

  /*!
  \brief do nonlinear iteration, e.g. Newton like iteration

  In each step, the linearised system based on the (partially
  intermediate) velocities and accelerations is assembled and
  increments for the acceleration and pressure of the new timestep
  are calculated.
  The solutions (acceleration, velocity and pressure) of the new
  timestep "*np" are updated in each nonlinear iteration step using
  this increment.

  Convergence is checked by controlling residual- and increment- norms.

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
  \brief Write solution to file for IO and restart if we have a
  restart step.

  */
  void GenAlphaOutput();


  //@}


  //! @name service methods for nonlinear iteration


  /*!
  \brief Evaluate acceleration and velocity at the intermediate
  time level n+alpha_M and n+alpha_F

  Has to be called every time after the values at time level np have
  been updated!

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

  o Dirichlet conditions applied in here to be able to use the
    residual for the convergence check

  o the true residual is store in vector force_

  */
  void GenAlphaAssembleResidualAndMatrix(int          itnum);

  /*!
  \brief Solve linear problem.

  The counter itnum of the nonlinear iteration on input allows to
  do the solver initialisation only in the first step.

  */
  void GenAlphaCalcIncrement(
    int          itnum
    );

  /*!
  \brief update the current acceleration, velocity and pressure (*np)

  o loop all degrees of freedom

     -> if we have a velocity degree of freedom, update accelerations
        (in accnp) and velocities (in velnp)

     -> if we have a pressure dof, update pressure (in velnp)

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
    int whichinitialfield,
    int whichfunc
    );

  /*!
  \brief do a restart

  */
  void ReadRestart(int step);

  //@}


  //! @name verification (tests)

  /*!
  \brief error calculation for beltrami flow

  */
  void EvaluateErrorComparedToAnalyticalSol();

  //@}


  //! @name methods for turbulence models (dynamic Smagorinsky)

  void ApplyFilterForDynamicComputationOfCs();

   //@}

protected:

private:
  // input arguments
  RefCountPtr<DRT::Discretization> discret_;
  LINALG::Solver&                  solver_;
  ParameterList&                   params_;
  IO::DiscretizationWriter&        output_;

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

  // parameter for linearisation
  bool                             newton_;

  // write restart every uprestart_ steps --- restartstep_ is
  // the counter used for this purpose
  int                              restartstep_;
  int                              uprestart_;
  int                              writestep_;
  int                              upres_;

  // (fine-scale) subgrid-viscosity flag
  int                              fssgv_;

  // turbulence statistics
  RefCountPtr<TurbulenceStatistics> turbulencestatistics_;

  // system matrices: standard (stabilized) matrix and (fine-scale)
  // subgrid-viscosity matrix
  RefCountPtr<Epetra_CrsMatrix>    sysmat_;
  RefCountPtr<Epetra_CrsMatrix>    sysmat_sv_;

  // scale-separation matrix for VM3-based solution approach
  RefCountPtr<Epetra_CrsMatrix>    scalesep_;

  // rhs: standard (stabilized) residual vector (rhs for the incremental form),
  // (fine-scale) subgrid-viscosity residual and scaling vector
  RefCountPtr<Epetra_Vector>       residual_;
  RefCountPtr<Epetra_Vector>       residual_sv_;
  RefCountPtr<Epetra_Vector>       sugrvisc_;

  // the force vector is a copy of the residual without dirichlet
  // forces applied
  RefCountPtr<Epetra_Vector>       force_;

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

  LINALG::MapExtractor      velpressplitter_;

  // connectivity generated from periodic boundary conditions
  // --------------------------------------------------------
  //
  // (the connection between master and slave nodes on this proc)
  map<int,vector<int> > mapmastertoslave_;

  // Vectors used for filtering (for dynamic Smagorinsky model)
  // --------------------------
  RefCountPtr<Epetra_MultiVector>  filtered_vel_;
  RefCountPtr<Epetra_MultiVector>  filtered_reynoldsstress_;
  RefCountPtr<Epetra_MultiVector>  filtered_modeled_subgrid_stress_;
  RefCountPtr<Epetra_MultiVector>  col_filtered_vel_;
  RefCountPtr<Epetra_MultiVector>  col_filtered_reynoldsstress_;
  RefCountPtr<Epetra_MultiVector>  col_filtered_modeled_subgrid_stress_;

  RefCountPtr<vector<double> >     averaged_LijMij_;
  RefCountPtr<vector<double> >     averaged_MijMij_;
  RefCountPtr<vector<double> >     planecoords_;
  // time measurement
  // ----------------
  RefCountPtr<Time> timedyntot_     ;
  RefCountPtr<Time> timedyninit_    ;
  RefCountPtr<Time> timedynloop_    ;
  RefCountPtr<Time> timenlnloop_    ;
  RefCountPtr<Time> timeeleloop_    ;
  RefCountPtr<Time> timeevaldirich_ ;
  RefCountPtr<Time> timeapplydirich_;
  RefCountPtr<Time> timesolver_     ;
  RefCountPtr<Time> timeout_        ;
  RefCountPtr<Time> timenonlinup_   ;

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
  double                   dtfilter_;
};




#endif  // not FLUID_GENALPHA_INTEGRATION_H

#endif  // #ifdef CCADISCRET
