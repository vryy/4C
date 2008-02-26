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
#include "drt_periodicbc.H"
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
#include "../drt_lib/linalg_systemmatrix.H"
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
  friend class FluidResultTest;

public:

  /*!

  \brief Standard Constructor (public)

  o trigger numbering of dofs

  o create empty system matrix and right hand side

  o allocate distributed solution vectors

  o generate maps to split velocity and pressure dofs

  */
  FluidGenAlphaIntegration(RCP<DRT::Discretization>  actdis,
                           LINALG::Solver&           solver,
                           ParameterList&            params,
                           IO::DiscretizationWriter& output);

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


  /*!
  \brief Do console output

  */

  void GenAlphaEchoToScreen(
    const string& what_to_print
    );


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
  void GenAlphaAssembleResidualAndMatrix();

  /*!
  \brief Solve linear problem.

  The counter itnum of the nonlinear iteration on input allows to
  do the solver initialisation only in the first step.

  */
  void GenAlphaCalcIncrement(const double nlnres);

  /*!
  \brief update the current acceleration, velocity and pressure (*np)

  o loop all degrees of freedom

     -> if we have a velocity degree of freedom, update accelerations
        (in accnp) and velocities (in velnp)

     -> if we have a pressure dof, update pressure (in velnp)

  */
  void GenAlphaNonlinearUpdate();

  /*!
  \brief check for convergence of nonlinear iteration

  */
  bool GenAlphaNonlinearConvergenceCheck(double& badestnlnnorm);

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
  RCP<DRT::Discretization> discret_;
  LINALG::Solver&                  solver_;
  ParameterList&                   params_;
  IO::DiscretizationWriter&        output_;

  // the dimension of the problem --- determines which dofs are pressure dofs
  int                              numdim_;

  // parallelism
  int                              myrank_;

  // max. simulation time and max. number of steps
  int                              endstep_;
  double                           endtime_;

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

  // parameters governing the nonlinear iteration
  int                              itenum_;
  int                              itemax_;
  double                           ittol_;

  // this variables are storing error norms
  double                           L2preresnorm_;
  double                           L2velresnorm_;
  double                           L2incvelnorm_;
  double                           L2incprenorm_;

  // skip the last element call (computation of residual) if
  // itemax_ is reached in order to get some speedup for turbulence
  // calculations
  bool                             skiplastelecall_;

  // write restart every uprestart_ steps --- restartstep_ is
  // the counter used for this purpose
  int                              restartstep_;
  int                              uprestart_;
  int                              writestep_;
  int                              upres_;

  // fine-scale subgrid-viscosity flag
  int                              fssgv_;

  // turbulence statistics
  RCP<TurbulenceStatistics> turbulencestatistics_;

  // (standard) system matrix
  RCP<LINALG::SparseMatrix>    sysmat_;

  // (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>       residual_;

  // the force vector is a copy of the residual without dirichlet
  // forces applied
  RCP<Epetra_Vector>       force_;

  // Nonlinear iteration increment vector contains accleeration
  // and pressure increment
  RCP<Epetra_Vector>       increment_;


  // toggle vector indicating which dofs have Dirichlet BCs
  RCP<Epetra_Vector>       dirichtoggle_;

  // a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>       zeros_;

  // the vector containing body and surface forces
  RCP<Epetra_Vector>       neumann_loads_;

  // acceleration and pressure at time n+1, n and intermediate time step n+alpha_M
  // in all acc vectors, every component associated to a pressure dof is
  // completely meaningless --- the Navier Stokes equations are a differential
  // algebraic equation, there's no evolution equation for the pressure!
  RCP<Epetra_Vector>       accnp_;
  RCP<Epetra_Vector>       accn_;
  RCP<Epetra_Vector>       accam_;

  // velocities and pressures at time n+1, n and intermediate time step n+alpha_F
  // in velaf_, every component associated to a pressure dof is completely
  // meaningless since the pressure is treated implicitly.
  RCP<Epetra_Vector>       velnp_;
  RCP<Epetra_Vector>       veln_;
  RCP<Epetra_Vector>       velaf_;

  // vm3 solver
  RCP<VM3_Solver>       vm3_solver_;

  // only necessary for vm3 solver: fine-scale solution vector
  RCP<Epetra_Vector>    fsvelaf_;


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
  RCP<Epetra_MultiVector>  filtered_vel_;
  RCP<Epetra_MultiVector>  filtered_reynoldsstress_;
  RCP<Epetra_MultiVector>  filtered_modeled_subgrid_stress_;
  RCP<Epetra_MultiVector>  col_filtered_vel_;
  RCP<Epetra_MultiVector>  col_filtered_reynoldsstress_;
  RCP<Epetra_MultiVector>  col_filtered_modeled_subgrid_stress_;

  RCP<vector<double> >     averaged_LijMij_;
  RCP<vector<double> >     averaged_MijMij_;
  RCP<vector<double> >     planecoords_;
  // time measurement
  // ----------------
  RCP<Time> timedyntot_         ;
  RCP<Time> timedyninit_        ;
  RCP<Time> timedynloop_        ;
  RCP<Time> timenlnloop_        ;
  RCP<Time> timeeleloop_        ;
  RCP<Time> timeevaldirich_     ;
  RCP<Time> timeapplydirich_    ;
  RCP<Time> timesolver_         ;
  RCP<Time> timeout_            ;
  RCP<Time> timenonlinup_       ;
  RCP<Time> timesparsitypattern_;


  RCP<TimeMonitor> tm0_ref_ ;
  RCP<TimeMonitor> tm1_ref_ ;
  RCP<TimeMonitor> tm2_ref_ ;
  RCP<TimeMonitor> tm3_ref_ ;
  RCP<TimeMonitor> tm4_ref_ ;
  RCP<TimeMonitor> tm5_ref_ ;
  RCP<TimeMonitor> tm6_ref_ ;
  RCP<TimeMonitor> tm7_ref_ ;
  RCP<TimeMonitor> tm8_ref_ ;
  RCP<TimeMonitor> tm9_ref_ ;

#ifdef PERF
  RefCountPtr<Time>        timeeleassemble_      ;
  RefCountPtr<Time>        timeelegetdoflocation_;

  RefCountPtr<Time>        timeelederxy2_    ;
  RefCountPtr<Time>        timeelederxy_     ;
  RefCountPtr<Time>        timeeletau_       ;
  RefCountPtr<Time>        timeelegalerkin_  ;
  RefCountPtr<Time>        timeelepspg_      ;
  RefCountPtr<Time>        timeelesupg_      ;
  RefCountPtr<Time>        timeelecstab_     ;
  RefCountPtr<Time>        timeelevstab_     ;
  RefCountPtr<Time>        timeelecrossrey_  ;
  RefCountPtr<Time>        timeeleintertogp_ ;
  RefCountPtr<Time>        timeeleseteledata_;
  RefCountPtr<Time>        timeeletdextras_  ;

  RefCountPtr<Time>        timeelegetvelnp_  ;
  RefCountPtr<Time>        timeeleinitsmag_  ;
  RefCountPtr<Time>        timeeleinitstab_  ;
  RefCountPtr<Time>        timeelesysmat_    ;
#endif

  double                   dtele_;
  double                   dtsolve_;
  double                   dtfilter_;


};



#endif  // not FLUID_GENALPHA_INTEGRATION_H

#endif  // #ifdef CCADISCRET
