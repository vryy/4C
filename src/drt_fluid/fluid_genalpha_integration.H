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
which might be used for turbulence calculations, only the fully implicit
method would provide consistency with the implicit treatment of the
pressure.

The method includes additional features required by the FSI adapter
class.


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

  o trigger numbering of dofs (including periodic boundaries)

  o create empty system matrix and right hand side

  o allocate distributed solution vectors

  o generate maps to split velocity and pressure dofs

  o set alpha_F and alpha_M

  o calculate and set gamma parameter for generalised alpha method
    according to the second order necessary condition

  */
  FluidGenAlphaIntegration(
    RCP<DRT::Discretization>  actdis,
    LINALG::Solver&           solver,
    ParameterList&            params,
    IO::DiscretizationWriter& output,
    bool                      alefluid);

   /*!
  \brief Destructor

  */
  virtual ~FluidGenAlphaIntegration();

  //! @name loops and iterations

  /*!
  \brief Do time integration from this step/time to endstep/endtime

  o Do the timeloop

  o call summary for time monitor

  */
  void GenAlphaTimeloop();

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
  In special cases, the residual check in the last iteration could be
  suppressed

  */
  void DoGenAlphaPredictorCorrectorIteration();


  //@}


  //! @name service functions for time loop


  /*!
  \brief Increase time step counter and time
  */
  void GenAlphaIncreaseTimeAndStep()
    {
      step_++;
      time_+=dt_;

      return;
    }

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

         Update ALE displacements in case of ale

         In case of time dependent subscales, the time update
         includes an element call.
  */
  void GenAlphaTimeUpdate();


  /*!
  \brief Include this result into the averaging procedure for
  turbulence statistics

  */
  void GenAlphaTakeSample();


  /*!
  \brief Write solution to file for IO and restart if we have a
  restart step.

  */
  void GenAlphaOutput();


  /*!
  \brief Do console output, what is specified by input string

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

  o For a dynamic Smagorinsky model, the velocity field is
    filtered using a box filter on element poatches adjacent
    to the fluid nodes. The filtered fields are used to
    compute the eddy viscosity --- hex8 only!

  o Zeros out systemmatrix

  o Calls elements for assembly

  o Zeros out increment_ (solution vector)

  o Applies Dirichlet conditions to incremental linear system

  o Zeros out residual

  o Adds Neumann loads to the residual

  o Dirichlet conditions applied in here to be able to use the
    residual for the convergence check

  o the true residual is stored in vector force_ --- used for
    Dirichlet load computations

  */
  void GenAlphaAssembleResidualAndMatrix();

  /*!
  \brief Solve linear problem.

  Enforces to set up new a preconditioner in special cases.

  Includes dynamic update of tolerances for linear solver according
  to Michaels ad-hoc scheme.

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

  badestnlnnorm is used for the dynamic setting of the linear solve TOL

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
  /*!
    o Apply the box filter to compute averaged velocities and
      reynoldsstresses (this is an element call!)

    o Use filtered quantities to compute LijMij, MijMij and
      thus Cs

    o For turbulent channel flows, averaging of LijMij and
      MijMij is performed in all homogeneous planes.

    o negative values of Cs are clipped

  */
  void ApplyFilterForDynamicComputationOfCs();

   //@}

  //! @name access methods for composite algorithms
  /// monolithic FSI needs to access the linear fluid problem

  /*!
  \brief return the current increment for an initial guess

  */
  Teuchos::RCP<Epetra_Vector> InitialGuess() { return increment_; }

  /*!
  \brief return the residual vector, i.e. the rhs of the
  incremental fluid problem

  */
  Teuchos::RCP<Epetra_Vector> Residual()     { return residual_; }

  /*!
  \brief return the forces on Dirichlet nodes --- correct
  density scaling has tobe provided by the fluid implementation

  */
  Teuchos::RCP<Epetra_Vector> TrueResidual() { return force_; }

  /*!
  \brief return current velocity and pressure solutions

  */
  Teuchos::RCP<Epetra_Vector> Velnp() { return velnp_; }

  /*!
  \brief return previous velocity and pressure solutions

  */
  Teuchos::RCP<Epetra_Vector> Veln() { return veln_; }

  /*!
  \brief return pointer to gridvelocities at intermediate stage.
  They are set from outside (strange, isn't it?)

  */
  Teuchos::RCP<Epetra_Vector> GridVel() { return gridvelaf_; }

  /*!
  \brief Allow access to the vector marking the Dirichlet dofs

  */
  Teuchos::RCP<Epetra_Vector> Dirichlet() { return dirichtoggle_; }

  /*!
  \brief the inverse Dirich toggle vector for the FSI method

  */
  Teuchos::RCP<Epetra_Vector> InvDirichlet() { return invdirtoggle_; }

  /*!
  \brief Allow access to current grid displacements

  */
  Teuchos::RCP<Epetra_Vector> Dispnp() { return dispnp_; }

  /*!
  \brief Allow access to previous grid displacements

  */
  Teuchos::RCP<Epetra_Vector> Dispn() { return dispn_; }

  /*!
  \brief return the system matrix

  */
  Teuchos::RCP<LINALG::SparseMatrix> SysMat() { return sysmat_; }

  /*!
  \brief return the fluid discretisation

  */
  Teuchos::RCP<DRT::Discretization> Discretization()
    { return discret_; }

  /*!
  \brief provide the dofrowmap of the velocity dofs

  */
  Teuchos::RCP<const Epetra_Map> VelocityRowMap()
    { return velpressplitter_.OtherMap(); }

  /*!
  \brief provide the dofrowmap of the pressure dofs

  */
  Teuchos::RCP<const Epetra_Map> PressureRowMap()
    { return velpressplitter_.CondMap(); }

  /*
  \brief use increment from external call to update velnp
  (monolithic FSI)
  */
  void ExternIncrementOfVelnp(Teuchos::RCP<const Epetra_Vector> inc)
    {
      // copy inc from external Newton to increment_
      increment_->Update(1.0, *inc, 0.);

      // update (u,p) trial
      velnp_->Update(1.0,*increment_,1.0);
      return;
    }

  /*
  \brief get scaling factor for residual --- just the density
  for this generalized alpha implementation
  */
  double ResidualScaling() const { return density_; }
  /*
  \brief return current time
  */
  double Time() const { return time_; }
  /*
  \brief return current time step number
  */
  int Step() const { return step_; }
  /*
  \brief return the nonlinear itemax
  */
  int Itemax() const { return params_.get<int>("max nonlin iter steps"); }
  /*
  \brief set the nonlinear itemax
  */
  void SetItemax(int itemax) { params_.set<int>("max nonlin iter steps", itemax); }

  /*!
  \brief determine grid velocity according to genalpha time integration

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
    This is the linear solve as required for steepest descent FSI.
   */
  void LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);


protected:

private:
  //! @name input arguments of the constructor
  //
  //! the discretization
  RCP<DRT::Discretization>  discret_;
  //! solver object
  LINALG::Solver&           solver_ ;
  //! parameterlist including time params, stabilization params and turbulence sublist
  ParameterList&            params_ ;
  //! io writer/reader
  IO::DiscretizationWriter& output_ ;
  //! decides whether we use an Eulerian or an ALE formulation
  bool                      alefluid_;
  //@}

  //! the dimension of the problem --- determines which dofs
  //are pressure dofs
  int                       numdim_ ;

  //! parallelism --- my processor id
  int                       myrank_ ;

  //! fluid density
  double                    density_;

  //! @name time integration parameters
  //
  //! max. number of time steps
  int                       endstep_;
  //! max. simulation time
  double                    endtime_;
  //! current time
  double                    time_;
  //! current time step
  int                       step_;
  //! current time step size
  double                    dt_;

  //! parameter for generalised alpha --- alpha "momentum"
  double                    alphaM_;
  //! parameter for generalised alpha --- alpha "force"
  double                    alphaF_;
  //! parameter for generalised alpha --- computed from the latter
  //in order to achieve 2nd order accuracy
  double                    gamma_;
  //@}

  //! @name parameters for nonlinear iteration
  //
  //! parameter for linearisation
  bool                      newton_;

  //! parameters governing the nonlinear iteration --- iteration counter
  int                       itenum_      ;
  //! parameters governing the nonlinear iteration --- upper bound for
  //! number of nonlinear iterations
  int                       itemax_      ;
  //! parameters governing the nonlinear iteration --- convergence criterion
  double                    ittol_       ;

  // these variables are storing error norms
  //! L2 error norm --- residual associated with pressure dofs, scaled by density
  double                    L2preresnorm_;
  //! L2 error norm --- residual associated with acceleration dofs, scaled by density
  double                    L2velresnorm_;
  //! L2 error norm --- normalized increment in iteration step (velocity dofs)
  double                    L2incvelnorm_;
  //! L2 error norm --- normalized increment in iteration step (pressure dofs)
  double                    L2incprenorm_;

  //! skip the last element call (computation of residual) if
  //!  itemax_ is reached in order to get some speedup for turbulence
  //!  calculations
  bool                      skiplastelecall_;
  //@}

  //! @name IO features

  // write restart every uprestart_ steps --- restartstep_ is
  // the counter used for this purpose

  //! write restart every uprestart_ steps
  int                       uprestart_;
  //! count steps till uprestart_
  int                       restartstep_;
  //! write result every upres_ steps
  int                       upres_;
  //! count steps till upres_
  int                       writestep_;

  //! parameters for sampling/dumping period
  int samstart_;
  int samstop_;
  int dumperiod_;

  //! flag for special flow: currently channel flow or flow in a lid-driven cavity
  string special_flow_;

  //! turbulence statistics for turbulent channel flow example
  RCP<TurbulenceStatistics> turbulencestatistics_;

  //! turbulence statistics for lid-driven cavity example
  RCP<TurbulenceStatisticsLdc> turbulencestatistics_ldc_;

  //@}


  //! @name vectors stating the linear system
  //        ---------------------------------
  //
  //! (standard) system matrix
  RCP<LINALG::SparseMatrix> sysmat_;

  //! (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>        residual_;

  //! toggle vector indicating which dofs have Dirichlet BCs
  RCP<Epetra_Vector>        dirichtoggle_;

  //! and its inverse
  RCP<Epetra_Vector>        invdirtoggle_;

  //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>        zeros_;

  //! the vector containing body and surface forces
  RCP<Epetra_Vector>        neumann_loads_;

  //! Nonlinear iteration increment vector contains accleeration
  //! and pressure increment
  RCP<Epetra_Vector>        increment_;
  //@}


  //! @name Force computation on Dirichlet boundaries
  // -----------------------------------------
  //
  //! the force vector is a copy of the residual without dirichlet
  //! conditions applied --- probably, for all inner nodes,
  //! "equilibrium" would be a better name - but for the Dirichlet
  //! nodes, this vector contains the consistent (normed) force on
  //! boundary nodes
  RCP<Epetra_Vector>        force_;
  //@}


  //! @name vectors for the time integration process
  //        ----------------------------------------
  //
  // acceleration and pressure at time n+1, n and intermediate time step n+alpha_M
  // in all acc vectors, every component associated to a pressure dof is
  // completely meaningless --- the Navier Stokes equations are a differential
  // algebraic equation, there's no evolution equation for the pressure!

  //! acceleration at time n+1
  //! in accnp_, every component associated to a pressure dof is
  //! completely meaningless --- the Navier Stokes equations are a differential
  //! algebraic equation, there's no evolution equation for the pressure!
  RCP<Epetra_Vector>        accnp_;
  //! acceleration at time n
  //! in  accn_, every component associated to a pressure dof is
  //! completely meaningless --- the Navier Stokes equations are a differential
  //! algebraic equation, there's no evolution equation for the pressure!
  RCP<Epetra_Vector>        accn_;
  //! acceleration at intermediate time step n+alpha_M
  //! in accam_, every component associated to a pressure dof is
  //! completely meaningless --- the Navier Stokes equations are a differential
  //! algebraic equation, there's no evolution equation for the pressure!
  RCP<Epetra_Vector>        accam_;

  //! velocities and pressures at time n+1
  RCP<Epetra_Vector>        velnp_;
  //! velocities and pressures at time n
  RCP<Epetra_Vector>        veln_;
  //! velocities at intermediate time step n+alpha_F
  //! in velaf_, every component associated to a pressure dof is completely
  //! meaningless since the pressure is treated implicitly.
  RCP<Epetra_Vector>        velaf_;

  //! ALE related stuff: grid displacements, time n+1
  RCP<Epetra_Vector>        dispnp_   ;
  //! ALE related stuff: most recent grid displacements, time n
  RCP<Epetra_Vector>        dispn_    ;
  //! ALE related stuff: grid velocities, most recent time step n
  RCP<Epetra_Vector>        gridveln_ ;
  //! ALE related stuff: intermediate grid velocities at time n+af
  RCP<Epetra_Vector>        gridvelaf_;
  //@}

  //! @name vm3 solver
  //        ----------
  //
  //! vm3 solver
  RCP<VM3_Solver>           vm3_solver_;

  //! coarse- and fine-scale solution vector
  RCP<Epetra_Vector>        csvelaf_;
  RCP<Epetra_Vector>        fsvelaf_;

  //! fine-scale subgrid-viscosity flag
  string                    fssgv_;

  //@}

  //! @name Vectors used for convergence check
  //         ----------------------------------

  //! extracts velocities and pressure dofs from dofmap
  LINALG::MapExtractor      velpressplitter_;
  //@}


  //! @name connectivity generated from periodic boundary conditions
  //         --------------------------------------------------------
  //
  //! the connection between master and slave nodes on this proc
  map<int,vector<int> >     mapmastertoslave_;
  //@}

  //! @name vectors used for filtering (for dynamic Smagorinsky model)
  //        --------------------------

  //! the box filtered velocities in nodes (3 vectors)
  RCP<Epetra_MultiVector>   filtered_vel_;
  //! the box filtered reynoldsstresses in nodes (9 vectors)
  RCP<Epetra_MultiVector>   filtered_reynoldsstress_;
  //! the modeled subgrid stress in nodes (9 vectors)
  RCP<Epetra_MultiVector>   filtered_modeled_subgrid_stress_;
  //! the filtered vel exported to column map
  RCP<Epetra_MultiVector>   col_filtered_vel_;
  //! the filtered reystress exported to column map
  RCP<Epetra_MultiVector>   col_filtered_reynoldsstress_;
  //! the modeled subgrid stresses exported to column map
  RCP<Epetra_MultiVector>   col_filtered_modeled_subgrid_stress_;

  //! for turbulent channel flow, LijMij in averaged in each hom. plane
  RCP<vector<double> >      averaged_LijMij_;
  //! for turbulent channel flow, MijMij in averaged in each hom. plane
  RCP<vector<double> >      averaged_MijMij_;
  //! the plane coordinates for the above mentioned averaging procedure
  RCP<vector<double> >      planecoords_;
  //@}

  //! @name time measurement
  //        ----------------

  //! time measurement --- total time
  RCP<Teuchos::Time>                 timedyntot_         ;
  //! time measurement --- time in constructor
  RCP<Teuchos::Time>                 timedyninit_        ;
  //! time measurement --- time spent in timeloop
  RCP<Teuchos::Time>                 timedynloop_        ;
  //! time measurement --- time spent in nonlinear iteration
  RCP<Teuchos::Time>                 timenlnloop_        ;
  //! time measurement --- time spent on the element loop
  RCP<Teuchos::Time>                 timeeleloop_        ;
  //! time measurement --- time spent to evaluate dirichlet conditions
  RCP<Teuchos::Time>                 timeevaldirich_     ;
  //! time measurement --- time spent to apply the dirichlet conditions
  //! to the system, i.e. to zero out rows and write ones one the
  //! diagonal elements
  RCP<Teuchos::Time>                 timeapplydirich_    ;
  //! time to solve the linear problem
  RCP<Teuchos::Time>                 timesolver_         ;
  //! time required for output
  RCP<Teuchos::Time>                 timeout_            ;
  //! time for the nonlinear vector update
  RCP<Teuchos::Time>                 timenonlinup_       ;
  //! time for the call to sysmat_->complete
  RCP<Teuchos::Time>                 timesparsitypattern_;

  //! used to add time to timedyntot
  RCP<Teuchos::TimeMonitor>          tm0_ref_ ;
  //! used to add time to timeevaldirich
  RCP<Teuchos::TimeMonitor>          tm1_ref_ ;
  //! used to add time to timedynloop
  RCP<Teuchos::TimeMonitor>          tm2_ref_ ;
  //! used to add time to timeeleloop
  RCP<Teuchos::TimeMonitor>          tm3_ref_ ;
  //! used to add time to timeapplydirich
  RCP<Teuchos::TimeMonitor>          tm4_ref_ ;
  //! used to add time to timesolver
  RCP<Teuchos::TimeMonitor>          tm5_ref_ ;
  //! used to add time to timenlnloop
  RCP<Teuchos::TimeMonitor>          tm6_ref_ ;
  //! used to add time to timedyninit
  RCP<Teuchos::TimeMonitor>          tm7_ref_ ;
  //! used to add time to timeout
  RCP<Teuchos::TimeMonitor>          tm8_ref_ ;
  //! used to add time to timenonlinup
  RCP<Teuchos::TimeMonitor>          tm9_ref_ ;

#ifdef PERF
  // very detailed time measurement
  RefCountPtr<Teuchos::Time>         timeeleassemble_      ;
  RefCountPtr<Teuchos::Time>         timeelegetdoflocation_;

  RefCountPtr<Teuchos::Time>         timeelederxy2_    ;
  RefCountPtr<Teuchos::Time>         timeelederxy_     ;
  RefCountPtr<Teuchos::Time>         timeeletau_       ;
  RefCountPtr<Teuchos::Time>         timeelegalerkin_  ;
  RefCountPtr<Teuchos::Time>         timeelepspg_      ;
  RefCountPtr<Teuchos::Time>         timeelesupg_      ;
  RefCountPtr<Teuchos::Time>         timeelecstab_     ;
  RefCountPtr<Teuchos::Time>         timeelevstab_     ;
  RefCountPtr<Teuchos::Time>         timeelecrossrey_  ;
  RefCountPtr<Teuchos::Time>         timeeleintertogp_ ;
  RefCountPtr<Teuchos::Time>         timeeleseteledata_;
  RefCountPtr<Teuchos::Time>         timeeletdextras_  ;

  RefCountPtr<Teuchos::Time>         timeelegetvelnp_  ;
  RefCountPtr<Teuchos::Time>         timeeleinitsmag_  ;
  RefCountPtr<Teuchos::Time>         timeeleinitstab_  ;
  RefCountPtr<Teuchos::Time>         timeelesysmat_    ;
#endif

  //! count cpu time for element call on proc
  double                    dtele_;
  //! count cpu time for solver call on proc
  double                    dtsolve_;
  //! count cpu time for filtering (in case of a dynamic Smagorinsky model) on proc
  double                    dtfilter_;

  //@}

};

#endif  // not FLUID_GENALPHA_INTEGRATION_H

#endif  // #ifdef CCADISCRET
