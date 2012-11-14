/*----------------------------------------------------------------------*/
/*!
\file fluid_genalpha_integration.H

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
  |     hierarchical basis"
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
*/
/*----------------------------------------------------------------------*/
#ifndef FLUID_GENALPHA_INTEGRATION_H
#define FLUID_GENALPHA_INTEGRATION_H

#include "../drt_fluid/drt_periodicbc.H"
#include "../drt_fluid/drt_transfer_turb_inflow.H"
#include "../drt_fluid/turbulence_statistic_manager.H"
#include "../drt_fluid/dyn_smag.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../drt_fluid/fluid_timint.H"
#include "../drt_fluid/fluidresulttest.H"

#include "../drt_io/io.H"

#include "../drt_fluid/fluid_MHD_evaluate.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_function.H"
#include "../drt_lib/drt_locsys.H"
#include "../drt_lib/drt_nodematchingoctree.H"
#include "../linalg/linalg_solver.H"
#include "../linalg/linalg_utils.H"
#include "../linalg/linalg_mapextractor.H"
#include "../linalg/linalg_sparsematrix.H"
#include "../drt_inpar/inpar_fluid.H"
#include "../drt_inpar/inpar_poroelast.H"

#ifdef PARALLEL
#include "Epetra_MpiComm.h"
#else
#include "Epetra_SerialComm.h"
#endif /* PARALLEL */

#include "Teuchos_RefCountPtr.hpp"
#include "Teuchos_ParameterList.hpp"

#include "Epetra_Vector.h"
#include "Epetra_Operator.h"
#include "Epetra_RowMatrix.h"

#include <ml_common.h>
#include <ml_include.h>
#include <ml_epetra_utils.h>
#include <ml_epetra.h>
#include <ml_epetra_operator.h>
#include <ml_MultiLevelPreconditioner.h>
#include <ml_agg_genP.h>
#include <ml_operator.h>
#include <MLAPI_Error.h>
#include <MLAPI_CompObject.h>
#include <MLAPI_TimeObject.h>
#include <MLAPI_Operator.h>
#include <MLAPI_Operator_Utils.h>
#include <MLAPI_MultiVector.h>
#include <MLAPI_InverseOperator.h>
#include <MLAPI_Expressions.h>
#include <MLAPI_BaseOperator.h>
#include <MLAPI_Workspace.h>
#include <MLAPI_Aggregation.h>
#include <MLAPI_Eig.h>

#include <ctime>
#include <cstdlib>
#include <iostream>
#include <Teuchos_TimeMonitor.hpp>

using namespace std;
using namespace Teuchos;
using namespace MLAPI;

namespace FLD {

// forwarded declarations
class TurbulenceStatisticManager;

class FluidGenAlphaIntegration : public TimInt
{
  //! @name Friends
  friend class FluidResultTest;
  friend class TurbulenceStatisticManager;

  //@}
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
    const Teuchos::RCP<DRT::Discretization>&    actdis,
    const Teuchos::RCP<LINALG::Solver>&             solver,
    const Teuchos::RCP<ParameterList>&              params,
    const Teuchos::RCP<IO::DiscretizationWriter>&   output,
    bool                        alefluid,
    RCP<map<int,vector<int> > > pbcmapmastertoslave);

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
  \brief preparation of time step by performing several procedures

  */
  void GenAlphaPrepareTimeStep();

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
  \brief Predict velocity and pressure of the new timestep. Choices for
   the predictor are a constant predictor for the velocity and the pressure
   or one of the following: zero-acceleration predictor, constant-
   acceleration predictor, constant-increment predictor.

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
  \brief statistics time sample and output of solution and statistics

  */
  void GenAlphaStatisticsAndOutput();


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
    const INPAR::FLUID::InitialField initfield,
    const int whichfunc
    );

  /// set initial porosity field
  void SetInitialPorosityField(const INPAR::POROELAST::InitialField initfield,const int startfuncno)
  { dserror("not implemented!");return; };

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
  Teuchos::RCP<const Epetra_Vector> InitialGuess() { return increment_; }

  /*!
  \brief return the residual vector, i.e. the rhs of the
  incremental fluid problem

  */
  Teuchos::RCP<Epetra_Vector> Residual()     { return residual_; }

  /*!
  \brief return the forces on Dirichlet nodes --- correct
  density scaling has to be provided by the fluid implementation

  */
  Teuchos::RCP<const Epetra_Vector> TrueResidual() { return force_; }

  /*!
  \brief return current velocity and pressure solutions

  */
  Teuchos::RCP<const Epetra_Vector> Velnp() { return velnp_; }

  /*!
  \brief return velocity and pressure solutions at
         intermediate time step n+alpha_F

         We have to ensure, that the current velnp value is used for
         determination of an intermediate velocity before it is provided to
         the outer world
  */
  Teuchos::RCP<const Epetra_Vector> Velaf()
    {
      // to be sure: (re-)compute current value of velaf_
      GenAlphaComputeIntermediateSol();
      return velaf_;
    }

  /*!
  \brief return acceleration values at
         intermediate time step n+alpha_M
  */
  Teuchos::RCP<const Epetra_Vector> Accam()
    {
      // to be sure: (re-)compute current value of accam_
      // assuming that accnp_ is already up to date!
      GenAlphaComputeIntermediateSol();
      return accam_;
    }

  /*!
  \brief return previous velocity and pressure solutions

  */
  Teuchos::RCP<const Epetra_Vector> Veln() { return veln_; }

  /*!
  \brief return pointer to gridvelocities at intermediate stage.
  They are set from outside (strange, isn't it?)

  */
  Teuchos::RCP<const Epetra_Vector> GridVel() { return gridvelaf_; }

  /*!
  \brief provide access to the Dirichlet map

  */
  const Teuchos::RCP<const LINALG::MapExtractor> DirichMaps() { return dbcmaps_; }

  /*!
  \brief Expand the Dirichlet DOF set

  The method expands the DOF set (map) which contains the DOFs
  subjected to Dirichlet boundary conditions. For instance, the method is
  called by the staggered FSI in which the velocities on the FSI
  interface are prescribed by the other fields.

  */
  void AddDirichCond(const Teuchos::RCP<const Epetra_Map> maptoadd);

  /*!
  /brief Shorten the Dirichlet DOF set

  The method shortens the DOF set (map) which contains the DOFs
  subjected to Dirichlet boundary conditions. For instance, the method is
  called by the monolithic FSI where fluid boundary conditions on the
  FSI interface have to be removed if it is a fluid split.

  */
  virtual void RemoveDirichCond(const Teuchos::RCP<const Epetra_Map> maptoremove);

  /*!
  \brief Extract the Dirichlet toggle vector based on Dirichlet BC maps

  This method provides backward compatability only. Formerly, the Dirichlet conditions
  were handled with the Dirichlet toggle vector. Now, they are stored and applied
  with maps, ie #dbcmaps_. Eventually, this method will be removed.

  */
  const Teuchos::RCP<const Epetra_Vector> Dirichlet();

  /*!
  \brief Extract the Inverse Dirichlet toggle vector based on Dirichlet BC maps

  This method provides backward compatability only. Formerly, the Dirichlet conditions
  were handled with the Dirichlet toggle vector. Now, they are stored and applied
  with maps, ie #dbcmaps_. Eventually, this method will be removed.

  */
  const Teuchos::RCP<const Epetra_Vector> InvDirichlet();

  /*!
  \brief Allow access to current grid displacements

  */
  Teuchos::RCP<const Epetra_Vector> Dispnp() { return dispnp_; }

  /*!
  \brief Allow access to previous grid displacements

  */
  Teuchos::RCP<const Epetra_Vector> Dispn() { return dispn_; }

  /*!
  \brief return the system matrix if it is a plain matrix

  */
  Teuchos::RCP<LINALG::SparseMatrix> SystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::SparseMatrix>(sysmat_); }

  /*!
  \brief return the system matrix if it is a block matrix

  */
  Teuchos::RCP<LINALG::BlockSparseMatrixBase> BlockSystemMatrix() { return Teuchos::rcp_dynamic_cast<LINALG::BlockSparseMatrixBase>(sysmat_); }


  Teuchos::RCP<LINALG::BlockSparseMatrixBase> ShapeDerivatives() { return shapederivatives_; }


  /*!
  \brief provide the velocity-pressure splitter

  */
  LINALG::MapExtractor& VelPresSplitter() {return velpressplitter_;}

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

  /*!
  \brief use increment from external call to update velnp
  (monolithic FSI)
  */
  void ExternIncrementOfVelnp(Teuchos::RCP<const Epetra_Vector> dacc)
    {
      // copy dacc from external Newton to increment_
      increment_->Update(1.0, *dacc, 0.);

      // reset velocity, pressure and acceleration to t(n) values because dacc
      // contains the whole Newton step:
      //
      //        +-   n+1       n -+
      // dacc = | acc     - acc   |
      //        +-   (i+1)       -+
      //
      // along with the pressure values.
      //
      // dacc really is nothing but the sum of all solutions of the linear
      // problems of one time step.

      // Do not destroy Dirichlet values!
      dbcmaps_->InsertOtherVector(dbcmaps_->ExtractOtherVector(veln_), velnp_);

      GenAlphaCalcInitialAccelerations();

      // normal Newton iteration update
      //
      // Note that here the increment_ is the whole Newton increment rather
      // than the latest increment contribution. However, the update remains
      // the same.
      GenAlphaNonlinearUpdate();
      return;
    }

  /*!
  \brief get scaling factor for residual --- just the density
  for this generalized alpha implementation
  */
  double ResidualScaling() const { return density_; }

  /// return time integration factor
  double TimIntParam() const;

  /*!
  \brief return current time
  */
  double Time() const { return time_; }
  /*!
  \brief return current time step number
  */
  int Step() const { return step_; }
  /// time step
  double Dt() const { return dt_; }
  /// time integration parameter
  double Gamma() const { return gamma_; }
  /*!
  \brief return the nonlinear itemax
  */
  int Itemax() const { return params_->get<int>("max nonlin iter steps"); }
  /*!
  \brief set the nonlinear itemax
  */
  void SetItemax(int itemax) { params_->set<int>("max nonlin iter steps", itemax); }

  /*!
  \brief return type of time integration scheme

  */
  const INPAR::FLUID::TimeIntegrationScheme TimIntScheme() const {return INPAR::FLUID::timeint_gen_alpha;}

  /*!
  \brief determine grid velocity according to genalpha time integration

  */
  void UpdateGridv();

  /*!
  \brief prepare AVM3-based scale separation

  */
  void AVM3Preparation();

  /*!
  \brief AVM3-based scale separation

  */
  void AVM3Separation();

  /*!
  \brief calculate lift&drag forces and angular momenta

  */
  void LiftDrag() const;

  /*!
  \brief calculate traction values at Dirichlet boundaries

  */
  Teuchos::RCP<Epetra_Vector> CalcStresses();

  /// integrate shape functions at nodes marked by condition
  /*!
    Needed for Mortar coupling at the FSI interface
   */
  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape(std::string condname);

  /// switch fluid field to block matrix
  void UseBlockMatrix(Teuchos::RCP<std::set<int> >     condelements,
                      const LINALG::MultiMapExtractor& domainmaps  ,
                      const LINALG::MultiMapExtractor& rangemaps   ,
                      bool                             splitmatrix);

  /// linear solve with prescribed dirichlet conditions and without history
  /*!
    This is the linear solve as required for steepest descent FSI.
   */
  void LinearRelaxationSolve(Teuchos::RCP<Epetra_Vector> relax);


  /*!
    Return (rotatory) transformation matrix of local co-ordinate systems if
    existent. Else return Teuchos::null.

   */
  Teuchos::RCP<const LINALG::SparseMatrix> GetLocSysTrafo() const
  {
    return (locsysman_ != Teuchos::null) ? locsysman_->Trafo() : Teuchos::null;
  }

  void SetScalarField(
     RCP<const Epetra_Vector> scalarnp,
     Teuchos::RCP<DRT::Discretization> scatradis,
     const int whichscalar = -1);

  /// provide access to turbulence statistics manager
  Teuchos::RCP<FLD::TurbulenceStatisticManager> TurbulenceStatisticManager();
  /// provide access to the box filter for dynamic Smagorinsky model
  Teuchos::RCP<FLD::DynSmagFilter> DynSmagFilter();

  //! @name methods for fsi
  /// Extrapolation of vectors from mid-point to end-point t_{n+1}
  Teuchos::RCP<Epetra_Vector> ExtrapolateEndPoint(
      Teuchos::RCP<Epetra_Vector> vecn, ///< vector at time level t_n
      Teuchos::RCP<Epetra_Vector> vecm ///< vector at time level of equilibrium
      );
  //@}




  void Timeloop() { GenAlphaTimeloop(); }
  Teuchos::RCP<const Epetra_Vector> RHS(){  return Residual(); }
  Teuchos::RCP<const Epetra_Vector> ConvectiveVel(){
    if (GridVel() == Teuchos::null)
      return Velnp(); // no moving mesh present
    else
    {
      // make an intermediate copy of velnp
      Teuchos::RCP<Epetra_Vector> convel = Teuchos::rcp(new Epetra_Vector(*(Velnp())));
      // now subtract the grid velocity
      convel->Update(-1.0,*(GridVel()),1.0);

      return convel;
    }
  }
  Teuchos::RCP<const Epetra_Map> DofRowMap()
  {
    const Epetra_Map* dofrowmap = discret_->DofRowMap();
    return Teuchos::rcp(dofrowmap, false);
  }
  Teuchos::RCP<const LINALG::MapExtractor> GetDBCMapExtractor()
  {
    return DirichMaps();
  }
  void PrepareTimeStep()
  {
    GenAlphaPrepareTimeStep();
  }
  void Evaluate(Teuchos::RCP<const Epetra_Vector> accstepinc)
  {
    if (accstepinc!=Teuchos::null)
    {
      ExternIncrementOfVelnp(accstepinc);
    }

    GenAlphaComputeIntermediateSol();
    GenAlphaAssembleResidualAndMatrix();
  }
  void Update()
  {
    Teuchos::RCP<Epetra_Vector> interfaceforcem = interface_->ExtractFSICondVector(TrueResidual());

    interfaceforcen_ = ExtrapolateEndPoint(interfaceforcen_,interfaceforcem);

    GenAlphaTimeUpdate();
  }
  void StatisticsAndOutput()
  {
    GenAlphaStatisticsAndOutput();
  }
  void Output()
  {
    GenAlphaOutput();
  }
  void NonlinearSolve()
  {
    DoGenAlphaPredictorCorrectorIteration();
  }
  Teuchos::RCP<const Epetra_Map> InnerVelocityRowMap()
  {
    return innervelmap_;
  }
  void SetMeshMap(Teuchos::RCP<const Epetra_Map> mm)
  {
    meshmap_.Setup(*dis_->DofRowMap(),mm,LINALG::SplitMap(*dis_->DofRowMap(),*mm));
  }
  double TimeScaling() const
  {
    const double dt = Dt();
    const double gamma = Gamma();
    return 1./(dt*dt*gamma);
  }
  Teuchos::RCP<Epetra_Vector> ExtractInterfaceForces()
  {
    Teuchos::RCP<Epetra_Vector> interfaceforcem = interface_->ExtractFSICondVector(TrueResidual());

    return ExtrapolateEndPoint(interfaceforcen_,interfaceforcem);
  }
  Teuchos::RCP<Epetra_Vector> ExtractInterfaceVelnp()
  {
    return interface_->ExtractFSICondVector(Velnp());
  }
  Teuchos::RCP<Epetra_Vector> ExtractInterfaceVeln()
  {
    return interface_->ExtractFSICondVector(Veln());
  }
  void ApplyInterfaceVelocities(Teuchos::RCP<Epetra_Vector> ivel)
  {
    // --------------------------------------------------
    // apply new Dirichlet values to velnp according to interface velocity
    interface_->InsertFSICondVector(ivel,velnp_);

    // adjust accnp according to new Dirichlet values of velnp
    //
    //                                  n+1     n
    //                               vel   - vel
    //       n+1      n  gamma-1.0      (0)
    //    acc    = acc * --------- + ------------
    //       (0)           gamma      gamma * dt
    //
    GenAlphaCalcInitialAccelerations();

    // this is very easy, but there are two dangers:
    // - We change ivel here. It must not be used afterwards.
    // - The algorithm must support the change of Dirichlet DOFs
    AddDirichCond(interface_->FSICondMap());
  }
  void ApplyMeshDisplacement(Teuchos::RCP<const Epetra_Vector> fluiddisp)
  {
    meshmap_.InsertCondVector(fluiddisp,dispnp_);

    // new grid velocity
    UpdateGridv();
  }
  void ApplyMeshVelocity(Teuchos::RCP<const Epetra_Vector> gridvel)
  {
    meshmap_.InsertCondVector(gridvel,gridvelaf_);
  }
  void DisplacementToVelocity(Teuchos::RCP<Epetra_Vector> fcx)
  {
    // We convert Delta d(n+1,i+1) to Delta a(n+1,i+1) here.
    //
    // Delta d(n+1,i+1) =   ( Delta u(n+1,i+1)              + u(n) ) * dt
    //
    //                  = ( ( Delta a(n+1,i+1) * gamma * dt + u(n) ) * dt
    //

    const double dt = Dt();
    const double gamma = Gamma();

    // get interface velocity at t(n)
    Teuchos::RCP<Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

    // reduce to Delta u(n+1,i+1)
    fcx->Update(-1.,*veln,1./dt);

    // reduce to Delta a(n+1,i+1)
    fcx->Scale(1./dt/gamma);
  }
  /// Conversion of displacement to velocity at the interface considering predictors or inhomogeneous DBCs
  void DisplacementToVelocity(
      Teuchos::RCP<Epetra_Vector> fcx,    ///< fluid displacement step increment
      Teuchos::RCP<Epetra_Vector> ddgpre, ///< displacement increment due to structural predictor
      Teuchos::RCP<Epetra_Vector> dugpre  ///< velocity increment due to fluid predictor
  );

  void VelocityToDisplacement(Teuchos::RCP<Epetra_Vector> fcx)
  {
    const double dt = Dt();
    const double gamma = Gamma();

    // get interface velocity at t(n)
    Teuchos::RCP<Epetra_Vector> veln = Interface()->ExtractFSICondVector(Veln());

    fcx->Update(1.,*veln,dt*gamma);
    fcx->Scale(dt);
  }

  /// Conversion of velocity to displacement at the interface considering predictors or inhomogeneous DBCs
  void VelocityToDisplacement(
      Teuchos::RCP<Epetra_Vector> fcx,    ///< fluid velocity step increment at interface
      Teuchos::RCP<Epetra_Vector> ddgpre, ///< displacement increment due to structural predictor
      Teuchos::RCP<Epetra_Vector> dugpre  ///< velocity increment due to fluid predictor
  );

  Teuchos::RCP<Epetra_Vector> IntegrateInterfaceShape()
  {
    return interface_->ExtractFSICondVector(IntegrateInterfaceShape("FSICoupling"));
  }
  void UseBlockMatrix(bool splitmatrix)
  {
    Teuchos::RCP<std::set<int> > condelements = Interface()->ConditionedElementMap(*Discretization());
    UseBlockMatrix(condelements,*Interface(),*Interface(),splitmatrix);
  }
  Teuchos::RCP<Epetra_Vector> RelaxationSolve(Teuchos::RCP<Epetra_Vector> ivel)
  {
    const Epetra_Map* dofrowmap = Discretization()->DofRowMap();
    Teuchos::RCP<Epetra_Vector> relax = LINALG::CreateVector(*dofrowmap,true);
    interface_->InsertFSICondVector(ivel,relax);
    LinearRelaxationSolve(relax);
    return ExtractInterfaceForces();
  }
  Teuchos::RCP<DRT::ResultTest> CreateFieldTest()
  {
    return Teuchos::rcp(new FluidResultTest(*this));
  }
  Teuchos::RCP<const Epetra_Vector> ExtractVelocityPart(Teuchos::RCP<const Epetra_Vector> velpres)
  {
     return (VelPresSplitter()).ExtractOtherVector(velpres);
  }
  void SetTimeLomaFields(
      RCP<const Epetra_Vector> scalarnp,
      const double             thermpressnp,
      RCP<const Epetra_Vector> scatraresidual,
      Teuchos::RCP<DRT::Discretization> scatradis,
      const int whichscalar)
  {
      SetScalarField(scalarnp,scatradis, whichscalar);
      return;
  };
  void Integrate() { GenAlphaTimeloop(); }
  virtual Teuchos::RCP<FLD::UTILS::MapExtractor>const& Interface() const { return interface_; }
  virtual Teuchos::RCP<LINALG::Solver>  LinearSolver() { return solver_; }
  virtual void PrepareSolve()        { return; }

  Teuchos::RCP<FLD::UTILS::MapExtractor>               interface_;
  Teuchos::RCP<Epetra_Vector>            interfaceforcen_;  /// interface force at old time level t_n
  LINALG::MapExtractor                   meshmap_;
  Teuchos::RCP<DRT::Discretization>      dis_;
  Teuchos::RefCountPtr<Epetra_Map>       innervelmap_;



private:

  //@}

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementGeneralFluidParameter();

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over a time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementTimeParameter();

  //! @name Set general parameter in class f3Parameter
  /*!

  \brief parameter (fix over all time step) are set in this method.
         Therefore, these parameter are accessible in the fluid element
         and in the fluid boundary element

  */
  void SetElementTurbulenceParameter();

protected:

private:
  //! @name input arguments of the constructor
  //
  //! project nullspace (useful for completely Dirichlet bounded problems)
  bool                        project_;
  //! integrated basis vector (pressure component)
  Teuchos::RCP<Epetra_Vector> w_;
  //! basis vector for nullspace (vector of ones in the pressure component)
  Teuchos::RCP<Epetra_Vector> c_;

  //! decides whether we use an Eulerian or an ALE formulation
  bool                        alefluid_;
  //@}

  //! the dimension of the problem --- determines which dofs
  //are pressure dofs
  int                         numdim_ ;

  //! fluid density
  double                      density_;

  //! @name time integration parameters
  //
  //! max. number of time steps
  int                         endstep_;
  //! max. simulation time
  double                      endtime_;
  //! current time step size
  double                      dt_;

  //! parameter for generalised alpha --- alpha "momentum"
  double                      alphaM_;
  //! parameter for generalised alpha --- alpha "force"
  double                      alphaF_;
  //! parameter for generalised alpha --- computed from the latter
  //in order to achieve 2nd order accuracy
  double                      gamma_;
  //@}

  //! @name parameters for nonlinear iteration
  //
  //! parameter for linearisation
  INPAR::FLUID::LinearisationAction   newton_;

  //! predictor used in nonlinear iteration
  string                      predictor_;

  //! in the case of a constant incremnt predictor we need velnm. This bool
  //! indicates whether it has been computed before or whether we need a
  //! different start predictor
  bool                        velnm_available_;

  //! parameters governing the nonlinear iteration --- iteration counter
  int                         itenum_      ;
  //! parameters governing the nonlinear iteration --- convergence criterion
  double                      ittol_       ;

  // these variables are storing error norms
  //! L2 error norm --- residual associated with pressure dofs, scaled by density
  double                      L2preresnorm_;
  //! L2 error norm --- residual associated with acceleration dofs, scaled by density
  double                      L2velresnorm_;
  //! L2 error norm --- normalized increment in iteration step (velocity dofs)
  double                      L2incvelnorm_;
  //! L2 error norm --- normalized increment in iteration step (pressure dofs)
  double                      L2incprenorm_;

  //! skip the last element call (computation of residual) if
  //!  itemax_ is reached in order to get some speedup for turbulence
  //!  calculations
  bool                        skiplastelecall_;
  //@}

  //! @name IO features

  // write restart every uprestart_ steps --- restartstep_ is
  // the counter used for this purpose

  //! compute and write traction vector
  int                         writestresses_;

  //! @name Turbulence statistics evaluation

  //! parameters for sampling/dumping period
  int                         samstart_;

  //! flag for special flow: currently channel flow, flow in
  //! a lid-driven cavity and flow around a square cylinder
  string                      special_flow_;

  //! flag for potential Neumann-type outflow stabilization
  string                      outflow_stab_;

  //! manager for turbulence statistics
  RCP<FLD::TurbulenceStatisticManager> statisticsmanager_;

  //@}


  //! @name vectors stating the linear system
  //        ---------------------------------
  //
  //! (standard) system matrix
  RCP<LINALG::SparseOperator> sysmat_;

  /// linearization with respect to mesh motion
  RCP<LINALG::BlockSparseMatrixBase> shapederivatives_;

  //! (standard) residual vector (rhs for the incremental form),
  RCP<Epetra_Vector>          residual_;

  //! maps for extracting Dirichlet and free DOF sets
  RCP<LINALG::MapExtractor>   dbcmaps_;

  //! a vector of zeros to be used to enforce zero dirichlet boundary conditions
  RCP<Epetra_Vector>          zeros_;

  //! the vector containing body and surface forces
  RCP<Epetra_Vector>          neumann_loads_;

  //! the vector containing potential Neumann-type outflow stabilization
  RCP<Epetra_Vector>          outflow_stabil_;

  //! Nonlinear iteration increment vector contains accleeration
  //! and pressure increment
  RCP<Epetra_Vector>          increment_;
  //@}


  //! @name Force computation on Dirichlet boundaries
  // -----------------------------------------
  //
  //! the force vector is a copy of the residual without dirichlet
  //! conditions applied --- probably, for all inner nodes,
  //! "equilibrium" would be a better name - but for the Dirichlet
  //! nodes, this vector contains the consistent (normed) force on
  //! boundary nodes
  RCP<Epetra_Vector>          force_;
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
  RCP<Epetra_Vector>          accnp_;
  //! acceleration at time n
  //! in  accn_, every component associated to a pressure dof is
  //! completely meaningless --- the Navier Stokes equations are a differential
  //! algebraic equation, there's no evolution equation for the pressure!
  RCP<Epetra_Vector>          accn_;
  //! acceleration at intermediate time step n+alpha_M
  //! in accam_, every component associated to a pressure dof is
  //! completely meaningless --- the Navier Stokes equations are a differential
  //! algebraic equation, there's no evolution equation for the pressure!
  RCP<Epetra_Vector>          accam_;

  //! velocities and pressures at time n+1
  RCP<Epetra_Vector>          velnp_;
  //! velocities and pressures at time n
  RCP<Epetra_Vector>          veln_;
  //! velocities at intermediate time step n+alpha_F
  //! in velaf_, every component associated to a pressure dof is completely
  //! meaningless since the pressure is treated implicitly.
  RCP<Epetra_Vector>          velaf_;

  //! velocities and pressure at time n-1
  //! required only for the constant increment predictor
  RCP<Epetra_Vector>          velnm_;

  //! scalar at time n+alpha_F/n+1
  //! required only for Neumann boundary conditions
  RCP<Epetra_Vector>          scaaf_;

  //! ALE related stuff: grid displacements, time n+1
  RCP<Epetra_Vector>          dispnp_   ;
  //! ALE related stuff: most recent grid displacements, time n
  RCP<Epetra_Vector>          dispn_    ;
  //! ALE related stuff: most recent grid displacements, time n-1
  RCP<Epetra_Vector>          dispnm_    ;
  //! ALE related stuff: grid velocities, most recent time step n
  RCP<Epetra_Vector>          gridveln_ ;
  //! ALE related stuff: intermediate grid velocities at time n+af
  RCP<Epetra_Vector>          gridvelaf_;
  //@}

  //! fine-scale subgrid-viscosity flag
  string                      fssgv_;

  //! only necessary for AVM3: scale-separation matrix
  RCP<LINALG::SparseMatrix>   Sep_;

  //! only necessary for AVM3: fine-scale solution vector
  RCP<Epetra_Vector>          fsvelaf_;

  //@}

  //! @name Vectors used for convergence check
  //        ----------------------------------

  //! extracts velocities and pressure dofs from dofmap
  LINALG::MapExtractor        velpressplitter_;
  //@}


  //! @name connectivity generated from periodic boundary conditions
  //        --------------------------------------------------------
  //
  //! the connection between master and slave nodes on this proc
  RCP<map<int,vector<int> > > pbcmapmastertoslave_;
  //@}


  //! @name Dirichlet BCs with local co-ordinate system
  //        -------------------------------------------
  //
  //! Dirichlet BCs with local co-ordinate system
  Teuchos::RCP<DRT::UTILS::LocsysManager> locsysman_;
  //@}


  //! @name Turbulent inflow conditions
  //        ---------------------------
  //
  //! a manager doing the transfer of boundary data for
  //! turbulent inflow profiles from a separate (periodic) domain
  Teuchos::RCP<TransferTurbulentInflowCondition> turbulent_inflow_condition_;
  //@}


  //! @name Mixed-hybrid Dirichlet condition
  //        --------------------------------
  //
  //! object for a redistributed evaluation of the mixed-hybrid Dirichlet condition
  RCP<FluidMHDEvaluate> MHD_evaluator_;
  //@}


  //! @name Dynamic Smagorinsky model: methods and variables
  //        -------------------------

  //! one instance of the filter object
  RCP<FLD::DynSmagFilter>                      DynSmag_;

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

  //! count cpu time for element call on proc
  double                    dtele_;
  //! count cpu time for solver call on proc
  double                    dtsolve_;
  //! count cpu time for filtering (in case of a dynamic Smagorinsky model) on proc
  double                    dtfilter_;
  //@}
};

}

#endif  // not FLUID_GENALPHA_INTEGRATION_H

