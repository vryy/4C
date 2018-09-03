/*----------------------------------------------------------------------*/
/*!
\file xfluid_levelset_coupling_algorithm.cpp

\brief Basis of xfluid-levelset coupling.

\level 3

<pre>
\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>
 */
/*----------------------------------------------------------------------*/

#include "../drt_levelset/levelset_algorithm.H"
#include "../drt_levelset/levelset_timint_ost.H"
#include "../drt_scatra/scatra_timint_ost.H"

#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_fluid/fluid_utils_mapextractor.H"  //needed?
#include "../drt_fluid_xfluid/xfluid.H"

#include "../drt_lib/drt_discret_xfem.H"

#include "../drt_mat/newtonianfluid.H"
#include "../drt_mat/matlist.H"

#include "xfluid_levelset_coupling_algorithm.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
XFLUIDLEVELSET::Algorithm::Algorithm(const Epetra_Comm& comm, const Teuchos::ParameterList& prbdyn,
    const Teuchos::ParameterList& solverparams,
    const std::map<std::string, int>& dofset_coupling_map)
    : ScaTraFluidCouplingAlgorithm(comm, prbdyn, false, "scatra", solverparams),
      dt_(0.0),
      maxtime_(0.0),
      stepmax_(0),
      itmax_(0),
      ittol_(1.0),
      upres_(-1),
      write_center_of_mass_(false),
      smoothedgradphitype_(INPAR::TWOPHASE::smooth_grad_phi_l2_projection),
      scalesmoothedgradients_(false),
      velnpi_(Teuchos::null),
      phinpi_(Teuchos::null),
      prbdyn_(prbdyn),
      dofset_coupling_map_(dofset_coupling_map)
{
  // Needs to stay emtpy
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
XFLUIDLEVELSET::Algorithm::~Algorithm() { return; }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Init(
    const Teuchos::ParameterList& prbdyn,        ///< parameter list for global problem
    const Teuchos::ParameterList& scatradyn,     ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
    const std::string& disname,                  ///< name of scalar transport discretization
    const bool isale                             ///< ALE flag
)
{
  // call Setup() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Init(prbdyn, scatradyn, solverparams, disname, isale);


  // TODO: Combine TWOPHASE and XFLUIDLEVELSET. Create a Parent class, TWOFLUIDCOUPLING or use the
  // existing ScaTraFluidCouplingAlgorithm.
  //        Derived classes will be FLUIDLEVELSET and XFLUIDLEVELSET.
  //        This could also incorporate ELCH and other FLUID-ScaTra coupling schemes.

  // ScaTra Field is not given an initial velocity field. Thus it has to be instantiated before the
  // ScaTra field is being solved for.

  // time-step length, maximum time and maximum number of steps
  dt_ = prbdyn_.get<double>("TIMESTEP");
  maxtime_ = prbdyn_.get<double>("MAXTIME");
  stepmax_ = prbdyn_.get<int>("NUMSTEP");

  // Output specific criterions
  write_center_of_mass_ = DRT::INPUT::IntegralValue<bool>(prbdyn_, "WRITE_CENTER_OF_MASS");

  // (preliminary) maximum number of iterations and tolerance for outer iteration
  ittol_ = prbdyn_.get<double>("CONVTOL");
  itmax_ = prbdyn_.get<int>("ITEMAX");

  upres_ = prbdyn_.get<int>("RESULTSEVRY");

  // TODO: Put this into condition manager
  smoothedgradphitype_ = DRT::INPUT::IntegralValue<INPAR::TWOPHASE::SmoothGradPhi>(
      prbdyn.sublist("SURFACE TENSION"), "SMOOTHGRADPHI");
  scalesmoothedgradients_ = DRT::INPUT::IntegralValue<bool>(
      prbdyn.sublist("SURFACE TENSION"), "SCALE_SMOOTHED_GRADIENTS");

  // Instantiate vectors contatining outer loop increment data
  fsvelincnorm_.reserve(itmax_);
  fspressincnorm_.reserve(itmax_);
  fsphiincnorm_.reserve(itmax_);

  return;
}

void XFLUIDLEVELSET::Algorithm::DoAlgorithmSpecificInit()
{
  // TODO: set as member, first initialized in the dyn!
  const int nds_fluid_proxy_in_scatra = dofset_coupling_map_["vel_fluid_proxy_in_scatra"];
  const int nds_xfluid = dofset_coupling_map_["vel_in_xfluid"];

  // TODO: do we need this after similar to the fluid-proxy in the fluid after redistributing scatra
  // due to particles?
  const int nds_phi_scatra = dofset_coupling_map_["phi_in_scatra"];
  const int nds_scatra_proxy_in_fluid = dofset_coupling_map_["phi_scatra_proxy_in_fluid"];



  // we expect a level-set algorithm here, which potentially carries a particle algorithm
  Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField(), true);

  // check whether fluid and scatra discret still have the same maps
  // they may change due a modified ghosting required, i.e., for particle level-set methods

  Teuchos::RCP<DRT::Discretization> fluiddis = FluidField()->Discretization();
  Teuchos::RCP<DRT::Discretization> scatradis = ScaTraField()->Discretization();

  const Epetra_Map* scatraelecolmap = scatradis->ElementColMap();
  const Epetra_Map* fluidelecolmap = fluiddis->ElementColMap();

  // TODO: do we need this?
  fluiddis->ReplaceDofSet(
      nds_scatra_proxy_in_fluid, scatradis->GetDofSetProxy(nds_phi_scatra), false);

  if (not scatraelecolmap->PointSameAs(*fluidelecolmap))
  {
    if (Comm().MyPID() == 0)
      std::cout << "----- Adapt fluid ghosting to scatra ghosting ------" << std::endl;

    // adapt fluid ghosting to scatra ghosting
    fluiddis->ExtendedGhosting(*scatraelecolmap, true, true, true, false);

    Teuchos::RCP<DRT::DiscretizationXFEM> xfluiddis =
        Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(fluiddis, true);

    // Need to call InitialFillComplete again to get the correct Map after the Extended Ghosting of
    // the fluid dofset. calls AssignDegreesOfFreedom for all non-proxy dofsets and the nds
    // initialdofsets
    std::vector<int> nds;
    nds.push_back(nds_xfluid);
    xfluiddis->InitialFillComplete(nds);

    // replace the fluid-dofset proxy in the scatra discretization
    scatradis->ReplaceDofSet(
        nds_fluid_proxy_in_scatra, xfluiddis->GetInitialDofSetProxy(nds_xfluid), false);

    // does not have an effect (as wanted!), as the additional fluid proxy,
    // which has been set in scatra, has its own AssignDegreesOfFreedom that does nothing!
    scatradis->FillComplete(true, false, false);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Setup()
{
  // call Setup() in base class
  ADAPTER::ScaTraFluidCouplingAlgorithm::Setup();

  const int restart = DRT::Problem::Instance()->Restart();

  if (restart) Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField(), true)->CreateInitialState();

  // Fluid-Scatra Iteration vectors are initialized
  velnpi_ = Teuchos::rcp(new Epetra_Vector(FluidField()->StdVelnp()->Map()), true);
  velnpi_->Update(1.0, *FluidField()->StdVelnp(), 0.0);
  phinpi_ = Teuchos::rcp(new Epetra_Vector(ScaTraField()->Phinp()->Map()), true);
  phinpi_->Update(1.0, *ScaTraField()->Phinp(), 0.0);

  // do not change the following order of calls
  // this cannot be done in Init() for some reason


  if (!restart)  // otherwise SetScaTraValuesInFluid would overwrite the restarted vectors in the
                 // cond-manager
  {
    // 1.) set reconstructed normals and curvature values of the initial scatra field
    // to allow to obtain a right transport velocity field in the next call when setting fluid
    // values in the scatra field.
    SetScaTraValuesInFluid();

    // 2.) set the interface transport velocity in the scatra field

    SetFluidValuesInScaTra(true);
  }

  return;
}


/*---------------------------------------------------------------------------------------*
| public: algorithm for a instationary XTPF problem                         winter 10/14 |
*----------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::TimeLoop()
{
  OutputInitialField();

  // time loop
  while (NotFinished())
  {
    IncrementTimeAndStep();

    // prepare time step
    PrepareTimeStep();

    // do outer iteration loop for particular type of algorithm
    OuterLoop();

    // update for next time step
    TimeUpdate();

    // write output to files
    Output();

  }  // time loop

  return;
}

/*------------------------------------------------------------------------------------------------*
 | public: algorithm for a stationary XTPF problem                                   winter 10/14 |
 *------------------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SolveStationaryProblem()
{
  if (Comm().MyPID() == 0)
  {
    printf(
        "------Stationary-Xfluid-LevelSet-XFEM------  time step "
        "----------------------------------------\n");
  }

  // check time integration schemes of single fields
  // remark: this was already done in ScaTraFluidCouplingAlgorithm() before
  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary)
    dserror("Scatra time integration scheme is not stationary");

  // write Scatra output (fluid output has been already called in FluidField()->Integrate();
  ScaTraField()->Output();

  // Give Scatra Values to fluid.
  // Needed for curvature etc..
  SetScaTraValuesInFluid();

  // run the simulation, calls the xfluid-"integrate()" routine
  FluidField()->Integrate();

  //  // solve level set equation
  if (Comm().MyPID() == 0)
    std::cout << "/!\\ warning === Level-set field not solved for Fluid_XFEM_LevelSet problems"
              << std::endl;


  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::OuterLoop()
{
  int itnum = 0;
  bool stopnonliniter = false;

  PrintToScreen(
      "\n****************************************\n          OUTER ITERATION "
      "LOOP\n****************************************\n");

  if (Comm().MyPID() == 0)
  {
    printf("TIME: %11.4E/%11.4E  DT = %11.4E  %s  STEP = %4d/%4d\n", Time(), maxtime_, dt_,
        ScaTraField()->MethodTitle().c_str(), Step(), stepmax_);
  }

  // initially solve scalar transport equation
  // (values for intermediate time steps were calculated at the end of PerpareTimeStep)

  // We always do s - f - s or more cyles s - [ f - s] -... -[f - s]
  // Set relevant Fluid values in ScaTra field. (first transport with predicted fluid velocity = vel
  // from last time-step)
  // TODO: issue (the fluid solution and its interface position is not the same as the last computed
  // for the scatra!!!) if we would use the transport velocity directly from the interface, the
  // positions do not match!

  bool scatra_init = false;

  // initial time-derivative needed for particles restart
  const int restart = DRT::Problem::Instance()->Restart();
  if (Step() == restart + 1) scatra_init = true;

  SetFluidValuesInScaTra(scatra_init);

  DoScaTraField();

  // Prepare variables for convergence check.
  PrepareOuterIteration();

  while (stopnonliniter == false)
  {
    itnum++;

    // Set relevant ScaTra values in Fluid field.
    SetScaTraValuesInFluid();

    // solve fluid flow equations
    DoFluidField();

    // Set relevant Fluid values in ScaTra field.
    SetFluidValuesInScaTra(false);

    // solve scalar transport equation
    DoScaTraField();

    // check convergence and stop iteration loop if convergence is achieved
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::DoFluidField()
{
  PrintToScreen(
      "\n****************************************\n              FLUID "
      "SOLVER\n****************************************\n");

  // Solve the Fluid field.
  FluidField()->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::DoScaTraField()
{
  PrintToScreen(
      "\n****************************************\n        SCALAR TRANSPORT "
      "SOLVER\n****************************************\n");

  // Solve the ScaTra field.
  ScaTraField()->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::PrintToScreen(const std::string str)
{
  if (Comm().MyPID() == 0) std::cout << str;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::TimeUpdate()
{
  // update scalar
  ScaTraField()->Update();

  // update fluid
  FluidField()->Update();

  return;
}

/*---------------------------------------------------------------------------------------*
| Prepares values and variables needed in the outer iteration                            |
*----------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::PrepareOuterIteration()
{
  //  //Update phi for outer loop convergence check
  phinpi_->Update(1.0, *ScaTraField()->Phinp(), 0.0);
  velnpi_->Update(1.0, *FluidField()->StdVelnp(), 0.0);

  // Clear the vectors containing the data for the partitioned increments
  fsvelincnorm_.clear();
  fspressincnorm_.clear();
  fsphiincnorm_.clear();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::PrepareTimeStep()
{
  // prepare scalar transport time step
  // (+ computation of initial scalar time derivative in first time step)
  ScaTraField()->PrepareTimeStep();

  // prepare fluid time step, among other things, predict velocity field
  FluidField()->PrepareTimeStep();

  // synchronicity check between algorithm and fields
  SynchronicityTimeCheck();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SynchronicityTimeCheck()
{
  if (FluidField()->Time() != Time())
    dserror("Time in Fluid time integration differs from time in two phase flow algorithm");
  if (ScaTraField()->Time() != Time())
    dserror("Time in ScaTra time integration differs from time in two phase flow algorithm");
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SetScaTraValuesInFluid()
{
  const Teuchos::RCP<const Epetra_Vector>& scatra_phinp = ScaTraField()->Phinp();

  // compute smoothed geometric quantities based on scatra discretization
  Teuchos::RCP<Epetra_MultiVector> scatra_smoothedgradphi = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> scatra_nodalcurvature = Teuchos::null;

  // currently vectors based on the cutter dis stored in the conditionmanager
  Teuchos::RCP<Epetra_Vector> fluid_phinp = Teuchos::null;
  Teuchos::RCP<Epetra_MultiVector> fluid_smoothedgradphi = Teuchos::null;
  Teuchos::RCP<Epetra_Vector> fluid_nodalcurvature = Teuchos::null;

  Teuchos::RCP<FLD::XFluid> xfluid = Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField(), true);

  // TODO: Get smoothedgradphitype from condition manager
  // get access to the fluid node based vectors
  xfluid->WriteAccess_GeometricQuantities(fluid_phinp, fluid_smoothedgradphi, fluid_nodalcurvature);

  // are smoothed geometric quantities required?
  bool require_smoothedgradphi = (fluid_smoothedgradphi == Teuchos::null) ? false : true;
  bool require_nodalcurvature = (fluid_nodalcurvature == Teuchos::null) ? false : true;

  // set level set in fluid field
  ComputeGeometricQuantities(require_smoothedgradphi, require_nodalcurvature, scatra_phinp,
      smoothedgradphitype_, scatra_smoothedgradphi, scatra_nodalcurvature);

#if (0)
  {
    std::cout << "scatra_phinp " << *scatra_phinp << std::endl;

    if (require_smoothedgradphi)
      std::cout << "scatra_smoothedgradphi " << *scatra_smoothedgradphi << std::endl;

    if (require_nodalcurvature)
      std::cout << "scatra_nodalcurvature " << *scatra_nodalcurvature << std::endl;
  }
#endif

  // map geometric quantities from scatra dof-based to fluid node based vectors
  CopyGeometricQuantities(require_smoothedgradphi, require_nodalcurvature, scatra_phinp,
      scatra_smoothedgradphi, scatra_nodalcurvature, fluid_phinp, fluid_smoothedgradphi,
      fluid_nodalcurvature);

  // export row to col vectors
  xfluid->ExportGeometricQuantities();

  return;
}

/*----------------------------------------------------------------------*
 | Set relevant values from ScaTra field in the Fluid field.            |
 *----------------------------------------------------------------------*/
Teuchos::RCP<Epetra_MultiVector> XFLUIDLEVELSET::Algorithm::GetSmoothedLevelSetGradient(
    const Teuchos::RCP<const Epetra_Vector>& phinp, INPAR::TWOPHASE::SmoothGradPhi smoothedgradphi)
{
  bool normalize_gradients = scalesmoothedgradients_;

  switch (smoothedgradphi)
  {
    case INPAR::TWOPHASE::smooth_grad_phi_l2_projection:
    {
      return ScaTraField()->ReconstructGradientAtNodesL2Projection(phinp, normalize_gradients);
      break;
    }
    case INPAR::TWOPHASE::smooth_grad_phi_superconvergent_patch_recovery_3D:
    {
      return ScaTraField()->ReconstructGradientAtNodesPatchRecon(phinp, 3, normalize_gradients);
      break;
    }
    case INPAR::TWOPHASE::smooth_grad_phi_superconvergent_patch_recovery_2Dz:
    {
      return ScaTraField()->ReconstructGradientAtNodesPatchRecon(phinp, 2, normalize_gradients);
      break;
    }
    case INPAR::TWOPHASE::smooth_grad_phi_meanvalue:
    {
      return ScaTraField()->ReconstructGradientAtNodesMeanAverage(phinp, normalize_gradients);
      break;
    }
    default:
      dserror("The chosen smoothing is not as of yet supported!");
      break;
  }

  return Teuchos::null;
}

/*----------------------------------------------------------------------*
 | Set relevant values from Fluid field in the ScaTra field.            |
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::ComputeGeometricQuantities(const bool require_smoothedgradphi,
    const bool require_nodalcurvature, const Teuchos::RCP<const Epetra_Vector>& phinp,
    const INPAR::TWOPHASE::SmoothGradPhi smoothedgradphitype,
    Teuchos::RCP<Epetra_MultiVector>& smoothedgradphi, Teuchos::RCP<Epetra_Vector>& nodalcurvature)
{
  smoothedgradphi = Teuchos::null;
  nodalcurvature = Teuchos::null;

  // compute the smoothed geometric quantities
  if (require_smoothedgradphi)
    smoothedgradphi = GetSmoothedLevelSetGradient(phinp, smoothedgradphitype);

  if (require_nodalcurvature)
    nodalcurvature =
        Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField(), true)
            ->GetNodalCurvature(phinp, GetSmoothedLevelSetGradient(phinp, smoothedgradphitype));

  return;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::CopyGeometricQuantities(const bool require_smoothedgradphi,
    const bool require_nodalcurvature, const Teuchos::RCP<const Epetra_Vector>& scatra_phinp,
    const Teuchos::RCP<const Epetra_MultiVector>& scatra_smoothedgradphi,
    const Teuchos::RCP<const Epetra_Vector>& scatra_nodalcurvature,
    const Teuchos::RCP<Epetra_Vector>& fluid_phinp,
    const Teuchos::RCP<Epetra_MultiVector>& fluid_smoothedgradphi,
    const Teuchos::RCP<Epetra_Vector>& fluid_nodalcurvature)
{
  // this is simply a transform method between quantities w.r.t scatra dis to quantities in the
  // condition manager

  ///-------------------------------------------
  // copy the level set values

  if (fluid_phinp == Teuchos::null) dserror("fluid_phinp null pointer");
  if (scatra_phinp == Teuchos::null) dserror("scatra_phinp null pointer");

  // safety checks for matching dof-based maps
  if (!fluid_phinp->Map().PointSameAs(scatra_phinp->Map()))
    dserror("Unequal maps: fluid_phinp  -- not PointSameAs --  scatra_phinp");

  fluid_phinp->Update(1.0, *scatra_phinp, 0.0);


  ///-------------------------------------------
  // copy the nodal gradients
  if (require_smoothedgradphi)
  {
    if (fluid_smoothedgradphi == Teuchos::null) dserror("fluid_smoothedgradphi null pointer");
    if (scatra_smoothedgradphi == Teuchos::null) dserror("scatra_smoothedgradphi null pointer");

    if (!fluid_smoothedgradphi->Map().PointSameAs(scatra_smoothedgradphi->Map()))
    {
      dserror("Unequal maps: fluid_smoothedgradphi  -- not PointSameAs --  scatra_smoothedgradphi");
    }

    fluid_smoothedgradphi->Update(1.0, *scatra_smoothedgradphi, 0.0);
  }

  ///-------------------------------------------
  // copy the nodal curvatures

  if (require_nodalcurvature)
  {
    if (fluid_nodalcurvature == Teuchos::null) dserror("fluid_nodalcurvature null pointer");
    if (scatra_nodalcurvature == Teuchos::null) dserror("scatra_nodalcurvature null pointer");

    if (!fluid_nodalcurvature->Map().PointSameAs(scatra_nodalcurvature->Map()))
      dserror("Unequal maps: fluid_nodalcurvature  -- not PointSameAs --  scatra_nodalcurvature");

    fluid_nodalcurvature->Update(1.0, *scatra_nodalcurvature, 0.0);
  }

  return;
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::SetFluidValuesInScaTra(bool init)
{
  if (Comm().MyPID() == 0) std::cout << "\n------  SetFluidValuesInScaTra ---------\n";


  Teuchos::RCP<SCATRA::LevelSetAlgorithm> levelsetalgo =
      Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField(), true);

  Teuchos::RCP<FLD::XFluid> xfluid = Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField(), true);

  levelsetalgo->SetVelocityField(
      xfluid->GetTransportVelocity(), Teuchos::null, Teuchos::null, Teuchos::null, false,
      init  // TODO: what does this init??? how to set initial accelerations in the scatra?
  );

  return;
}



/*------------------------------------------------------------------------------------------------*
 | Fluid - ScaTra partitioned convergence check                                                   |
 |                                                                                                |
 |    The increment between outer iterations is checked, and if lower than given tolerance        |
 |    the loop is exited. However, at least 2 solutions of each field is required to perform a    |
 |    convergence check.                                                                          |
 |                                                                                   winter 02/15 |
 *------------------------------------------------------------------------------------------------*/
bool XFLUIDLEVELSET::Algorithm::ConvergenceCheck(int itnum)
{
  // define flags for Fluid and ScaTra convergence check
  bool fluidstopnonliniter = false;
  bool scatrastopnonliniter = false;

  bool notconverged = false;

  if (itmax_ <= 0) dserror("Set iterations to something reasonable!!!");

  double fsvelincnorm = 1.0;
  double fspressincnorm = 1.0;
  double fsphiincnorm = 1.0;
  // Get increment for outer loop of Fluid and ScaTra
  GetOuterLoopIncFluid(fsvelincnorm, fspressincnorm, itnum);
  GetOuterLoopIncScaTra(fsphiincnorm, itnum);

  fsvelincnorm_[itnum - 1] = fsvelincnorm;
  fspressincnorm_[itnum - 1] = fspressincnorm;
  fsphiincnorm_[itnum - 1] = fsphiincnorm;

  if (Comm().MyPID() == 0)
  {
    printf(
        "\n|+------ TWO PHASE FLOW CONVERGENCE CHECK:  time step %2d, outer iteration %2d ------+|",
        Step(), itnum);
    printf(
        "\n|- iter/itermax -|----tol-[Norm]---|-- fluid-inc --|-- press inc --|-- levset inc --|");
  }

  for (int k_itnum = 0; k_itnum < itnum; k_itnum++)
  {
    if (k_itnum == 0)
    {
      if (Comm().MyPID() == 0)
      {
        printf("\n|     %2d/%2d      | %10.3E [L2] |       -       |       -       |   %10.3E   |",
            (k_itnum + 1), itmax_, ittol_, fsphiincnorm_[k_itnum]);
      }  // end if processor 0 for output
    }
    else
    {
      if (Comm().MyPID() == 0)
      {
        printf("\n|     %2d/%2d      | %10.3E [L2] |  %10.3E   |  %10.3E   |   %10.3E   |",
            (k_itnum + 1), itmax_, ittol_, fsvelincnorm_[k_itnum], fspressincnorm_[k_itnum],
            fsphiincnorm_[k_itnum]);
      }  // end if processor 0 for output
    }
  }
  if (Comm().MyPID() == 0)
    printf(
        "\n|+---------------------------------------------------------------------------------+|"
        "\n");


  if ((fsvelincnorm <= ittol_) and (fspressincnorm <= ittol_) and itnum > 1)
    fluidstopnonliniter = true;

  if ((fsphiincnorm <= ittol_)) scatrastopnonliniter = true;


  // If tolerance or number of maximum iterations are reached
  if ((fluidstopnonliniter and scatrastopnonliniter) or (itnum >= itmax_))
  {
    notconverged = true;
  }

  if (Comm().MyPID() == 0)
  {
    if ((itnum == stepmax_) and (notconverged == true))
    {
      printf("|+---------------- not converged ----------------------+|");
      printf("\n|+-----------------------------------------------------+|\n");
    }
  }

  return notconverged;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Output()
{
  if (Step() % upres_ ==
      0)  // Only perform output for given RESULTSEVRY in Control Algo section of input.
  {
    FluidField()->Output();
    ScaTraField()->Output();
  }

  if (write_center_of_mass_)
  {
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->MassCenterUsingSmoothing();
  }

  return;
}

/*------------------------------------------------------------------------------------------------*
 | protected: output of initial field                                                winter 07/14 |
 *------------------------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::OutputInitialField()
{
  if (Step() == 0)
  {
    // output fluid initial state
    if (FluidField()->TimIntScheme() != INPAR::FLUID::timeint_stationary) FluidField()->Output();

    // output Levelset function initial state
    if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_stationary) ScaTraField()->Output();
  }

  if (write_center_of_mass_)
  {
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->MassCenterUsingSmoothing();
  }

  return;
}

/*----------------------------------------------------------------------*
 | perform result test                                     winter 06/14 |
 *----------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::TestResults()
{
  // perform result tests if required
  DRT::Problem::Instance()->AddFieldTest(FluidField()->CreateFieldTest());
  // DRT::Problem::Instance()->TestAll(Comm());

  if (ScaTraField()->MethodName() != INPAR::SCATRA::timeint_gen_alpha)
  {
    Teuchos::rcp_dynamic_cast<SCATRA::LevelSetAlgorithm>(ScaTraField())->TestResults();
  }
  else
  {
    dserror("Unknown time integration for Level Set field in Two Phase Flow problems.");
  }

  return;
}

/* -------------------------------------------------------------------------------*
 | Restart a X-two phase problem                                              winter|
 * -------------------------------------------------------------------------------*/
void XFLUIDLEVELSET::Algorithm::Restart(int step)
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "##################################################" << std::endl;
    std::cout << "#                                                #" << std::endl;
    std::cout << "#     Restart of T(wo) P(hase) F(low) problem    #" << std::endl;
    std::cout << "#                                                #" << std::endl;
    std::cout << "##################################################" << std::endl;

    std::cout << "##########################################################################"
              << std::endl;
    std::cout << "#                                                                        #"
              << std::endl;
    std::cout << "#     WARNING: ONLY RESTART FROM XFLUID AND SCATRA ALLOWED FOR NOW!!!    #"
              << std::endl;
    std::cout << "#                                                                        #"
              << std::endl;
    std::cout << "##########################################################################"
              << std::endl;
  }

  // condition-manager restart has been called alredy in the xfluid's init!

  FluidField()->ReadRestart(step);
  ScaTraField()->ReadRestart(step);

  SetTimeStep(FluidField()->Time(), step);

  // TODO: issue: this however will overwrite the scatra solution (the xfluid cutter-state however
  // is behind the scatra state!) Needed for particle restart, try to comment!
  //  SetFluidValuesInScaTra(true);

  return;
}


void XFLUIDLEVELSET::Algorithm::GetOuterLoopIncFluid(
    double& fsvelincnorm, double& fspressincnorm, int itnum)
{
  Teuchos::RCP<const Epetra_Vector> velnpip =
      FluidField()->StdVelnp();  // Contains Fluid and Pressure

  // Extract velocity and pressure components.
  Teuchos::RCP<const LINALG::MapExtractor> velpresspliter =
      Teuchos::rcp_dynamic_cast<FLD::XFluid>(FluidField())->VelPresSplitterStd();

  Teuchos::RCP<const Epetra_Vector> onlyvel = velpresspliter->ExtractOtherVector(velnpip);
  Teuchos::RCP<const Epetra_Vector> onlypress = velpresspliter->ExtractCondVector(velnpip);

  Teuchos::RCP<Epetra_Vector> onlyveli = Teuchos::rcp(new Epetra_Vector(onlyvel->Map()), true);
  Teuchos::RCP<Epetra_Vector> onlypressi = Teuchos::rcp(new Epetra_Vector(onlypress->Map()), true);

  if (itnum > 1)
  {
    onlyveli = velpresspliter->ExtractOtherVector(velnpi_);
    onlypressi = velpresspliter->ExtractCondVector(velnpi_);
  }

  double velnormL2 = 1.0;
  double pressnormL2 = 1.0;

  onlyvel->Norm2(&velnormL2);
  onlypress->Norm2(&pressnormL2);

  if (velnormL2 < 1e-5) velnormL2 = 1.0;
  if (pressnormL2 < 1e-5) pressnormL2 = 1.0;

  double fsvelnormL2 = 1.0;
  double fspressnormL2 = 1.0;

  // compute increment and L2-norm of increment
  //-----------------------------------------------------
  Teuchos::RCP<Epetra_Vector> incvel = Teuchos::rcp(new Epetra_Vector(onlyvel->Map()), true);
  incvel->Update(1.0, *onlyvel, -1.0, *onlyveli, 0.0);
  incvel->Norm2(&fsvelnormL2);

  Teuchos::RCP<Epetra_Vector> incpress = Teuchos::rcp(new Epetra_Vector(onlypress->Map()), true);
  incpress->Update(1.0, *onlypress, -1.0, *onlypressi, 0.0);
  incpress->Norm2(&fspressnormL2);
  //-----------------------------------------------------

  fsvelincnorm = fsvelnormL2 / velnormL2;
  fspressincnorm = fspressnormL2 / pressnormL2;

#if DEBUG
  //-------------------------
  std::cout << "fsvelnormL2: " << fsvelnormL2 << std::endl;
  std::cout << "velnormL2: " << velnormL2 << std::endl << std::endl;

  std::cout << "fspressnormL2: " << fspressnormL2 << std::endl;
  std::cout << "pressnormL2: " << pressnormL2 << std::endl << std::endl;
//-------------------------
#endif

  velnpi_->Update(1.0, *velnpip, 0.0);
}

void XFLUIDLEVELSET::Algorithm::GetOuterLoopIncScaTra(double& fsphiincnorm, int itnum)
{
  Teuchos::RCP<const Epetra_Vector> phinpip = ScaTraField()->Phinp();

  double phinormL2 = 1.0;

  phinpip->Norm2(&phinormL2);
  if (phinormL2 < 1e-5) phinormL2 = 1.0;
  double fsphinormL2 = 1.0;

  // compute increment and L2-norm of increment
  //-----------------------------------------------------
  Teuchos::RCP<Epetra_Vector> incphi = Teuchos::rcp(new Epetra_Vector(phinpip->Map(), true));
  incphi->Update(1.0, *phinpip, -1.0, *phinpi_, 0.0);
  incphi->Norm2(&fsphinormL2);
  //-----------------------------------------------------

  fsphiincnorm = fsphinormL2 / phinormL2;

#if DEBUG
  //-------------------------
  std::cout << "fsphinormL2: " << fsphinormL2 << std::endl;
  std::cout << "phinormL2: " << phinormL2 << std::endl << std::endl;
//-------------------------
#endif


  phinpi_->Update(1.0, *phinpip, 0.0);
}
