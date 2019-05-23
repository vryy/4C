/*----------------------------------------------------------------------*/
/*!
\file elch_moving_boundary_algorithm.cpp

\brief Basis of all ELCH algorithms with moving boundaries

\level 2

\maintainer Christoph Schmidt
            schmidt@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
*/
/*----------------------------------------------------------------------*/


#include "elch_moving_boundary_algorithm.H"
#include "../drt_inpar/inpar_elch.H"
#include "../drt_fluid/fluid_utils_mapextractor.H"
#include "../linalg/linalg_utils.H"
#include "../drt_io/io.H"

#include "../drt_scatra/scatra_timint_elch.H"
//#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::MovingBoundaryAlgorithm::MovingBoundaryAlgorithm(const Epetra_Comm& comm,
    const Teuchos::ParameterList& elchcontrol, const Teuchos::ParameterList& scatradyn,
    const Teuchos::ParameterList& solverparams)
    : ScaTraFluidAleCouplingAlgorithm(comm, scatradyn, "FSICoupling", solverparams),
      pseudotransient_(false),
      molarvolume_(elchcontrol.get<double>("MOLARVOLUME")),
      idispn_(Teuchos::null),
      idispnp_(Teuchos::null),
      iveln_(Teuchos::null),
      itmax_(elchcontrol.get<int>("MOVBOUNDARYITEMAX")),
      ittol_(elchcontrol.get<double>("MOVBOUNDARYCONVTOL")),
      theta_(elchcontrol.get<double>("MOVBOUNDARYTHETA")),
      elch_params_(elchcontrol)
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ELCH::MovingBoundaryAlgorithm::~MovingBoundaryAlgorithm() { return; }



/*----------------------------------------------------------------------*
| Setup                                                     rauch 08/16 |
*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Init(
    const Teuchos::ParameterList& prbdyn,        ///< parameter list for global problem
    const Teuchos::ParameterList& scatradyn,     ///< parameter list for scalar transport subproblem
    const Teuchos::ParameterList& solverparams,  ///< parameter list for scalar transport solver
    const std::string& disname,                  ///< name of scalar transport discretization
    const bool isale                             ///< ALE flag
)
{
  // call setup in base class
  ADAPTER::ScaTraFluidAleCouplingAlgorithm::Init(prbdyn, scatradyn, solverparams, disname, isale);

  // safety check
  if (!ScaTraField()->Discretization()->GetCondition("ScaTraFluxCalc"))
    dserror(
        "Scalar transport discretization must have boundary condition for flux calculation at FSI "
        "interface!");

  pseudotransient_ = (DRT::INPUT::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elch_params_,
                          "MOVINGBOUNDARY") == INPAR::ELCH::elch_mov_bndry_pseudo_transient);

  return;
}


/*----------------------------------------------------------------------*
| Setup                                                     rauch 08/16 |
*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Setup()
{
  // call init in base class
  ADAPTER::ScaTraFluidAleCouplingAlgorithm::Setup();

  // set pointers
  idispn_ = FluidField()->ExtractInterfaceVeln();
  idispnp_ = FluidField()->ExtractInterfaceVeln();
  iveln_ = FluidField()->ExtractInterfaceVeln();

  idispn_->PutScalar(0.0);
  idispnp_->PutScalar(0.0);
  iveln_->PutScalar(0.0);

  // calculate normal flux vector field only at FSICoupling boundaries (no output to file)
  if (pseudotransient_ or (theta_ < 0.999))
  {
    SolveScaTra();  // set-up trueresidual_
  }

  // transfer moving mesh data
  ScaTraField()->ApplyMeshMovement(AleField()->Dispnp(), 2);

  // initialize the multivector for all possible cases
  fluxn_ = ScaTraField()->CalcFluxAtBoundary(false);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::TimeLoop()
{
  // safety checks
  CheckIsInit();
  CheckIsSetup();

  // provide information about initial field (do not do for restarts!)
  if (Step() == 0)
  {
    FluidField()->StatisticsAndOutput();
    if (AlgoParameters().get<int>("RESTARTEVRY") != 0)
      FluidField()->DiscWriter()->WriteVector("idispn", idispnp_);
    AleField()->Output();
  }

  // prepare scatra field
  ScaTraField()->PrepareTimeLoop();

  if (not pseudotransient_)
  {
    // transfer convective velocity = fluid velocity - grid velocity
    ScaTraField()->SetVelocityField(FluidField()->ConvectiveVel(),  // = velnp - grid velocity
        FluidField()->Hist(), Teuchos::null, Teuchos::null, 1);
  }

  // transfer moving mesh data
  ScaTraField()->ApplyMeshMovement(AleField()->Dispnp(), 2);

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    PrepareTimeStep();

    Teuchos::RCP<Epetra_Vector> incr(FluidField()->ExtractInterfaceVeln());
    incr->PutScalar(0.0);
    double incnorm(0.0);
    int iter(0);
    bool stopiter(false);

    // ToDo
    // improve this convergence test
    // (better check increment of ivel ????, test relative value etc.)
    while (stopiter == false)  // do at least one step
    {
      iter++;

      /// compute interface displacement and velocity
      ComputeInterfaceVectors(idispnp_, iveln_);

      // save guessed value before solve
      incr->Update(1.0, *idispnp_, 0.0);

      // solve nonlinear Navier-Stokes system on a deforming mesh
      SolveFluidAle();

      // solve transport equations for ion concentrations and electric potential
      SolveScaTra();

      /// compute interface displacement and velocity
      ComputeInterfaceVectors(idispnp_, iveln_);

      // compare with value after solving
      incr->Update(-1.0, *idispnp_, 1.0);

      // compute L2 norm of increment
      incr->Norm2(&incnorm);

      if (Comm().MyPID() == 0)
      {
        std::cout << "After outer iteration " << iter << " of " << itmax_
                  << ":  ||idispnpinc|| = " << incnorm << std::endl;
      }
      if (incnorm < ittol_)
      {
        stopiter = true;
        if (Comm().MyPID() == 0) std::cout << "   || Outer iteration loop converged! ||\n\n\n";
      }
      if (iter == itmax_)
      {
        stopiter = true;
        if (Comm().MyPID() == 0)
          std::cout << "   || Maximum number of iterations reached: " << itmax_ << " ||\n\n\n";
      }
    }

    double normidsinp;
    idispnp_->Norm2(&normidsinp);
    std::cout << "norm of isdispnp = " << normidsinp << std::endl;

    // update all single field solvers
    Update();

    // compute error for problems with analytical solution
    ScaTraField()->EvaluateErrorComparedToAnalyticalSol();

    // write output to screen and files
    Output();

  }  // time loop

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::PrepareTimeStep()
{
  IncrementTimeAndStep();

  // screen output
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "*************************************************************************"
              << std::endl;
    std::cout << "  MOVING-BOUNDARY ALGORITHM FOR ELECTROCHEMISTRY  ---  STEP = " << std::setw(4)
              << Step() << "/" << std::setw(4) << NStep() << std::endl;
    std::cout << "*************************************************************************"
              << std::endl
              << std::endl;
  }

  FluidField()->PrepareTimeStep();
  AleField()->PrepareTimeStep();

  // prepare time step
  /* remark: initial velocity field has been transferred to scalar transport field in constructor of
   * ScaTraFluidCouplingMovingBoundaryAlgorithm (initialvelset_ == true). Time integration schemes,
   * such as the one-step-theta scheme, are thus initialized correctly.
   */
  ScaTraField()->PrepareTimeStep();

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::SolveFluidAle()
{
  // screen output
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "*********************" << std::endl;
    std::cout << "  FLUID-ALE SOLVER   " << std::endl;
    std::cout << "*********************" << std::endl;
  }

  // solve nonlinear Navier-Stokes system on a moving mesh
  FluidAleNonlinearSolve(idispnp_, iveln_, pseudotransient_);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::SolveScaTra()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "************************" << std::endl;
    std::cout << "       ELCH SOLVER      " << std::endl;
    std::cout << "************************" << std::endl;
  }

  switch (FluidField()->TimIntScheme())
  {
    case INPAR::FLUID::timeint_npgenalpha:
    case INPAR::FLUID::timeint_afgenalpha:
      dserror("ConvectiveVel() not implemented for Gen.Alpha versions");
      break;
    case INPAR::FLUID::timeint_one_step_theta:
    case INPAR::FLUID::timeint_bdf2:
    {
      if (not pseudotransient_)
      {
        // transfer convective velocity = fluid velocity - grid velocity
        ScaTraField()->SetVelocityField(FluidField()->ConvectiveVel(),  // = velnp - grid velocity
            FluidField()->Hist(), Teuchos::null, Teuchos::null, 1);
      }
    }
    break;
    default:
      dserror("Time integration scheme not supported");
      break;
  }

  // transfer moving mesh data
  ScaTraField()->ApplyMeshMovement(AleField()->Dispnp(), 2);

  // solve coupled electrochemistry equations
  ScaTraField()->Solve();

  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Update()
{
  FluidField()->Update();
  AleField()->Update();
  ScaTraField()->Update();

  // perform time shift of interface displacement
  idispn_->Update(1.0, *idispnp_, 0.0);
  // perform time shift of interface mass flux vectors
  fluxn_->Update(1.0, *fluxnp_, 0.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  FluidField()->StatisticsAndOutput();
  // additional vector needed for restarts:
  int uprestart = AlgoParameters().get<int>("RESTARTEVRY");
  if ((uprestart != 0) && (FluidField()->Step() % uprestart == 0))
  {
    FluidField()->DiscWriter()->WriteVector("idispn", idispnp_);
  }

#if 0
  //ToDo
  // for visualization only...
  Teuchos::RCP<Epetra_Vector> idispnpfull = LINALG::CreateVector(*(FluidField()->DofRowMap()),true);
  (FluidField()->Interface()).AddFSICondVector(idispnp_,idispnpfull);
  FluidField()->DiscWriter()->WriteVector("idispnfull",idispnpfull);
#endif

  // now the other physical fiels
  ScaTraField()->Output();
  AleField()->Output();

  return;
}


void ELCH::MovingBoundaryAlgorithm::ComputeInterfaceVectors(
    Teuchos::RCP<Epetra_Vector> idispnp, Teuchos::RCP<Epetra_Vector> iveln)
{
  // calculate normal flux vector field at FSI boundaries (no output to file)
  fluxnp_ = ScaTraField()->CalcFluxAtBoundary(false);

  // access discretizations
  Teuchos::RCP<DRT::Discretization> fluiddis = FluidField()->Discretization();
  Teuchos::RCP<DRT::Discretization> scatradis = ScaTraField()->Discretization();

  // no support for multiple reactions at the interface !
  // id of the reacting species
  int reactingspeciesid = 0;

  const Epetra_BlockMap& ivelmap = iveln->Map();

  // loop over all local nodes of fluid discretization
  for (int lnodeid = 0; lnodeid < fluiddis->NumMyRowNodes(); lnodeid++)
  {
    // Here we rely on the fact that the scatra discretization
    // is a clone of the fluid mesh. => a scatra node has the same
    // local (and global) ID as its corresponding fluid node!

    // get the processor's local fluid node with the same lnodeid
    DRT::Node* fluidlnode = fluiddis->lRowNode(lnodeid);
    // get the degrees of freedom associated with this fluid node
    std::vector<int> fluidnodedofs = fluiddis->Dof(0, fluidlnode);

    if (ivelmap.MyGID(fluidnodedofs[0]))  // is this GID (implies: node) relevant for iveln_?
    {
      // determine number of space dimensions (numdof - pressure dof)
      const int numdim = ((int)fluidnodedofs.size()) - 1;
      // number of dof per node in ScaTra
      int numscatradof = scatradis->NumDof(0, scatradis->lRowNode(lnodeid));

      std::vector<double> Values(numdim);
      for (int index = 0; index < numdim; ++index)
      {
        const int pos = lnodeid * numscatradof + reactingspeciesid;
        // interface growth has opposite direction of metal ion mass flow -> minus sign !!
        Values[index] = (-molarvolume_) * (theta_ * (((*fluxnp_)[index])[pos]) +
                                              (1.0 - theta_) * (((*fluxn_)[index])[pos]));
      }

      // now insert only the first numdim entries (pressure dof is not inserted!)
      int error = iveln_->ReplaceGlobalValues(numdim, &Values[0], &fluidnodedofs[0]);
      if (error > 0) dserror("Could not insert values into vector iveln_: error %d", error);
    }
  }

  // have to compute an approximate displacement from given interface velocity
  // id^{n+1} = id^{n} + \delta t vel_i
  idispnp->Update(1.0, *idispn_, 0.0);
  idispnp->Update(Dt(), *iveln_, 1.0);

  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::ReadRestart(int step)
{
  ScaTraFluidCouplingAlgorithm::ReadRestart(step);

  AleField()->ReadRestart(step);  // add reading of ALE restart data

  // finally read isdispn which was written to the fluid restart data
  IO::DiscretizationReader reader(FluidField()->Discretization(), step);
  reader.ReadVector(idispn_, "idispn");
  // read same result into vector isdispnp_ as a 'good guess'
  reader.ReadVector(idispnp_, "idispn");

  return;
}
