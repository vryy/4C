/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all ELCH algorithms with moving boundaries

\level 2
*/
/*----------------------------------------------------------------------*/


#include "4C_elch_moving_boundary_algorithm.hpp"

#include "4C_fluid_utils_mapextractor.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_elch.hpp"
#include "4C_io.hpp"
#include "4C_lib_discret.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_scatra_timint_elch.hpp"

FOUR_C_NAMESPACE_OPEN

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
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Init()
{
  // call setup in base class
  ADAPTER::ScaTraFluidAleCouplingAlgorithm::Init();

  // safety check
  if (!ScaTraField()->discretization()->GetCondition("ScaTraFluxCalc"))
  {
    FOUR_C_THROW(
        "Scalar transport discretization must have boundary condition for flux calculation at FSI "
        "interface!");
  }

  pseudotransient_ = (CORE::UTILS::IntegralValue<INPAR::ELCH::ElchMovingBoundary>(elch_params_,
                          "MOVINGBOUNDARY") == INPAR::ELCH::elch_mov_bndry_pseudo_transient);
}


/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::Setup()
{
  // call init in base class
  ADAPTER::ScaTraFluidAleCouplingAlgorithm::Setup();

  // set pointers
  idispn_ = fluid_field()->extract_interface_veln();
  idispnp_ = fluid_field()->extract_interface_veln();
  iveln_ = fluid_field()->extract_interface_veln();

  idispn_->PutScalar(0.0);
  idispnp_->PutScalar(0.0);
  iveln_->PutScalar(0.0);

  // calculate normal flux vector field only at FSICoupling boundaries (no output to file)
  if (pseudotransient_ or (theta_ < 0.999))
  {
    solve_sca_tra();  // set-up trueresidual_
  }

  // transfer moving mesh data
  ScaTraField()->ApplyMeshMovement(ale_field()->Dispnp());

  // initialize the multivector for all possible cases
  fluxn_ = ScaTraField()->CalcFluxAtBoundary(false);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::TimeLoop()
{
  // safety checks
  check_is_init();
  check_is_setup();

  // provide information about initial field (do not do for restarts!)
  if (Step() == 0)
  {
    fluid_field()->StatisticsAndOutput();
    if (algo_parameters().get<int>("RESTARTEVRY") != 0)
      fluid_field()->DiscWriter()->WriteVector("idispn", idispnp_);
    ale_field()->Output();
  }

  // prepare scatra field
  ScaTraField()->prepare_time_loop();

  if (not pseudotransient_)
  {
    // transfer convective velocity = fluid velocity - grid velocity
    ScaTraField()->set_velocity_field(fluid_field()->ConvectiveVel(),  // = velnp - grid velocity
        fluid_field()->Hist(), Teuchos::null, Teuchos::null);
  }

  // transfer moving mesh data
  ScaTraField()->ApplyMeshMovement(ale_field()->Dispnp());

  // time loop
  while (NotFinished())
  {
    // prepare next time step
    prepare_time_step();

    auto incr = fluid_field()->extract_interface_veln();
    incr->PutScalar(0.0);
    double incnorm = 0.0;
    int iter = 0;
    bool stopiter = false;

    // ToDo
    // improve this convergence test
    // (better check increment of ivel ????, test relative value etc.)
    while (!stopiter)  // do at least one step
    {
      iter++;

      /// compute interface displacement and velocity
      compute_interface_vectors(idispnp_, iveln_);

      // save guessed value before solve
      incr->Update(1.0, *idispnp_, 0.0);

      // solve nonlinear Navier-Stokes system on a deforming mesh
      solve_fluid_ale();

      // solve transport equations for ion concentrations and electric potential
      solve_sca_tra();

      /// compute interface displacement and velocity
      compute_interface_vectors(idispnp_, iveln_);

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
    update();

    // compute error for problems with analytical solution
    ScaTraField()->evaluate_error_compared_to_analytical_sol();

    // write output to screen and files
    output();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::prepare_time_step()
{
  increment_time_and_step();

  // screen output
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "*************************************************************************"
              << std::endl;
    std::cout << "  MOVING-BOUNDARY ALGORITHM FOR ELECTROCHEMISTRY  ---  STEP = " << std::setw(4)
              << Step() << "/" << std::setw(4) << n_step() << std::endl;
    std::cout << "*************************************************************************"
              << std::endl
              << std::endl;
  }

  fluid_field()->prepare_time_step();
  ale_field()->prepare_time_step();

  // prepare time step
  /* remark: initial velocity field has been transferred to scalar transport field in constructor of
   * ScaTraFluidCouplingMovingBoundaryAlgorithm (initialvelset_ == true). Time integration schemes,
   * such as the one-step-theta scheme, are thus initialized correctly.
   */
  ScaTraField()->prepare_time_step();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::solve_fluid_ale()
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
  fluid_ale_nonlinear_solve(idispnp_, iveln_, pseudotransient_);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::solve_sca_tra()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "************************" << std::endl;
    std::cout << "       ELCH SOLVER      " << std::endl;
    std::cout << "************************" << std::endl;
  }

  switch (fluid_field()->TimIntScheme())
  {
    case INPAR::FLUID::timeint_npgenalpha:
    case INPAR::FLUID::timeint_afgenalpha:
      FOUR_C_THROW("ConvectiveVel() not implemented for Gen.Alpha versions");
      break;
    case INPAR::FLUID::timeint_one_step_theta:
    case INPAR::FLUID::timeint_bdf2:
    {
      if (not pseudotransient_)
      {
        // transfer convective velocity = fluid velocity - grid velocity
        ScaTraField()->set_velocity_field(
            fluid_field()->ConvectiveVel(),  // = velnp - grid velocity
            fluid_field()->Hist(), Teuchos::null, Teuchos::null);
      }
    }
    break;
    default:
      FOUR_C_THROW("Time integration scheme not supported");
      break;
  }

  // transfer moving mesh data
  ScaTraField()->ApplyMeshMovement(ale_field()->Dispnp());

  // solve coupled electrochemistry equations
  ScaTraField()->Solve();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::update()
{
  fluid_field()->Update();
  ale_field()->Update();
  ScaTraField()->Update();

  // perform time shift of interface displacement
  idispn_->Update(1.0, *idispnp_, 0.0);
  // perform time shift of interface mass flux vectors
  fluxn_->Update(1.0, *fluxnp_, 0.0);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::output()
{
  // Note: The order is important here! In here control file entries are
  // written. And these entries define the order in which the filters handle
  // the discretizations, which in turn defines the dof number ordering of the
  // discretizations.
  fluid_field()->StatisticsAndOutput();
  // additional vector needed for restarts:
  int uprestart = algo_parameters().get<int>("RESTARTEVRY");
  if ((uprestart != 0) && (fluid_field()->Step() % uprestart == 0))
  {
    fluid_field()->DiscWriter()->WriteVector("idispn", idispnp_);
  }

  // now the other physical fiels
  ScaTraField()->check_and_write_output_and_restart();
  ale_field()->Output();
}


void ELCH::MovingBoundaryAlgorithm::compute_interface_vectors(
    Teuchos::RCP<Epetra_Vector> idispnp, Teuchos::RCP<Epetra_Vector> iveln)
{
  // calculate normal flux vector field at FSI boundaries (no output to file)
  fluxnp_ = ScaTraField()->CalcFluxAtBoundary(false);

  // access discretizations
  Teuchos::RCP<DRT::Discretization> fluiddis = fluid_field()->discretization();
  Teuchos::RCP<DRT::Discretization> scatradis = ScaTraField()->discretization();

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
    CORE::Nodes::Node* fluidlnode = fluiddis->lRowNode(lnodeid);
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
      int error = iveln_->ReplaceGlobalValues(numdim, Values.data(), fluidnodedofs.data());
      if (error > 0) FOUR_C_THROW("Could not insert values into vector iveln_: error %d", error);
    }
  }

  // have to compute an approximate displacement from given interface velocity
  // id^{n+1} = id^{n} + \delta t vel_i
  idispnp->Update(1.0, *idispn_, 0.0);
  idispnp->Update(Dt(), *iveln_, 1.0);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::read_restart(int step)
{
  ScaTraFluidCouplingAlgorithm::read_restart(step);

  ale_field()->read_restart(step);  // add reading of ALE restart data

  // finally read isdispn which was written to the fluid restart data
  IO::DiscretizationReader reader(
      fluid_field()->discretization(), GLOBAL::Problem::Instance()->InputControlFile(), step);
  reader.ReadVector(idispn_, "idispn");
  // read same result into vector isdispnp_ as a 'good guess'
  reader.ReadVector(idispnp_, "idispn");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ELCH::MovingBoundaryAlgorithm::TestResults()
{
  auto* problem = GLOBAL::Problem::Instance();
  problem->AddFieldTest(fluid_field()->CreateFieldTest());
  problem->AddFieldTest(ale_field()->CreateFieldTest());
  problem->AddFieldTest(ScaTraField()->create_sca_tra_field_test());
  problem->TestAll(ScaTraField()->discretization()->Comm());
}
FOUR_C_NAMESPACE_CLOSE
