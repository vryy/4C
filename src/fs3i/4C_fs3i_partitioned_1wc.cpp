/*----------------------------------------------------------------------*/
/*! \file
\brief Algorithmic routines for partitioned solution approaches to
       fluid-structure-scalar-scalar interaction (FS3I) specifically
       related to one-way-coupled problem configurations

\level 2



*----------------------------------------------------------------------*/


#include "4C_fs3i_partitioned_1wc.hpp"

#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fsi_monolithic.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I1Wc::PartFS3I1Wc(const Epetra_Comm& comm) : PartFS3I(comm) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I1Wc::init()
{
  FS3I::PartFS3I::init();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I1Wc::setup()
{
  FS3I::PartFS3I::setup();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I1Wc::Timeloop()
{
  check_is_init();
  check_is_setup();

  // prepare time loop
  fsi_->PrepareTimeloop();
  SetFSISolution();

  // calculate inital time derivative, when restart was done from a part. FSI simulation
  if (Global::Problem::Instance()->restart() and
      Core::UTILS::IntegralValue<int>(
          Global::Problem::Instance()->FS3IDynamicParams(), "RESTART_FROM_PART_FSI"))
  {
    scatravec_[0]->ScaTraField()->prepare_first_time_step();
    scatravec_[1]->ScaTraField()->prepare_first_time_step();
  }

  // output of initial state
  constexpr bool force_prepare = true;
  fsi_->prepare_output(force_prepare);
  fsi_->output();
  ScatraOutput();

  while (NotFinished())
  {
    increment_time_and_step();
    set_struct_scatra_solution();
    DoFSIStep();
    SetFSISolution();
    do_scatra_step();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I1Wc::DoFSIStep()
{
  fsi_->prepare_time_step();
  fsi_->TimeStep(fsi_);
  constexpr bool force_prepare = false;
  fsi_->prepare_output(force_prepare);
  fsi_->update();
  fsi_->output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I1Wc::do_scatra_step()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n GAS TRANSPORT SOLVER \n***********************\n"
              << std::endl;
    std::cout << "+- step/max -+- abs-res-tol [norm] -+-- scal-res --+- rel-inc-tol [norm] -+-- "
                 "scal-inc --+"
              << std::endl;
  }

  // first scatra field is associated with fluid, second scatra field is
  // associated with structure

  bool stopnonliniter = false;
  int itnum = 0;

  prepare_time_step();

  while (stopnonliniter == false)
  {
    scatra_evaluate_solve_iter_update();
    itnum++;
    if (scatra_convergence_check(itnum)) break;
  }

  UpdateScatraFields();
  ScatraOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I1Wc::prepare_time_step()
{
  check_is_init();
  check_is_setup();

  // set mesh displacement field for present time step
  set_mesh_disp();

  // set velocity fields from fluid and structure solution
  // for present time step
  set_velocity_fields();

  // prepare time step for both fluid- and structure-based scatra field
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->prepare_time_step();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I1Wc::scatra_convergence_check(const int itnum)
{
  const Teuchos::ParameterList& fs3idyn = Global::Problem::Instance()->FS3IDynamicParams();
  Inpar::ScaTra::SolverType scatra_solvtype =
      Core::UTILS::IntegralValue<Inpar::ScaTra::SolverType>(fs3idyn, "SCATRA_SOLVERTYPE");

  double conresnorm(0.0);
  scatrarhs_->Norm2(&conresnorm);
  double incconnorm(0.0);
  scatraincrement_->Norm2(&incconnorm);

  switch (scatra_solvtype)
  {
    case Inpar::ScaTra::solvertype_linear_incremental:
    {
      // print the screen info
      if (Comm().MyPID() == 0)
      {
        printf("\n+-------------------+-------------------+\n");
        printf("| norm of residual  | norm of increment |\n");
        printf("+-------------------+-------------------+\n");
        printf("|    %10.3E     |    %10.3E     |\n", conresnorm, incconnorm);
        printf("+-------------------+-------------------+\n\n");
      }
      return true;
    }
    break;
    case Inpar::ScaTra::solvertype_nonlinear:
    {
      // some input parameters for the scatra fields
      const Teuchos::ParameterList& scatradyn =
          Global::Problem::Instance()->scalar_transport_dynamic_params();
      const int itemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
      const double ittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
      const double abstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");

      double connorm(0.0);
      // set up vector of absolute concentrations
      Teuchos::RCP<Epetra_Vector> con = Teuchos::rcp(new Epetra_Vector(scatraincrement_->Map()));
      Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Phinp();
      Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Phinp();
      setup_coupled_scatra_vector(con, scatra1, scatra2);
      con->Norm2(&connorm);

      // care for the case that nothing really happens in the concentration field
      if (connorm < 1e-5) connorm = 1.0;

      // print the screen info
      if (Comm().MyPID() == 0)
      {
        printf("|  %3d/%3d   |   %10.3E [L_2 ]  | %10.3E   |   %10.3E [L_2 ]  | %10.3E   |\n",
            itnum, itemax, abstolres, conresnorm, ittol, incconnorm / connorm);
      }

      // this is the convergence check
      // We always require at least one solve. We test the L_2-norm of the
      // current residual. Norm of residual is just printed for information
      if (conresnorm <= abstolres and incconnorm / connorm <= ittol)
      {
        if (Comm().MyPID() == 0)
        {
          // print 'finish line'
          printf(
              "+------------+----------------------+--------------+----------------------+---------"
              "-----+\n\n");
        }
        return true;
      }
      // warn if itemax is reached without convergence, but proceed to
      // next timestep...
      else if (itnum == itemax)
      {
        if (Comm().MyPID() == 0)
        {
          printf("+---------------------------------------------------------------+\n");
          printf("|            >>>>>> not converged in itemax steps!              |\n");
          printf("+---------------------------------------------------------------+\n");
        }
        // yes, we stop the iteration
        return true;
      }
      else
        return false;
    }
    break;
    default:
      FOUR_C_THROW("Illegal ScaTra solvertype in FS3I");
      break;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
