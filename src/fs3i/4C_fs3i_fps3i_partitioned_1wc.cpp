/*----------------------------------------------------------------------*/
/*! \file
\brief H-file associated with algorithmic routines for partitioned
       solution approaches to fluid-porous-structure-scalar-scalar interaction
       (FPS3I) specifically related to one-way-coupled problem
       configurations

\level 3


*----------------------------------------------------------------------*/
#include "4C_fs3i_fps3i_partitioned_1wc.hpp"

#include "4C_adapter_fld_poro.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_fem_condition_selector.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fluid_implicit_integration.hpp"
#include "4C_fluid_result_test.hpp"
#include "4C_fluid_utils.hpp"
#include "4C_fpsi_monolithic.hpp"
#include "4C_fpsi_monolithic_plain.hpp"
#include "4C_fpsi_utils.hpp"
#include "4C_fsi_monolithic.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_scatra.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_io_control.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_linear_solver_method_linalg.hpp"
#include "4C_poroelast_monolithic.hpp"
#include "4C_scatra_algorithm.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |  Constructor                                           hemmler 07/14 |
 *----------------------------------------------------------------------*/
FS3I::PartFpS3I1Wc::PartFpS3I1Wc(const Epetra_Comm& comm) : PartFPS3I(comm)
{
  // keep constructor empty
  return;
}


/*----------------------------------------------------------------------*
 |  Init                                                    rauch 09/16 |
 *----------------------------------------------------------------------*/
void FS3I::PartFpS3I1Wc::Init()
{
  FS3I::PartFPS3I::Init();



  return;
}


/*----------------------------------------------------------------------*
 |  Setup                                                   rauch 09/16 |
 *----------------------------------------------------------------------*/
void FS3I::PartFpS3I1Wc::setup()
{
  FS3I::PartFPS3I::setup();

  // add proxy of fluid degrees of freedom to scatra discretization
  if (scatravec_[0]->ScaTraField()->discretization()->AddDofSet(
          fpsi_->fluid_field()->discretization()->GetDofSetProxy()) != 1)
    FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");

  // check if scatra field has 2 discretizations, so that coupling is possible
  if (scatravec_[1]->ScaTraField()->discretization()->AddDofSet(
          fpsi_->poro_field()->structure_field()->discretization()->GetDofSetProxy()) != 1)
    FOUR_C_THROW("unexpected dof sets in structure field");
  // check if scatra field has 3 discretizations, so that coupling is possible
  if (scatravec_[1]->ScaTraField()->discretization()->AddDofSet(
          fpsi_->poro_field()->fluid_field()->discretization()->GetDofSetProxy()) != 2)
    FOUR_C_THROW("unexpected dof sets in structure field");

  // Note: in the scatra fields we have now the following dof-sets:
  // fluidscatra dofset 0: fluidscatra dofset
  // fluidscatra dofset 1: fluid dofset
  // structscatra dofset 0: structscatra dofset
  // structscatra dofset 1: structure dofset
  // structscatra dofset 2: porofluid dofset

  return;
}


/*----------------------------------------------------------------------*
 |  Timeloop                                              hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFpS3I1Wc::Timeloop()
{
  check_is_init();
  check_is_setup();

  // write FPSI solution into scatra field
  SetFPSISolution();

  // output of initial state
  ScatraOutput();

  fpsi_->PrepareTimeloop();

  while (NotFinished())
  {
    increment_time_and_step();

    DoFPSIStep();  // TODO: One could think about skipping the very costly FSI/FPSI calculation for
                   // the case that it is stationary at some point (Thon)
    SetFPSISolution();  // write FPSI solution into scatra field
    do_scatra_step();
  }
}


/*----------------------------------------------------------------------*
 |  FPSI step                                             hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFpS3I1Wc::DoFPSIStep()
{
  fpsi_->prepare_time_step();
  fpsi_->SetupNewton();
  fpsi_->TimeStep();

  constexpr bool force_prepare = false;
  fpsi_->prepare_output(force_prepare);
  fpsi_->update();
  fpsi_->output();
}


/*----------------------------------------------------------------------*
 |  Scatra step                                           hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFpS3I1Wc::do_scatra_step()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n SCALAR TRANSPORT SOLVER \n***********************\n";
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


/*----------------------------------------------------------------------*
 |  Prepare time step                                     hemmler 07/14 |
 *----------------------------------------------------------------------*/
void FS3I::PartFpS3I1Wc::prepare_time_step()
{
  check_is_init();
  check_is_setup();

  // prepare time step for both fluid- and poro-based scatra field
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->prepare_time_step();
  }
}


/*----------------------------------------------------------------------*
 |  Check of convergence of scatra solver                 hemmler 07/14 |
 *----------------------------------------------------------------------*/
bool FS3I::PartFpS3I1Wc::scatra_convergence_check(const int itnum)
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
      FOUR_C_THROW("Illegal ScaTra solvertype in FPS3I");
      break;
  }
  return false;
}

FOUR_C_NAMESPACE_CLOSE
