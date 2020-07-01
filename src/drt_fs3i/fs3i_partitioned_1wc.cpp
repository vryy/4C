/*----------------------------------------------------------------------*/
/*! \file
\brief Algorithmic routines for partitioned solution approaches to
       fluid-structure-scalar-scalar interaction (FS3I) specifically
       related to one-way-coupled problem configurations

\level 2



*----------------------------------------------------------------------*/


#include "fs3i_partitioned_1wc.H"

#include "../drt_fsi/fsi_monolithic.H"
#include "../drt_scatra/scatra_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"
#include "../drt_inpar/inpar_scatra.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_fluid_fsi.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I_1WC::PartFS3I_1WC(const Epetra_Comm& comm) : PartFS3I(comm) {}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_1WC::Init()
{
  FS3I::PartFS3I::Init();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_1WC::Setup()
{
  FS3I::PartFS3I::Setup();
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_1WC::Timeloop()
{
  CheckIsInit();
  CheckIsSetup();

  // prepare time loop
  fsi_->PrepareTimeloop();
  SetFSISolution();

  // calculate inital time derivative, when restart was done from a part. FSI simulation
  if (DRT::Problem::Instance()->Restart() and
      DRT::INPUT::IntegralValue<int>(
          DRT::Problem::Instance()->FS3IDynamicParams(), "RESTART_FROM_PART_FSI"))
  {
    scatravec_[0]->ScaTraField()->PrepareFirstTimeStep();
    scatravec_[1]->ScaTraField()->PrepareFirstTimeStep();
  }

  // output of initial state
  fsi_->PrepareOutput();
  fsi_->Output();
  ScatraOutput();

  while (NotFinished())
  {
    IncrementTimeAndStep();
    SetStructScatraSolution();
    DoFSIStep();
    SetFSISolution();
    DoScatraStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_1WC::DoFSIStep()
{
  fsi_->PrepareTimeStep();
  fsi_->TimeStep(fsi_);
  fsi_->PrepareOutput();
  fsi_->Update();
  fsi_->Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_1WC::DoScatraStep()
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

  PrepareTimeStep();

  while (stopnonliniter == false)
  {
    ScatraEvaluateSolveIterUpdate();
    itnum++;
    if (ScatraConvergenceCheck(itnum)) break;
  }

  UpdateScatraFields();
  ScatraOutput();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_1WC::PrepareTimeStep()
{
  CheckIsInit();
  CheckIsSetup();

  // set mesh displacement field for present time step
  SetMeshDisp();

  // set velocity fields from fluid and structure solution
  // for present time step
  SetVelocityFields();

  // prepare time step for both fluid- and structure-based scatra field
  for (unsigned i = 0; i < scatravec_.size(); ++i)
  {
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatra = scatravec_[i];
    scatra->ScaTraField()->PrepareTimeStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I_1WC::ScatraConvergenceCheck(const int itnum)
{
  const Teuchos::ParameterList& fs3idyn = DRT::Problem::Instance()->FS3IDynamicParams();
  INPAR::SCATRA::SolverType scatra_solvtype =
      DRT::INPUT::IntegralValue<INPAR::SCATRA::SolverType>(fs3idyn, "SCATRA_SOLVERTYPE");

  double conresnorm(0.0);
  scatrarhs_->Norm2(&conresnorm);
  double incconnorm(0.0);
  scatraincrement_->Norm2(&incconnorm);

  switch (scatra_solvtype)
  {
    case INPAR::SCATRA::solvertype_linear_incremental:
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
    case INPAR::SCATRA::solvertype_nonlinear:
    {
      // some input parameters for the scatra fields
      const Teuchos::ParameterList& scatradyn =
          DRT::Problem::Instance()->ScalarTransportDynamicParams();
      const int itemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
      const double ittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
      const double abstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");

      double connorm(0.0);
      // set up vector of absolute concentrations
      Teuchos::RCP<Epetra_Vector> con = Teuchos::rcp(new Epetra_Vector(scatraincrement_->Map()));
      Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField()->Phinp();
      Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField()->Phinp();
      SetupCoupledScatraVector(con, scatra1, scatra2);
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
      dserror("Illegal ScaTra solvertype in FS3I");
      break;
  }
  return false;
}
