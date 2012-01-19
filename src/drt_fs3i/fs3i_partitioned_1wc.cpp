#ifdef CCADISCRET

#include "fs3i_partitioned_1wc.H"

#include "../drt_fsi/fsi_monolithic_nox.H"
#include "../drt_scatra/passive_scatra_algorithm.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I_1WC::PartFS3I_1WC(const Epetra_Comm& comm)
  :PartFS3I(comm)
{
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_1WC::Timeloop()
{
  // output of initial state
  ScatraOutput();

  fsi_->PrepareTimeloop();

  while (fsi_->NotFinished())
  {
    DoFsiStep();
    SetFSISolution();
    DoScatraStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I_1WC::ScatraConvergenceCheck(const int itnum)
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  int myrank = comm.MyPID();

  // some input parameters for the scatra fields
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const int itemax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
  const double ittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");
  const double abstolres = scatradyn.sublist("NONLINEAR").get<double>("ABSTOLRES");

  //----------------------------------------------------- compute norms
  double conresnorm(0.0);
  scatrarhs_->Norm2(&conresnorm);
  // set up vector of absolute concentrations
  double connorm_L2(0.0);
  Teuchos::RCP<Epetra_Vector> con = rcp(new Epetra_Vector(scatraincrement_->Map()));
  Teuchos::RCP<const Epetra_Vector> scatra1 = scatravec_[0]->ScaTraField().Phinp();
  Teuchos::RCP<const Epetra_Vector> scatra2 = scatravec_[1]->ScaTraField().Phinp();
  SetupCoupledScatraVector(con,scatra1,scatra2);
  con->Norm2(&connorm_L2);

  // care for the case that nothing really happens in the concentration field
  if (connorm_L2 < 1e-5)
  {
    connorm_L2 = 1.0;
  }

  // absolute tolerance for deciding if residual is (already) zero
  // prevents additional solver calls that will not improve the residual anymore

  //-------------------------------------------------- output to screen
  /* special case of very first iteration step:
      - solution increment is not yet available
      - do not perform a solver call when the initial residuals are < EPS14*/
  if (itnum == 0)
  {
    if (myrank == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |\n",itnum,itemax,ittol,conresnorm);
    }
    // abort iteration, when there's nothing more to do
    if (conresnorm < abstolres)
    {
      // print 'finish line'
      if (myrank == 0)
      {
        printf("+------------+-------------------+--------------+\n");
      }
      return true;
    }
    else
      return false;
  }
  /* ordinary case later iteration steps:
     - solution increment can be printed
     - convergence check should be done*/
  else
  {
    // print the screen info
    if (myrank == 0)
    {
      printf("|  %3d/%3d   | %10.3E[L_2 ]  | %10.3E   |\n",itnum,itemax,ittol,conresnorm);
    }

    // this is the convergence check
    // We always require at least one solve. We test the L_2-norm of the
    // current residual. Norm of residual is just printed for information
    if (conresnorm <= ittol)
    {
      if (myrank == 0)
      {
        // print 'finish line'
        printf("+------------+-------------------+--------------+\n");
      }
      return true;
    }

    // abort iteration, when there's nothing more to do! -> more robustness
    else if (conresnorm < abstolres)
    {
      // print 'finish line'
      if (myrank == 0)
      {
        printf("+------------+-------------------+--------------+\n");
      }
      return true;
    }

    // warn if itemax is reached without convergence, but proceed to
    // next timestep...
    else if (itnum == itemax)
    {
      if (myrank == 0)
      {
        printf("+-----------------------------------------------+\n");
        printf("| >>>>>>> scatra not converged in itemax steps! |\n");
        printf("+-----------------------------------------------+\n");
      }
      // yes, we stop the iteration
      return true;
    }
    else
      return false;
  }
}

#endif
