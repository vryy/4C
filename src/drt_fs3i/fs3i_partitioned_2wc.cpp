#ifdef CCADISCRET

#include "fs3i_partitioned_2wc.H"

#include "../drt_fsi/fsi_monolithic_nox.H"
#include "../drt_scatra/passive_scatra_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FS3I::PartFS3I_2WC::PartFS3I_2WC(const Epetra_Comm& comm)
  :PartFS3I(comm)
{
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FS3I::PartFS3I_2WC::Timeloop()
{
  // output of initial state for scalars
  ScatraOutput();

  //fsi_->PrepareTimeloop();

  while (fsi_->NotFinished())
  {
    //DoFsiStep();
    SetFSISolution();
    DoScatraStep();
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool FS3I::PartFS3I_2WC::ScatraConvergenceCheck(int itnum)
{
#ifdef PARALLEL
  const Epetra_Comm& comm = scatravec_[0]->ScaTraField().Discretization()->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // define flags for convergence check for scatra fields
  bool scatra1stopnonliniter = false;
  bool scatra2stopnonliniter = false;

  // get tolerance and maximum number of iterations for convergence check
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  const int    itmax = scatradyn.sublist("NONLINEAR").get<int>("ITEMAX");
  const double ittol = scatradyn.sublist("NONLINEAR").get<double>("CONVTOL");

  // convergence check of scatra fields
  if (comm.MyPID() == 0)
  {
    cout<<"\n****************************************\n  CONVERGENCE CHECK FOR ITERATION STEP\n****************************************\n";
    cout<<"\n****************************************\n   FLUID-BASED SCALAR TRANSPORT CHECK\n****************************************\n";
  }
  scatra1stopnonliniter = scatravec_[0]->ScaTraField().ConvergenceCheck(itnum,itmax,ittol);

  if (comm.MyPID() == 0) cout<<"\n****************************************\n STRUCTURE-BASED SCALAR TRANSPORT CHECK\n****************************************\n";
  scatra2stopnonliniter = scatravec_[1]->ScaTraField().ConvergenceCheck(itnum,itmax,ittol);

  if (scatra1stopnonliniter == true and scatra2stopnonliniter == true) return true;
  else                                                                 return false;
}

#endif
