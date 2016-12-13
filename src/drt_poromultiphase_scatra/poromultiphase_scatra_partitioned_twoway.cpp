/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_partitioned_twoway.cpp

 \brief two-way coupled partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/




#include "poromultiphase_scatra_partitioned_twoway.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_poromultiphase/poromultiphase_base.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_scatra/scatra_timint_implicit.H"

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::PoroMultiPhaseScaTraPartitionedTwoWay(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams):
    PoroMultiPhaseScaTraPartitioned(comm, globaltimeparams),
    scaincnp_(Teuchos::null),
    structincnp_(Teuchos::null),
    fluidincnp_(Teuchos::null),
    itmax_(-1),
    ittol_(-1)
{

}

/*----------------------------------------------------------------------*
 | initialization                                            vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::Init(
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& poroparams,
    const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams,
    const Teuchos::ParameterList& scatraparams,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    const std::string& scatra_disname,
    bool isale,
    int nds_disp,
    int nds_vel,
    int nds_solidpressure,
    int ndsporofluid_scatra)
{
  //call base class
  POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitioned::Init(
      globaltimeparams,
      algoparams,
      poroparams,
      structparams,
      fluidparams,
      scatraparams,
      struct_disname,
      fluid_disname,
      scatra_disname,
      isale,
      nds_disp,
      nds_vel,
      nds_solidpressure,
      ndsporofluid_scatra);

  // read input variables
  itmax_ = algoparams.get<int>("ITEMAX");
  ittol_ = algoparams.get<double>("TOLINC_GLOBAL");

  // initialize increment vectors
  scaincnp_ = Teuchos::rcp(new Epetra_Vector(*(scatra_->ScaTraField()->Phinp())));
  structincnp_ = (Teuchos::rcp(new Epetra_Vector(*(poromulti_->StructDispnp()))));
  fluidincnp_ = (Teuchos::rcp(new Epetra_Vector(*(poromulti_->FluidPhinp()))));

}

/*----------------------------------------------------------------------*
 | time loop                                                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::Timeloop()
{
  PrepareTimeLoop();

  while (NotFinished())
  {
    PrepareTimeStep();

    Solve();

    UpdateAndOutput();

  }
  return;

}

/*----------------------------------------------------------------------*
 | prepare one time step                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::PrepareTimeStep(bool printheader)
{
  // the global control routine has its own time_ and step_ variables, as well as the single fields
  // keep them in sync!
  IncrementTimeAndStep();

  if(printheader)
    PrintHeader();

  SetPoroSolution();
  scatra_->ScaTraField()->PrepareTimeStep();
  // set structure-based scalar transport values
  SetScatraSolution();

  poromulti_-> PrepareTimeStep();
  SetPoroSolution();

  return;
}

/*----------------------------------------------------------------------*
 | prepare the time loop                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::PrepareTimeLoop()
{
  // set structure-based scalar transport values
  SetScatraSolution();
  poromulti_->PrepareTimeLoop();
  // initial output for scatra field
  SetPoroSolution();
  scatra_->ScaTraField()->Output();
  return;
}

/*----------------------------------------------------------------------*
 | setup the system if necessary                             vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::SetupSystem()
{
  poromulti_->SetupSystem();
  return;
}

/*----------------------------------------------------------------------*
 | update fields and output results                         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::UpdateAndOutput()
{
  poromulti_-> UpdateAndOutput();

  scatra_->ScaTraField()->Update();
  scatra_->ScaTraField()->EvaluateErrorComparedToAnalyticalSol();
  scatra_->ScaTraField()->Output();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::Solve()
{
  int  itnum = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // store scalar from first solution for convergence check (like in
    // elch_algorithm: use current values)
    scaincnp_->Update(1.0,*scatra_->ScaTraField()->Phinp(),0.0);
    structincnp_->Update(1.0,*poromulti_->StructDispnp(),0.0);
    fluidincnp_->Update(1.0,*poromulti_->FluidPhinp(),0.0);

    // set structure-based scalar transport values
    SetScatraSolution();

    // solve structural system
    DoPoroStep();

    // set mesh displacement and velocity fields
    SetPoroSolution();

    // solve scalar transport equation
    DoScatraStep();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}


/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::DoPoroStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n POROUS MEDIUM SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  poromulti_-> TimeStep();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_->ScaTraField()->Solve();

}

/*----------------------------------------------------------------------*
 | convergence check for both fields (scatra & poro) (copied from tsi)
 *----------------------------------------------------------------------*/
bool POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::ConvergenceCheck(int itnum)
{

  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  //    | scalar increment |_2
  //  -------------------------------- < Tolerance
  //     | scalar+1 |_2

  // variables to save different L2 - Norms
  // define L2-norm of incremental scalar and scalar
  double scaincnorm_L2(0.0);
  double scanorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);
  double fluidincnorm_L2(0.0);
  double fluidnorm_L2(0.0);

  // build the current scalar increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  scaincnp_->Update(1.0,*(scatra_->ScaTraField()->Phinp()),-1.0);
  structincnp_->Update(1.0,*(poromulti_->StructDispnp()),-1.0);
  fluidincnp_->Update(1.0,*(poromulti_->FluidPhinp()),-1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  scatra_->ScaTraField()->Phinp()->Norm2(&scanorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  poromulti_->StructDispnp()->Norm2(&dispnorm_L2);
  fluidincnp_->Norm2(&fluidincnorm_L2);
  poromulti_->FluidPhinp()->Norm2(&fluidnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (scanorm_L2 < 1e-6) scanorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;
  if (fluidnorm_L2 < 1e-6) fluidnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"\n";
    std::cout<<"***********************************************************************************\n";
    std::cout<<"    OUTER ITERATION STEP    \n";
    std::cout<<"***********************************************************************************\n";
    printf("+--------------+------------------------+--------------------+--------------------+--------------------+\n");
    printf("|-  step/max  -|-  tol      [norm]     -|--  scalar-inc      --|--  disp-inc      --|--  fluid-inc      --|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         | %10.3E         | %10.3E         |",
         itnum,itmax_,ittol_,scaincnorm_L2/scanorm_L2,dispincnorm_L2/dispnorm_L2,fluidincnorm_L2/fluidnorm_L2);
    printf("\n");
    printf("+--------------+------------------------+--------------------+--------------------+--------------------+\n");
  }

  // converged
  if ((scaincnorm_L2/scanorm_L2 <= ittol_) and
      (dispincnorm_L2/dispnorm_L2 <= ittol_) and
      (fluidincnorm_L2/fluidnorm_L2 <= ittol_)
      )
  {
    stopnonliniter = true;
    if (Comm().MyPID()==0 )
    {
      printf("\n");
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                       |\n", itnum,itmax_);
      printf("+--------------+------------------------+--------------------+--------------------+\n");
    }
  }

  // warn if itemax is reached without convergence, but proceed to next
  // timestep
  if ((itnum==itmax_) and
       ( (scaincnorm_L2/scanorm_L2 > ittol_) or (dispincnorm_L2/dispnorm_L2 > ittol_) or (fluidincnorm_L2/fluidnorm_L2 > ittol_) )
     )
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) )
    {
      printf("|     >>>>>> not converged in itemax steps!                                       |\n");
      printf("+--------------+------------------------+--------------------+--------------------+\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}


