/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_partitioned_twoway.cpp

 \brief two-way coupled solution algorithm
        for porous multiphase flow through elastic medium problems

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/

#include "poromultiphase_partitioned_twoway.H"

#include "../drt_adapter/ad_porofluidmultiphase_wrapper.H"
#include "../drt_adapter/ad_str_wrapper.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 | constructor                                              vuong 08/16 |
 *----------------------------------------------------------------------*/
POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::PoroMultiPhasePartitionedTwoWay(
    const Epetra_Comm& comm,
    const Teuchos::ParameterList& globaltimeparams):
    PoroMultiPhasePartitioned(comm, globaltimeparams),
    phiincnp_(Teuchos::null),
    dispincnp_(Teuchos::null),
    ittol_(0.0),
    itmax_(0),
    itnum_(0)
{

}

/*----------------------------------------------------------------------*
 | initialization                                            vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::Init(
    const Teuchos::ParameterList& globaltimeparams,
    const Teuchos::ParameterList& algoparams,
    const Teuchos::ParameterList& structparams,
    const Teuchos::ParameterList& fluidparams,
    const std::string& struct_disname,
    const std::string& fluid_disname,
    bool isale,
    int nds_disp,
    int nds_vel,
    int nds_solidpressure,
    int ndsporofluid_scatra)
{
  //call base class
  POROMULTIPHASE::PoroMultiPhasePartitioned::Init(
      globaltimeparams,
      algoparams,
      structparams,
      fluidparams,
      struct_disname,
      fluid_disname,
      isale,
      nds_disp,
      nds_vel,
      nds_solidpressure,
      ndsporofluid_scatra);

  // initialize increment vectors
  phiincnp_ = LINALG::CreateVector(*fluid_->DofRowMap(0),true);
  dispincnp_ = LINALG::CreateVector(*structure_->DofRowMap(0),true);

  // Get the parameters for the ConvergenceCheck
  itmax_ = algoparams.get<int>("ITEMAX");
  ittol_ = algoparams.get<double>("CONVTOL");
}

/*----------------------------------------------------------------------*
 | time loop                                                 vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::Timeloop()
{
  // prepare the loop
  PrepareTimeLoop();

  //time loop
  while (NotFinished())
  {
    PrepareTimeStep();

    TimeStep();

    UpdateAndOutput();

  }

  return;

}

/*----------------------------------------------------------------------*
 | prepare one time step                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();
  if(printheader)
    PrintHeader();

  SetSolidPressure(fluid_->SolidPressure());
  //NOTE: the predictor of the structure is called in here
  structure_-> PrepareTimeStep();

  SetStructSolution(structure_->Dispnp(),structure_->Velnp());
  fluid_->PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 | prepare the time loop                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::PrepareTimeLoop()
{
  //initial output
  structure_-> PrepareOutput();
  structure_-> Output();
  SetStructSolution(structure_->Dispnp(),structure_->Velnp());
  fluid_->PrepareTimeLoop();

  return;
}

/*----------------------------------------------------------------------*
 | setup the system if necessary                             vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::SetupSystem()
{
  return;
}

/*----------------------------------------------------------------------*
 | Outer Timeloop  without relaxation                        vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::OuterLoop()
{
  // reset counter
  itnum_ = 0;
  bool stopnonliniter = false;

  if (Comm().MyPID()==0)
  {
    std::cout<<"\n****************************************\n          OUTER ITERATION LOOP\n****************************************\n";
  }

  while (stopnonliniter==false)
  {
    // increment number of iteration
    itnum_++;

    // update the states to the last solutions obtained
    IterUpdateStates();

    // 1.) solve structural system
    DoStructStep();
    // 2.) set disp and vel states in scatra field
    SetStructSolution(structure_->Dispnp(),structure_->Velnp());

    // 1.) solve scalar transport equation
    DoFluidStep();

    // 2.) set solid pressure state in structure field
    SetSolidPressure(fluid_->SolidPressure());


    // check convergence for all fields
    // stop iteration loop if converged
    stopnonliniter = ConvergenceCheck(itnum_);
  }

  return;
}

/*----------------------------------------------------------------------*
 | convergence check for both fields (copied form tsi)       vuong 08/16 |
 *----------------------------------------------------------------------*/
bool POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::ConvergenceCheck(int itnum)
{

  // convergence check based on the scalar increment
  bool stopnonliniter = false;

  //    | scalar increment |_2
  //  -------------------------------- < Tolerance
  //    | scalar state n+1 |_2
  //
  // AND
  //
  //    | scalar increment |_2
  //  -------------------------------- < Tolerance , with n := global length of vector
  //        dt * sqrt(n)
  //
  // The same is checked for the structural displacements.
  //

  // variables to save different L2 - Norms
  // define L2-norm of increments
  double phiincnorm_L2(0.0);
  double phinorm_L2(0.0);
  double dispincnorm_L2(0.0);
  double dispnorm_L2(0.0);

  // build the current scalar increment Inc T^{i+1}
  // \f Delta T^{k+1} = Inc T^{k+1} = T^{k+1} - T^{k}  \f
  phiincnp_->Update(1.0,*(fluid_->Phinp()),-1.0);
  dispincnp_->Update(1.0,*(structure_->Dispnp()),-1.0);

  // build the L2-norm of the scalar increment and the scalar
  phiincnp_->Norm2(&phiincnorm_L2);
  fluid_->Phinp()->Norm2(&phinorm_L2);
  dispincnp_->Norm2(&dispincnorm_L2);
  structure_->Dispnp()->Norm2(&dispnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (phinorm_L2 < 1e-6) phinorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"\n";
    std::cout<<"***********************************************************************************\n";
    std::cout<<"    OUTER ITERATION STEP    \n";
    std::cout<<"***********************************************************************************\n";
    printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
    printf("|-  step/max  -|-  tol      [norm]  -|-  scalar-inc  -|-  disp-inc      -|-  scalar-rel-inc  -|-  disp-rel-inc  -|\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]   |  %10.3E    |  %10.3E      |  %10.3E        |  %10.3E      |",
         itnum,itmax_,ittol_,phiincnorm_L2/Dt()/sqrt(phiincnp_->GlobalLength()),dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()),phiincnorm_L2/phinorm_L2,dispincnorm_L2/dispnorm_L2);
    printf("\n");
    printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
  }

  // converged
  if ( ((phiincnorm_L2/phinorm_L2) <= ittol_) and ((dispincnorm_L2/dispnorm_L2) <= ittol_) and ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))<=ittol_) and ((phiincnorm_L2/Dt()/sqrt(phiincnp_->GlobalLength()))<=ittol_) )
  {
    stopnonliniter = true;
    if (Comm().MyPID()==0 )
    {
      printf("|  Outer Iteration loop converged after iteration %3d/%3d !                                                      |\n", itnum,itmax_);
      printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
    }
  }

  // stop if itemax is reached without convergence
  // timestep
  if ( (itnum==itmax_) and (((phiincnorm_L2/phinorm_L2) > ittol_) or ((dispincnorm_L2/dispnorm_L2) > ittol_) or ((dispincnorm_L2/Dt()/sqrt(dispincnp_->GlobalLength()))>ittol_) or (phiincnorm_L2/Dt()/sqrt(phiincnp_->GlobalLength()))>ittol_ ) )
  {
    stopnonliniter = true;
    if ((Comm().MyPID()==0) )
    {
      printf("|     >>>>>> not converged in itemax steps!                                                                      |\n");
      printf("+--------------+---------------------+----------------+------------------+--------------------+------------------+\n");
      printf("\n");
      printf("\n");
    }
    dserror("The partitioned solver did not converge in ITEMAX steps!");
  }

  return stopnonliniter;
}

/*----------------------------------------------------------------------*
 | Solve structure filed                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::DoStructStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  structure_-> Solve();

  return;
}

/*----------------------------------------------------------------------*
 | Solve poro fluid field                                    vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::DoFluidStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n****************************\n  PORO MULTIPHASE FLUID SOLVER \n****************************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  fluid_->Solve();

  return;
}

/*----------------------------------------------------------------------*
 | update fields and output results                         vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::UpdateAndOutput()
{
  // prepare the output
  structure_-> PrepareOutput();

  // update single fields
  structure_-> Update();
  fluid_->Update();

  // evaluate error if desired
  fluid_->EvaluateErrorComparedToAnalyticalSol();

  // output single fields
  structure_-> Output();
  fluid_->Output();
}

/*----------------------------------------------------------------------*
 | update the current states in every iteration             vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASE::PoroMultiPhasePartitionedTwoWay::IterUpdateStates()
{
  // store last solutions (current states).
  // will be compared in ConvergenceCheck to the solutions,
  // obtained from the next Struct and Scatra steps.
  phiincnp_ ->Update(1.0,*fluid_->Phinp(),0.0);
  dispincnp_->Update(1.0,*structure_->Dispnp(),0.0);

  return;
}  // IterUpdateStates()
