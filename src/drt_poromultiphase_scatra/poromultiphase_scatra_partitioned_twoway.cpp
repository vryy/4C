/*----------------------------------------------------------------------*/
/*!
 \file poromultiphase_scatra_partitioned_twoway.cpp

 \brief two-way coupled partitioned algorithm for scalar transport within multiphase porous medium

   \level 3

   \maintainer  Lena Yoshihara
                yoshihara@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
 *----------------------------------------------------------------------*/




#include "poromultiphase_scatra_partitioned_twoway.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

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
  ittol_ = algoparams.sublist("PARTITIONED").get<double>("CONVTOL");

  // initialize increment vectors
  scaincnp_    = Teuchos::rcp(new Epetra_Vector(*(ScatraAlgo()->ScaTraField()->Discretization()->DofRowMap())));
  structincnp_ = Teuchos::rcp(new Epetra_Vector(*(PoroField()->StructDofRowMap())));
  fluidincnp_  = (Teuchos::rcp(new Epetra_Vector(*(PoroField()->FluidDofRowMap()))));

}

/*----------------------------------------------------------------------*
 | setup the monolithic fluid-structure system (called in               |
 | poromultiphase_scatra_dyn.cpp)                      kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::SetupSystem()
{
  PoroField()->SetupSystem();
  return;
}

/*----------------------------------------------------------------------*
 | setup the solver if necessary                        kremheller 03/17 |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::SetupSolver()
{
  PoroField()->SetupSolver();
  return;
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
    std::cout<<"\n";
    std::cout<<"********************************************************************************" <<
        "***************************************************************\n";
    std::cout<<"* PARTITIONED OUTER ITERATION LOOP ----- MULTIPORO  <-------> SCATRA         " <<
        "                                                                 *\n";
    std::cout<<"* STEP: " << std::setw(5) << std::setprecision(4) << std::scientific << Step() << "/"
        << std::setw(5) << std::setprecision(4) << std::scientific << NStep() << ", Time: "
        << std::setw(11) << std::setprecision(4) << std::scientific << Time() << "/"
        << std::setw(11) << std::setprecision(4) << std::scientific << MaxTime() << ", Dt: "
        << std::setw(11) << std::setprecision(4) << std::scientific << Dt() <<
        "                                                                           *"<< std::endl;
  }

  while (stopnonliniter==false)
  {
    itnum++;

    // store scalar from first solution for convergence check (like in
    // elch_algorithm: use current values)
    scaincnp_->Update(1.0,*ScatraAlgo()->ScaTraField()->Phinp(),0.0);
    structincnp_->Update(1.0,*PoroField()->StructDispnp(),0.0);
    fluidincnp_->Update(1.0,*PoroField()->FluidPhinp(),0.0);

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

  // Newton-Raphson iteration
  PoroField()-> TimeStep();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROMULTIPHASESCATRA::PoroMultiPhaseScaTraPartitionedTwoWay::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout<<"\n";
    std::cout<<"*****************************************************************************************************************\n";
    std::cout<<"TRANSPORT SOLVER   \n";
    std::cout<<"*****************************************************************************************************************\n";

  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  ScatraAlgo()->ScaTraField()->Solve();

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
  scaincnp_->Update(1.0,*(ScatraAlgo()->ScaTraField()->Phinp()),-1.0);
  structincnp_->Update(1.0,*(PoroField()->StructDispnp()),-1.0);
  fluidincnp_->Update(1.0,*(PoroField()->FluidPhinp()),-1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  ScatraAlgo()->ScaTraField()->Phinp()->Norm2(&scanorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  PoroField()->StructDispnp()->Norm2(&dispnorm_L2);
  fluidincnp_->Norm2(&fluidincnorm_L2);
  PoroField()->FluidPhinp()->Norm2(&fluidnorm_L2);

  // care for the case that there is (almost) zero scalar
  if (scanorm_L2 < 1e-6) scanorm_L2 = 1.0;
  if (dispnorm_L2 < 1e-6) dispnorm_L2 = 1.0;
  if (fluidnorm_L2 < 1e-6) fluidnorm_L2 = 1.0;

  // print the incremental based convergence check to the screen
  if (Comm().MyPID()==0 )
  {
    std::cout<<"                                                                                                                                              *\n";
    std::cout<<"+------------------------------------------------------------------------------------------------------+                                      *\n";
    std::cout<<"| PARTITIONED OUTER ITERATION STEP ----- MULTIPORO  <-------> SCATRA                                   |                                      *\n";
    printf("+--------------+------------------------+--------------------+--------------------+--------------------+                                      *\n");
    printf("|-  step/max  -|-  tol      [norm]     -|--  scalar-inc    --|--  disp-inc      --|--  fluid-inc     --|                                      *\n");
    printf("|   %3d/%3d    |  %10.3E[L_2 ]      | %10.3E         | %10.3E         | %10.3E         |",
         itnum,itmax_,ittol_,scaincnorm_L2/scanorm_L2,dispincnorm_L2/dispnorm_L2,fluidincnorm_L2/fluidnorm_L2);
    printf("                                      *\n");
    printf("+--------------+------------------------+--------------------+--------------------+--------------------+                                      *\n");
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
      printf("* MULTIPORO  <-------> SCATRA Outer Iteration loop converged after iteration %3d/%3d !                                                        *\n", itnum,itmax_);
      printf("***********************************************************************************************************************************************\n");
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
      printf("* MULTIPORO  <-------> SCATRA Outer Iteration loop not converged in itemax steps                                                              *\n");
      printf("***********************************************************************************************************************************************\n");
      printf("\n");
      printf("\n");
    }
  }

  return stopnonliniter;
}


