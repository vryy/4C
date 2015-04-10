/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_2wc.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/

#include "poro_scatra_part_2wc.H"
#include "poro_base.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Part_2WC::PORO_SCATRA_Part_2WC(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : PORO_SCATRA_Part(comm, timeparams),
    scaincnp_(Teuchos::rcp(new Epetra_Vector(*(ScaTraField()->Phinp())))),
    structincnp_(Teuchos::rcp(new Epetra_Vector(*(PoroField()->StructureField()()->Dispnp())))),
    fluidincnp_(Teuchos::rcp(new Epetra_Vector(*(PoroField()->FluidField()()->Velnp()))))
{
  if(comm.MyPID()==0)
    std::cout<<"\n Create PORO_SCATRA_Part_2WC algorithm ... \n"<<std::endl;

  // the problem is two way coupled, thus each discretization must know the other discretization
  AddDofSets();

  const Teuchos::ParameterList& params = DRT::Problem::Instance()->PoroScatraControlParams();
  // Get the parameters for the ConvergenceCheck
  itmax_ = params.get<int>("ITEMAX"); // default: =10
  ittol_ = params.get<double>("TOLINC_GLOBAL"); // default: =1e-6
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::Timeloop()
{
  //InitialCalculations();

  while (NotFinished())
  {
    PrepareTimeStep();

    Solve();

    PrepareOutput();

    Update();

    Output();

  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    SetScatraSolution();
    SetPoroSolution();

    PoroField()->ReadRestart(restart);
    ScaTraField()->ReadRestart(restart);

    //in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if(PoroField()->HasSubmeshes())
      AddDofSets(true);

    // the variables need to be set on other field
    SetScatraSolution();
    SetPoroSolution();

    //second restart needed due to two way coupling.
    ScaTraField()->ReadRestart(restart);
    PoroField()->ReadRestart(restart);

    //in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if(PoroField()->HasSubmeshes())
      AddDofSets(true);

    SetTimeStep(PoroField()->Time(), restart);

//    // Material pointers to other field were deleted during ReadRestart().
//    // They need to be reset.
//    POROELAST::UTILS::SetMaterialPointersMatchingGrid(PoroField()->StructureField()->Discretization(),
//                                                      PoroField()->FluidField()->Discretization());

    // Material pointers to other field were deleted during ReadRestart().
    // They need to be reset.
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(PoroField()->StructureField()->Discretization(),
                                                      ScaTraField()->Discretization());
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(PoroField()->FluidField()->Discretization(),
                                                      ScaTraField()->Discretization());
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(PoroField()->StructureField()->Discretization(),
                                                      ScaTraField()->Discretization());
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::DoPoroStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n POROUS MEDIUM SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  PoroField()-> Solve();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n  TRANSPORT SOLVER \n***********************\n";
  }

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  ScaTraField()->Solve();

}

/*----------------------------------------------------------------------*/
//prepare time step                                        vuong 08/13  |
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::PrepareTimeStep()
{
  // the global control routine has its own time_ and step_ variables, as well as the single fields
  // keep them in sinc!
  IncrementTimeAndStep();

  SetPoroSolution();
  ScaTraField()->PrepareTimeStep();
  // set structure-based scalar transport values
  SetScatraSolution();

  PoroField()-> PrepareTimeStep();
  SetPoroSolution();
 // SetScatraSolution();
}


/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::PrepareOutput()
{
  PoroField()-> PrepareOutput();
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::Update()
{
  PoroField()-> Update();
  ScaTraField()->Update();

  ScaTraField()->EvaluateErrorComparedToAnalyticalSol();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::Output()
{
  PoroField()-> Output();
  ScaTraField()->Output();
}


/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::Solve()
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
    scaincnp_->Update(1.0,*ScaTraField()->Phinp(),0.0);
    structincnp_->Update(1.0,*PoroField()->StructureField()->Dispnp(),0.0);
    fluidincnp_->Update(1.0,*PoroField()->FluidField()->Velnp(),0.0);

    // set structure-based scalar transport values
    SetScatraSolution();

    // solve structural system
    DoPoroStep();

    // set mesh displacement and velocity fields
    SetPoroSolution();

    // solve scalar transport equation
    DoScatraStep();
    //ScatraEvaluateSolveIterUpdate();

    // check convergence for all fields and stop iteration loop if
    // convergence is achieved overall
    stopnonliniter = ConvergenceCheck(itnum);
  }

  return;
}

/*----------------------------------------------------------------------*
 | convergence check for both fields (scatra & poro) (copied form tsi)
 *----------------------------------------------------------------------*/
bool POROELAST::PORO_SCATRA_Part_2WC::ConvergenceCheck(int itnum)
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
  scaincnp_->Update(1.0,*(ScaTraField()->Phinp()),-1.0);
  structincnp_->Update(1.0,*(PoroField()->StructureField()->Dispnp()),-1.0);
  fluidincnp_->Update(1.0,*(PoroField()->FluidField()->Velnp()),-1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  ScaTraField()->Phinp()->Norm2(&scanorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  PoroField()->StructureField()->Dispnp()->Norm2(&dispnorm_L2);
  fluidincnp_->Norm2(&fluidincnorm_L2);
  PoroField()->FluidField()->Velnp()->Norm2(&fluidnorm_L2);

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


