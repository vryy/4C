/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_2wc.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/

#include "poro_scatra_part_2wc.H"
#include "poro_base.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../linalg/linalg_utils.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Part_2WC::PORO_SCATRA_Part_2WC(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : PORO_SCATRA_Part(comm, timeparams),
    scaincnp_(Teuchos::rcp(new Epetra_Vector(*(scatra_->ScaTraField().Phinp())))),
    structincnp_(Teuchos::rcp(new Epetra_Vector(*(poro_->StructureField()()->Dispnp())))),
    fluidincnp_(Teuchos::rcp(new Epetra_Vector(*(poro_->FluidField()()->Velnp()))))
{
  // the problem is two way coupled, thus each discretization must know the other discretization
  AddDofSets();

  const Teuchos::ParameterList& params = DRT::Problem::Instance()->PoroScatraControlParams();
  // Get the parameters for the ConvergenceCheck
  itmax_ = params.get<int>("ITEMAX"); // default: =10
  ittol_ = params.get<double>("CONVTOL"); // default: =1e-6
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::AddDofSets(bool replace)
{
  // the problem is two way coupled, thus each discretization must know the other discretization
  Teuchos::RCP<DRT::DofSet> structdofset = Teuchos::null;
  Teuchos::RCP<DRT::DofSet> scatradofset = Teuchos::null;

  //get discretizations
  Teuchos::RCP<DRT::Discretization> structdis = poro_->StructureField()->Discretization();
  Teuchos::RCP<DRT::Discretization> scatradis = scatra_->ScaTraField().Discretization();

  if(poro_->HasSubmeshes())
  {
    // build a proxy of the structure discretization for the fluid field (the structure disc. is the bigger one)
    structdofset = structdis->GetDofSetProxy(structdis->NodeColMap(),structdis->ElementColMap());
    // build a proxy of the fluid discretization for the structure field
    scatradofset = scatradis->GetDofSetProxy(scatradis->NodeColMap(),scatradis->ElementColMap());
  }
  else
  {
    // build a proxy of the structure discretization for the fluid field
    structdofset = structdis->GetDofSetProxy();
    // build a proxy of the fluid discretization for the structure field
    scatradofset = scatradis->GetDofSetProxy();
  }

  if(not replace)
  {
    // check if ScatraField has 2 discretizations, so that coupling is possible
    if (scatradis->AddDofSet(structdofset) != 1)
      dserror("unexpected dof sets in fluid field");
    if (structdis->AddDofSet(scatradofset)!=2)
      dserror("unexpected dof sets in structure field");
  }
  else
  {
    scatradis->ReplaceDofSet(1,structdofset);
    structdis->ReplaceDofSet(2,scatradofset);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::Timeloop()
{
  //InitialCalculations();

  while (NotFinished())
  {
    PrepareTimeStep();

    OuterLoop();

    UpdateAndOutput();

  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    SetScatraSolution();
    SetPoroSolution();

    poro_->ReadRestart(restart);
    scatra_->ScaTraField().ReadRestart(restart);

    //in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if(poro_->HasSubmeshes())
      AddDofSets(true);

    // the variables need to be set on other field
    SetScatraSolution();
    SetPoroSolution();

    //second restart needed due to two way coupling.
    scatra_->ScaTraField().ReadRestart(restart);
    poro_->ReadRestart(restart);

    //in case of submeshes, we need to rebuild the subproxies, also (they are reset during restart)
    if(poro_->HasSubmeshes())
      AddDofSets(true);

    SetTimeStep(poro_->Time(), restart);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::DoPoroStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n STRUCTURE SOLVER \n***********************\n";
  }

  // Newton-Raphson iteration
  poro_-> Solve();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
  scatra_->ScaTraField().Solve();

}

/*----------------------------------------------------------------------*/
//prepare time step
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::PrepareTimeStep()
{
  // the global control routine has its own time_ and step_ variables, as well as the single fields
  // keep them in sinc!
  IncrementTimeAndStep();

  //SetPoroSolution();
  scatra_->ScaTraField().PrepareTimeStep();
  // set structure-based scalar transport values
  SetScatraSolution();

  poro_-> PrepareTimeStep();
  SetPoroSolution();
 // SetScatraSolution();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::UpdateAndOutput()
{
  poro_-> PrepareOutput();

  poro_-> Update();
  scatra_->ScaTraField().Update();

  scatra_->ScaTraField().EvaluateErrorComparedToAnalyticalSol();

  poro_-> Output();
  scatra_->ScaTraField().Output();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_2WC::OuterLoop()
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
    scaincnp_->Update(1.0,*scatra_->ScaTraField().Phinp(),0.0);
    structincnp_->Update(1.0,*poro_->StructureField()->Dispnp(),0.0);
    fluidincnp_->Update(1.0,*poro_->FluidField()->Velnp(),0.0);

    // set structure-based scalar transport values
    SetScatraSolution();

   // if(itnum!=1)
   //   poro_->PreparePartitionStep();
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
  scaincnp_->Update(1.0,*(scatra_->ScaTraField().Phinp()),-1.0);
  structincnp_->Update(1.0,*(poro_->StructureField()->Dispnp()),-1.0);
  fluidincnp_->Update(1.0,*(poro_->FluidField()->Velnp()),-1.0);

  // build the L2-norm of the scalar increment and the scalar
  scaincnp_->Norm2(&scaincnorm_L2);
  scatra_->ScaTraField().Phinp()->Norm2(&scanorm_L2);
  structincnp_->Norm2(&dispincnorm_L2);
  poro_->StructureField()->Dispnp()->Norm2(&dispnorm_L2);
  fluidincnp_->Norm2(&fluidincnorm_L2);
  poro_->FluidField()->Velnp()->Norm2(&fluidnorm_L2);

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


