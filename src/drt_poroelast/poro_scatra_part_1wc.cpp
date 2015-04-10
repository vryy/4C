/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_1wc.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15251
 </pre>
 *----------------------------------------------------------------------*/


#include "poro_scatra_part_1wc.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "poro_base.H"

#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::DoPoroStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n POROUS MEDIUM SOLVER \n***********************\n";
  }
  //1)  solve the step problem. Methods obtained from poroelast->TimeLoop(sdynparams); --> sdynparams
  //      CUIDADO, aqui vuelve a avanzar el paso de tiempo. Hay que corregir eso.
  //2)  Newton-Raphson iteration
  PoroField()-> Solve();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }
  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  ScaTraField()->Solve();
}

/*----------------------------------------------------------------------*/
//prepare time step                                  rauch/vuong 04/15  |
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::PrepareTimeStep()
{
  IncrementTimeAndStep();

  ScaTraField()->PrepareTimeStep();
  PoroField()->PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::PrepareOutput()
{
  PoroField()-> PrepareOutput();
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::Update()
{
  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  PoroField()-> Update();
  ScaTraField()->Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  ScaTraField()->EvaluateErrorComparedToAnalyticalSol();
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::Output()
{
  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  PoroField()-> Output();
  ScaTraField()->Output();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Part_1WC_PoroToScatra::PORO_SCATRA_Part_1WC_PoroToScatra(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : PORO_SCATRA_Part_1WC(comm, timeparams)
{
  if(comm.MyPID()==0)
    std::cout<<"\n Create PORO_SCATRA_Part_1WC_PoroToScatra algorithm ... \n"<<std::endl;

  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSet> structdofset
    = PoroField()->StructureField()->Discretization()->GetDofSetProxy();

  // check if scatra field has 2 discretizations, so that coupling is possible
  if (ScaTraField()->Discretization()->AddDofSet(structdofset)!=1)
    dserror("unexpected dof sets in scatra field");
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC_PoroToScatra::Timeloop()
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
void POROELAST::PORO_SCATRA_Part_1WC_PoroToScatra::Solve()
{
  DoPoroStep(); // It has its own time and timestep variables, and it increments them by itself.
  SetPoroSolution();
  DoScatraStep(); // It has its own time and timestep variables, and it increments them by itself.
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC_PoroToScatra::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    PoroField()->ReadRestart(restart);
    SetPoroSolution();
    ScaTraField()->ReadRestart(restart);

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
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Part_1WC_ScatraToPoro::PORO_SCATRA_Part_1WC_ScatraToPoro(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : PORO_SCATRA_Part_1WC(comm, timeparams)
{
  if(comm.MyPID()==0)
    std::cout<<"\n Create PORO_SCATRA_Part_1WC_ScatraToPoro algorithm ... \n"<<std::endl;

  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSet> scatradofset
    = ScaTraField()->Discretization()->GetDofSetProxy();

  // check if structure field has 2 discretizations, so that coupling is possible
  if (PoroField()->StructureField()->Discretization()->AddDofSet(scatradofset)!=1)
    dserror("unexpected dof sets in structure field");
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC_ScatraToPoro::Timeloop()
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
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC_ScatraToPoro::Solve()
{
  DoScatraStep(); // It has its own time and timestep variables, and it increments them by itself.
  SetScatraSolution();
  DoPoroStep(); // It has its own time and timestep variables, and it increments them by itself.
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC_ScatraToPoro::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    ScaTraField()->ReadRestart(restart);
    SetScatraSolution();
    PoroField()->ReadRestart(restart);

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
