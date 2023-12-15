/*----------------------------------------------------------------------*/
/*! \file

 \brief partitioned one way coupled poroelasticity scalar transport interaction algorithms

\level 2

*/
/*---------------------------------------------------------------------*/


#include "baci_poroelast_scatra_part_1wc.H"

#include "baci_adapter_fld_poro.H"
#include "baci_adapter_scatra_base_algorithm.H"
#include "baci_adapter_str_fpsiwrapper.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_scatra_timint_implicit.H"

BACI_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WC::DoPoroStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n POROUS MEDIUM SOLVER \n***********************\n";
  }
  // 1)  solve the step problem. Methods obtained from poroelast->TimeLoop(sdynparams); -->
  // sdynparams
  //      CUIDADO, aqui vuelve a avanzar el paso de tiempo. Hay que corregir eso.
  // 2)  Newton-Raphson iteration
  PoroField()->Solve();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WC::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }
  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  ScaTraField()->Solve();
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WC::PrepareOutput()
{
  constexpr bool force_prepare = false;
  PoroField()->PrepareOutput(force_prepare);
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WC::Update()
{
  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  PoroField()->Update();
  ScaTraField()->Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  ScaTraField()->EvaluateErrorComparedToAnalyticalSol();
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WC::Output()
{
  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  PoroField()->Output();
  ScaTraField()->CheckAndWriteOutputAndRestart();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELASTSCATRA::PoroScatraPart1WCPoroToScatra::PoroScatraPart1WCPoroToScatra(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraPart1WC(comm, timeparams)
{
  if (comm.MyPID() == 0)
    std::cout << "\n Create PoroScatraPart1WCPoroToScatra algorithm ... \n" << std::endl;
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WCPoroToScatra::Timeloop()
{
  // InitialCalculations();

  while (NotFinished())
  {
    PrepareTimeStep();

    Solve();

    PrepareOutput();

    Update();

    Output();
  }
}

/*----------------------------------------------------------------------*/
// prepare time step                                  rauch/vuong 04/15  |
/*----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WCPoroToScatra::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();
  if (printheader) PrintHeader();

  PoroField()->PrepareTimeStep();
  SetPoroSolution();
  ScaTraField()->PrepareTimeStep();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WCPoroToScatra::Solve()
{
  DoPoroStep();  // It has its own time and timestep variables, and it increments them by itself.
  SetPoroSolution();
  DoScatraStep();  // It has its own time and timestep variables, and it increments them by itself.
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WCPoroToScatra::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    PoroField()->ReadRestart(restart);
    SetPoroSolution();
    ScaTraField()->ReadRestart(restart);

    SetTimeStep(PoroField()->Time(), restart);

    // Material pointers to other field were deleted during ReadRestart().
    // They need to be reset.
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(
        PoroField()->StructureField()->Discretization(), ScaTraField()->Discretization());
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(
        PoroField()->FluidField()->Discretization(), ScaTraField()->Discretization());
  }
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELASTSCATRA::PoroScatraPart1WCScatraToPoro::PoroScatraPart1WCScatraToPoro(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraPart1WC(comm, timeparams)
{
  if (comm.MyPID() == 0)
    std::cout << "\n Create PoroScatraPart1WCScatraToPoro algorithm ... \n" << std::endl;

  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSetInterface> scatradofset =
      ScaTraField()->Discretization()->GetDofSetProxy();

  // check if structure field has 2 discretizations, so that coupling is possible
  if (PoroField()->StructureField()->Discretization()->AddDofSet(scatradofset) != 1)
    dserror("unexpected dof sets in structure field");
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WCScatraToPoro::Timeloop()
{
  // InitialCalculations();

  while (NotFinished())
  {
    PrepareTimeStep();

    Solve();

    PrepareOutput();

    Update();

    Output();
  }
}

/*----------------------------------------------------------------------*/
// prepare time step                                  rauch/vuong 04/15  |
/*----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WCScatraToPoro::PrepareTimeStep(bool printheader)
{
  IncrementTimeAndStep();
  if (printheader) PrintHeader();

  ScaTraField()->PrepareTimeStep();
  SetScatraSolution();
  PoroField()->PrepareTimeStep();
}


/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WCScatraToPoro::Solve()
{
  DoScatraStep();  // It has its own time and timestep variables, and it increments them by itself.
  SetScatraSolution();
  DoPoroStep();  // It has its own time and timestep variables, and it increments them by itself.
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELASTSCATRA::PoroScatraPart1WCScatraToPoro::ReadRestart(int restart)
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  if (restart)
  {
    ScaTraField()->ReadRestart(restart);
    SetScatraSolution();
    PoroField()->ReadRestart(restart);

    SetTimeStep(PoroField()->Time(), restart);

    // Material pointers to other field were deleted during ReadRestart().
    // They need to be reset.
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(
        PoroField()->StructureField()->Discretization(), ScaTraField()->Discretization());
    POROELAST::UTILS::SetMaterialPointersMatchingGrid(
        PoroField()->FluidField()->Discretization(), ScaTraField()->Discretization());
  }
}

BACI_NAMESPACE_CLOSE
