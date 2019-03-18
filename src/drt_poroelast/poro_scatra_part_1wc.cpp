/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_1wc.cpp

 \brief partitioned one way coupled poroelasticity scalar transport interaction algorithms

\level 2

\maintainer Christoph Ager
            ager@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289 15249 </pre>
 *----------------------------------------------------------------------*/


#include "poro_scatra_part_1wc.H"

#include "poro_base.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "../drt_adapter/ad_str_fpsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_scatra/scatra_timint_implicit.H"

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraPart1WC::DoPoroStep()
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
void POROELAST::PoroScatraPart1WC::DoScatraStep()
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
void POROELAST::PoroScatraPart1WC::PrepareOutput() { PoroField()->PrepareOutput(); }

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 04/15  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraPart1WC::Update()
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
void POROELAST::PoroScatraPart1WC::Output()
{
  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  PoroField()->Output();
  ScaTraField()->Output();
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
POROELAST::PoroScatraPart1WCPoroToScatra::PoroScatraPart1WCPoroToScatra(
    const Epetra_Comm& comm, const Teuchos::ParameterList& timeparams)
    : PoroScatraPart1WC(comm, timeparams)
{
  if (comm.MyPID() == 0)
    std::cout << "\n Create PoroScatraPart1WCPoroToScatra algorithm ... \n" << std::endl;
}

/*----------------------------------------------------------------------*
 |                                                   rauch/vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraPart1WCPoroToScatra::Timeloop()
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
void POROELAST::PoroScatraPart1WCPoroToScatra::PrepareTimeStep(bool printheader)
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
void POROELAST::PoroScatraPart1WCPoroToScatra::Solve()
{
  DoPoroStep();  // It has its own time and timestep variables, and it increments them by itself.
  SetPoroSolution();
  DoScatraStep();  // It has its own time and timestep variables, and it increments them by itself.
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraPart1WCPoroToScatra::ReadRestart(int restart)
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
POROELAST::PoroScatraPart1WCScatraToPoro::PoroScatraPart1WCScatraToPoro(
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
void POROELAST::PoroScatraPart1WCScatraToPoro::Timeloop()
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
void POROELAST::PoroScatraPart1WCScatraToPoro::PrepareTimeStep(bool printheader)
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
void POROELAST::PoroScatraPart1WCScatraToPoro::Solve()
{
  DoScatraStep();  // It has its own time and timestep variables, and it increments them by itself.
  SetScatraSolution();
  DoPoroStep();  // It has its own time and timestep variables, and it increments them by itself.
}

/*----------------------------------------------------------------------*
 |                                                         vuong 08/13  |
 *----------------------------------------------------------------------*/
void POROELAST::PoroScatraPart1WCScatraToPoro::ReadRestart(int restart)
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
