/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra_part_1wc.cpp

 \brief

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *----------------------------------------------------------------------*/


#include "poro_scatra_part_1wc.H"

#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "poro_base.H"

#include "../drt_adapter/ad_str_fsiwrapper.H"
#include "../drt_adapter/ad_fld_poro.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
  //3)  calculate stresses, strains, energies
  //4)  update all single field solvers
  //5)  write output to screen and files
  PoroField()-> DoTimeStep();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::DoScatraStep()
{
  if (Comm().MyPID() == 0)
  {
    std::cout
        << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }
  // -------------------------------------------------------------------
  // prepare time step
  // -------------------------------------------------------------------
  ScaTraField()->PrepareTimeStep();

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  ScaTraField()->Solve();

  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  ScaTraField()->Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  ScaTraField()->EvaluateErrorComparedToAnalyticalSol();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  ScaTraField()->Output();
}

/*----------------------------------------------------------------------*/
//prepare time step
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC::PrepareTimeStep()
{
  IncrementTimeAndStep();
  //PrintHeader();

  //SetScatraSolution();
  //SetPoroSolution();

  //PrepareTimeStep of single fields is called in DoPoroStep and DoScatraStep
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Part_1WC_PoroToScatra::PORO_SCATRA_Part_1WC_PoroToScatra(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : PORO_SCATRA_Part_1WC(comm, timeparams)
{
  // build a proxy of the structure discretization for the scatra field
  Teuchos::RCP<DRT::DofSet> structdofset
    = PoroField()->StructureField()->Discretization()->GetDofSetProxy();

  // check if scatra field has 2 discretizations, so that coupling is possible
  if (ScaTraField()->Discretization()->AddDofSet(structdofset)!=1)
    dserror("unexpected dof sets in scatra field");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC_PoroToScatra::Timeloop()
{
  //InitialCalculations();

  while (NotFinished())
  {
    PrepareTimeStep();

    DoPoroStep(); // It has its own time and timestep variables, and it increments them by itself.
    SetPoroSolution();
    DoScatraStep(); // It has its own time and timestep variables, and it increments them by itself.
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
POROELAST::PORO_SCATRA_Part_1WC_ScatraToPoro::PORO_SCATRA_Part_1WC_ScatraToPoro(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
  : PORO_SCATRA_Part_1WC(comm, timeparams)
{
  // build a proxy of the scatra discretization for the structure field
  Teuchos::RCP<DRT::DofSet> scatradofset
    = ScaTraField()->Discretization()->GetDofSetProxy();

  // check if structure field has 2 discretizations, so that coupling is possible
  if (PoroField()->StructureField()->Discretization()->AddDofSet(scatradofset)!=1)
    dserror("unexpected dof sets in structure field");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void POROELAST::PORO_SCATRA_Part_1WC_ScatraToPoro::Timeloop()
{
  //InitialCalculations();

  while (NotFinished())
  {
    PrepareTimeStep();

    DoScatraStep(); // It has its own time and timestep variables, and it increments them by itself.
    SetScatraSolution();
    DoPoroStep(); // It has its own time and timestep variables, and it increments them by itself.
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
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
  }
}
