/*----------------------------------------------------------------------*/
/*!
 \file poro_scatra.cpp

 \brief  scalar transport in porous media

 <pre>
   Maintainer: Anh-Tu Vuong
               vuong@lnm.mw.tum.de
               http://www.lnm.mw.tum.de
               089 - 289-15264
 </pre>
 *-----------------------------------------------------------------------*/

/*
 *  Implementation of passive scalar transport in porous media.
 *  Done by Miguel Urrecha (miguel.urrecha@upm.es)
*/


#include "poro_scatra.H"

#include "../drt_scatra/passive_scatra_algorithm.H"
#include "../drt_inpar/inpar_scatra.H"
#include "poroelast_utils.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_scatra/scatra_utils_clonestrategy.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/ad_str_fsiwrapper.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_RefCountPtr.hpp>

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
PORO_SCATRA::PartPORO_SCATRA::PartPORO_SCATRA(const Epetra_Comm& comm,
    const Teuchos::ParameterList& timeparams)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  const Teuchos::ParameterList& scatradyn  = problem->ScalarTransportDynamicParams();

  // get the solver number used for ScalarTransport solver
  const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");

  //1.- Set time and step.
  dt_ = timeparams.get<double> ("TIMESTEP");
  numstep_ = timeparams.get<int> ("NUMSTEP");
  timemax_ = timeparams.get<double> ("MAXTIME");

  step_ = 0;
  time_ = 0.0;

  //2.- Setup discretizations.
  SetupDiscretizations(comm);

  //3.- Create the two uncoupled subproblems.
  poroelast_subproblem_ = Teuchos::rcp(new POROELAST::Monolithic(comm, timeparams));
  scatra_subproblem_ = Teuchos::rcp(new ADAPTER::ScaTraBaseAlgorithm(timeparams,true,"scatra",problem->SolverParams(linsolvernumber)));
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::Timeloop()
{

  // 1.- Output of initial state
  scatra_subproblem_->ScaTraField().Output();

  // 2.- Actual time loop
  while (NotFinished()) // Actual time loop
  {
    IncrementTimeAndStep(); // This is just for control, not needed at all (not the time variables that the "Do"functions take).
    DoPoroStep(); // It has its own time and timestep variables, and it increments them by itself.
    SetPoroSolution();
    DoScatraStep(); // It has its own time and timestep variables, and it increments them by itself.
  }

}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::DoPoroStep()
{
  // solve the step problem. Methods obtained from poroelast->TimeLoop(sdynparams); --> sdynparams
  poroelast_subproblem_-> PrepareTimeStep(); //			CUIDADO, aqui vuelve a avanzar el paso de tiempo. Hay que corregir eso.
  // Newton-Raphson iteration
  poroelast_subproblem_-> NewtonFull();
  // calculate stresses, strains, energies
  poroelast_subproblem_-> PrepareOutput();
  // update all single field solvers
  poroelast_subproblem_-> Update();
  // write output to screen and files
  poroelast_subproblem_-> Output();
}

void PORO_SCATRA::PartPORO_SCATRA::DoScatraStep()
{
  const Epetra_Comm& comm =
      DRT::Problem::Instance()->GetDis("structure")->Comm();

  if (comm.MyPID() == 0)
  {
    cout
        << "\n***********************\n TRANSPORT SOLVER \n***********************\n";
  }
  // -------------------------------------------------------------------
  // prepare time step
  // -------------------------------------------------------------------
  scatra_subproblem_->ScaTraField().PrepareTimeStep();

  // -------------------------------------------------------------------
  //                  solve nonlinear / linear equation
  // -------------------------------------------------------------------
  scatra_subproblem_->ScaTraField().Solve();

  // -------------------------------------------------------------------
  //                         update solution
  //        current solution becomes old solution of next timestep
  // -------------------------------------------------------------------
  scatra_subproblem_->ScaTraField().Update();

  // -------------------------------------------------------------------
  // evaluate error for problems with analytical solution
  // -------------------------------------------------------------------
  scatra_subproblem_->ScaTraField().EvaluateErrorComparedToAnalyticalSol();

  // -------------------------------------------------------------------
  //                         output of solution
  // -------------------------------------------------------------------
  scatra_subproblem_->ScaTraField().Output();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::IncrementTimeAndStep()
{
  step_ += 1;
  time_ += dt_;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::SetPoroSolution()
{
  SetMeshDisp();
  SetVelocityFields();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::SetVelocityFields()
{
  scatra_subproblem_->ScaTraField().SetVelocityField(
      poroelast_subproblem_->FluidField().ConvectiveVel(), //convective vel.
      Teuchos::null, //acceleration
      poroelast_subproblem_->FluidField().Velnp(), //velocity
      Teuchos::null, //fsvel
      Teuchos::null, //dofset
      poroelast_subproblem_->FluidField().Discretization()); //discretization
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::SetMeshDisp()
{
  scatra_subproblem_->ScaTraField().ApplyMeshMovement(
      poroelast_subproblem_->FluidField().Dispnp(),
      poroelast_subproblem_->FluidField().Discretization());
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::SetupDiscretizations(const Epetra_Comm& comm)
{
  // Scheme    : the structure discretization is received from the input. Then, an ale-fluid disc.is cloned from the struct. one.
  //  After that, an ale-scatra disc. is cloned from the fluid disc. already created.

  DRT::Problem* problem = DRT::Problem::Instance();

  //1.-Initialization.
  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");

  //1.2.-Set degrees of freedom in the str. discretization
  if (!structdis->Filled() or !structdis->HaveDofs())
    structdis->FillComplete();

  //2.- Access the fluid discretization, make sure it's empty, and fill it cloning the structural one.
  if (!fluiddis->Filled())
    fluiddis->FillComplete();

  if (structdis->NumGlobalNodes() == 0)
    dserror("Structure discretization is empty!");

  if (fluiddis->NumGlobalNodes()==0)
  {
    // create the fluid discretization
    DRT::UTILS::CloneDiscretization<POROELAST::UTILS::PoroelastCloneStrategy>(structdis,fluiddis);
  }
  else
  dserror("Structure AND Fluid discretization present. This is not supported.");

  //3.-Access the scatra discretization, make sure it's empty, and fill it by cloning the structural one.
  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

  if(!scatradis->Filled())
  scatradis->FillComplete();

  if (scatradis->NumGlobalNodes()==0)
  {
    // create the fluid scatra discretization
    DRT::UTILS::CloneDiscretization<SCATRA::ScatraFluidCloneStrategy>(structdis,fluiddis);
  }
  else
  dserror("Structure AND ScaTra discretization present. This is not supported.");
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::ReadRestart()
{
  // read restart information, set vectors and variables
  // (Note that dofmaps might have changed in a redistribution call!)
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    poroelast_subproblem_->ReadRestart(restart);
    scatra_subproblem_->ScaTraField().ReadRestart(restart);

    time_ = poroelast_subproblem_->FluidField().Time();
    step_ = poroelast_subproblem_->FluidField().Step();
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::SetupSystem()
{
  poroelast_subproblem_->SetupSystem();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PORO_SCATRA::PartPORO_SCATRA::TestResults(const Epetra_Comm& comm)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  problem->AddFieldTest(poroelast_subproblem_->StructureField()->CreateFieldTest());
  problem->AddFieldTest(poroelast_subproblem_->FluidField().CreateFieldTest());
  problem->AddFieldTest(scatra_subproblem_->CreateScaTraFieldTest());
  problem->TestAll(comm);
}

