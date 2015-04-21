/*----------------------------------------------------------------------*/
/*!
\file sti_dyn.cpp

\brief main control routine for monolithic scalar-thermo interaction

<pre>
Maintainer: Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include <Teuchos_TimeMonitor.hpp>

#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"

#include "../drt_scatra/scatra_resulttest.H"

#include "sti_algorithm.H"
#include "sti_clonestrategy.H"
#include "sti_dyn.H"

/*--------------------------------------------------------------------------------*
 | entry point for simulations of scalar-thermo interaction problems   fang 04/15 |
 *--------------------------------------------------------------------------------*/
void sti_dyn(
    const int& restartstep   //! time step for restart
    )
{
  // access global problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // access communicator
  const Epetra_Comm& comm = problem->GetDis("scatra")->Comm();

  // print logo to screen
  if(comm.MyPID() == 0)
  {
    std::cout <<  "          _ ._  _ , _ ._                                                                                         " << std::endl;
    std::cout <<  "        (_ ' ( `  )_  .__)          ____            _                 _____ _                                    " << std::endl;
    std::cout <<  "      ( (  (    )   `)  ) _)       / ___|  ___ __ _| |_ _ __ __ _    |_   _| |__   ___ _ __ _ __ ___   ___       " << std::endl;
    std::cout <<  "     (__ (_   (_ . _) _) ,__)      \\___ \\ / __/ _` | __| '__/ _` |_____| | | '_ \\ / _ \\ '__| '_ ` _ \\ / _ \\" << std::endl;
    std::cout <<  "         `~~`\\ ' . /`~~`            ___) | (_| (_| | |_| | | (_| |_____| | | | | |  __/ |  | | | | | | (_) |    " << std::endl;
    std::cout <<  "         ,::: ;   ; :::,           |____/ \\___\\__,_|\\__|_|  \\__,_|     |_| |_| |_|\\___|_|  |_| |_| |_|\\___/" << std::endl;
    std::cout <<  "        ':::::::::::::::'           ___       _                      _   _                                       " << std::endl;
    std::cout <<  "             /_ __ \\               |_ _|_ __ | |_ ___ _ __ __ _  ___| |_(_) ___  _ __                           " << std::endl;
    std::cout <<  "     ╔════════════════════╗         | || '_ \\| __/ _ \\ '__/ _` |/ __| __| |/ _ \\| '_ \\                       " << std::endl;
    std::cout <<  "     ║████░░░░░░░░░░░░░░░░╚╗        | || | | | ||  __/ | | (_| | (__| |_| | (_) | | | |                          " << std::endl;
    std::cout <<  "     ║████░░░░░░░░░░░░░░░░░║       |___|_| |_|\\__\\___|_|  \\__,_|\\___|\\__|_|\\___/|_| |_|                    " << std::endl;
    std::cout <<  "     ║████░░░░░░░░░░░░░░░░╔╝                                                                                     " << std::endl;
    std::cout <<  "     ╚════════════════════╝                                                                                      " << std::endl;
  }

  // access, prepare, and check scatra and thermo discretizations
  Teuchos::RCP<DRT::Discretization> scatradis = problem->GetDis("scatra");
  if(!scatradis->Filled() or !scatradis->HaveDofs())
    scatradis->FillComplete();
  if(scatradis->NumGlobalNodes() == 0)
    dserror("The scatra discretization must not be empty, since the thermo discretization needs to be cloned from it!");
  Teuchos::RCP<DRT::Discretization> thermodis = problem->GetDis("thermo");
  if(!thermodis->Filled())
    thermodis->FillComplete();
  if(thermodis->NumGlobalNodes() != 0)
    dserror("The thermo discretization must be empty, since it is cloned from the scatra discretization!");

  // clone thermo discretization from scatra discretization, using clone strategy for scatra-thermo interaction
  DRT::UTILS::CloneDiscretization<STI::ScatraThermoCloneStrategy>(scatradis,thermodis);

  // add proxy of scalar transport degrees of freedom to thermo discretization and vice versa
  if(thermodis->AddDofSet(scatradis->GetDofSetProxy()) != 1)
    dserror("Thermo discretization has illegal number of dofsets!");
  if(scatradis->AddDofSet(thermodis->GetDofSetProxy()) != 1)
    dserror("Scatra discretization has illegal number of dofsets!");

  // add material of scatra elements to thermo elements and vice versa
  for(int i=0; i<scatradis->NumMyColElements(); ++i)
  {
    DRT::Element* scatraele = scatradis->lColElement(i);
    DRT::Element* thermoele = thermodis->gElement(scatraele->Id());

    thermoele->AddMaterial(scatraele->Material());
    scatraele->AddMaterial(thermoele->Material());
  }

  // access scatra parameter list
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // extract and check ID of linear solver
  const int solver_id = scatradyn.get<int>("LINEAR_SOLVER");
  if(solver_id == -1)
    dserror("No linear solver specified in input file section 'SCALAR TRANSPORT DYNAMIC'!");

  // instantiate monolithic algorithm for scatra-thermo interaction
  Teuchos::RCP<STI::Algorithm> sti_algorithm = Teuchos::rcp(new STI::Algorithm(comm,scatradyn,DRT::Problem::Instance()->SolverParams(solver_id)));

  // read restart data if necessary
  if(restartstep)
    sti_algorithm->ReadRestart(restartstep);

  // solve the whole tsi problem
  sti_algorithm->TimeLoop();

  // summarize performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform result tests
  problem->AddFieldTest(Teuchos::rcp<DRT::ResultTest>(new SCATRA::ScaTraResultTest(sti_algorithm->Scatra())));
  problem->AddFieldTest(Teuchos::rcp<DRT::ResultTest>(new SCATRA::ScaTraResultTest(sti_algorithm->Thermo())));
  problem->TestAll(comm);

  return;
}  // sti_dyn()
