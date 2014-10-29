/*!----------------------------------------------------------------------*/
/*!
\file stru_ale_dyn.cpp
\brief Control routine for structure with ale problems.


<pre>
Maintainer: Philipp Farah
            farah@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15257
</pre>
*/

/*----------------------------------------------------------------------*
 |  headers                                                  mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "wear_dyn.H"
#include "wear_partitioned.H"

#include "../drt_ale/ale_utils_clonestrategy.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"

#include "../drt_inpar/inpar_wear.H"

#include <Teuchos_TimeMonitor.hpp>
#include <Epetra_MpiComm.h>

/*----------------------------------------------------------------------*
 | entry point for structure ale in DRT                      mgit 04/11 |
 *----------------------------------------------------------------------*/
void wear_dyn_drt(int restart)
{
  // create a communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("structure")->Comm();

  // ***********************************************************
  // Setup the problem
  // ***********************************************************
  const Teuchos::ParameterList& sdyn     = DRT::Problem::Instance()->StructuralDynamicParams();
  const Teuchos::ParameterList& wearpara = DRT::Problem::Instance()->WearParams();

  //check if quasistatic analysis
  if(sdyn.get<std::string>("DYNAMICTYP")!= "Statics")
    dserror ("ERROR: Structure with ale only for quasistatic analysis so in new sti so far.");

  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->GetDis("structure");
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs()) structdis->FillComplete();

  // access the ale discretization
  Teuchos::RCP<DRT::Discretization> aledis = Teuchos::null;
  aledis = DRT::Problem::Instance()->GetDis("ale");
  if (!aledis->Filled()) aledis->FillComplete();

  // we use the structure discretization as layout for the ale discretization
  if (structdis->NumGlobalNodes()==0)
    dserror("ERROR: Structure discretization is empty!");

  // clone ale mesh from structure discretization
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(structdis,aledis);
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->Evaluate(params);
  }
  else
    dserror("ERROR: Reading an ALE mesh from the input file is not supported for this problem type.");
  // ***********************************************************

  Teuchos::RCP<WEAR::Algorithm> stru_ale = Teuchos::null;

  // structure ale object
  if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearCoupAlgo>(wearpara,"WEAR_COUPALGO")
      ==  INPAR::CONTACT::wear_stagg)
  {
    stru_ale = Teuchos::rcp(new WEAR::Partitioned(comm));
  }
  else if (DRT::INPUT::IntegralValue<INPAR::CONTACT::WearCoupAlgo>(wearpara,"WEAR_COUPALGO")
      ==  INPAR::CONTACT::wear_iterstagg)
  {
    stru_ale = Teuchos::rcp(new WEAR::Partitioned(comm));
  }
  else
    dserror("ERROR: Chosen algorithm not supported");

  // read restart before joining the time loop
  if (restart!=0)
    stru_ale->ReadRestart(restart);

  // solve the whole problem
  stru_ale->TimeLoop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  DRT::Problem::Instance()->AddFieldTest(stru_ale->StructureField()->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

  return;
} // wear_dyn_drt()

/*----------------------------------------------------------------------*/
