/*!----------------------------------------------------------------------*/
/*!
\file stru_ale_dyn.cpp
\brief Control routine for structure with ale problems.


<pre>
Maintainer: Markus Gitterle
            gitterle@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15251
</pre>
*/

/*----------------------------------------------------------------------*
 | definitions                                                mgit 04/11 |
 *----------------------------------------------------------------------*/
#ifdef CCADISCRET

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

/*----------------------------------------------------------------------*
 |  headers                                                  mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "stru_ale_dyn.H"
#include "stru_ale_utils.H"
#include "stru_ale_algorithm.H"

#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_condition_utils.H"

#include <Teuchos_TimeMonitor.hpp>

/*----------------------------------------------------------------------*
 | entry point for structure ale in DRT                      mgit 04/11 |
 *----------------------------------------------------------------------*/
void stru_ale_dyn_drt(int disnumsf,int disnumaf,int restart)
{
  // create a communicator
#ifdef PARALLEL
  const Epetra_Comm& comm = DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm();
#else
  Epetra_SerialComm comm;
#endif

  //check if quasistatic analysis
  const Teuchos::ParameterList& sdyn = DRT::Problem::Instance()->StructuralDynamicParams();
  if(sdyn.get<string>("DYNAMICTYP")!= "Statics")
    dserror ("Structure with ale only for quasistatic analysis so in new sti so far.");
  
  // access the structure discretization, make sure it is filled
  Teuchos::RCP<DRT::Discretization> structdis = Teuchos::null;
  structdis = DRT::Problem::Instance()->Dis(disnumsf,0);
  // set degrees of freedom in the discretization
  if (!structdis->Filled() or !structdis->HaveDofs()) structdis->FillComplete();

  // access the ale discretization
  Teuchos::RCP<DRT::Discretization> aledis = Teuchos::null;
  aledis = DRT::Problem::Instance()->Dis(disnumaf,0);
  if (!aledis->Filled()) aledis->FillComplete();

  // we use the structure discretization as layout for the ale discretization
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");
  
  // duplication of structure discretization will follow
  if (aledis->NumGlobalNodes()==0)
  { 
    // fetch the desired material id for the thermo elements
    const int matid = -1;
    
    // create the ale discretization
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<STRU_ALE::UTILS::AleStructureCloneStrategy> > clonewizard
        = Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<STRU_ALE::UTILS::AleStructureCloneStrategy>() );

      clonewizard->CreateMatchingDiscretization(structdis,aledis,matid);
    }
  }
  
  // structure ale object
  Teuchos::RCP<STRU_ALE::Algorithm> stru_ale = Teuchos::rcp(new STRU_ALE::Algorithm(comm));

  // solve the whole problem
  stru_ale->TimeLoop();
  
  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  DRT::Problem::Instance()->AddFieldTest(stru_ale->StructureField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);

  return;
} // stru_ale_dyn_drt()

/*----------------------------------------------------------------------*/
#endif  // CCADISCRET
