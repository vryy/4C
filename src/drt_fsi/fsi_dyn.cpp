
#ifdef CCADISCRET

#include <string>
#include <vector>
#include <set>
#include <functional>

#include <Teuchos_TimeMonitor.hpp>

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"
#include "fsi_dirichletneumannslideale.H"
#include "fsi_robinneumann.H"
#include "fsi_robin.H"
#include "fsi_monolithicfluidsplit.H"
#include "fsi_monolithiclagrange.H"
#include "fsi_monolithicstructuresplit.H"
#include "fsi_monolithicxfem.H"
#include "fsi_partitionedmonolithic.H"
#include "fsi_lungmonolithic.H"
#include "fsi_lungmonolithic_structuresplit.H"
#include "fsi_lungmonolithic_fluidsplit.H"
#include "fsi_constrmonolithic_fluidsplit.H"
#include "fsi_constrmonolithic_structuresplit.H"
#include "fsi_mortarmonolithic_structuresplit.H"
#include "fsi_fluidfluidmonolithic_structuresplit.H"
#include "fsi_mortarmonolithic_fluidsplit.H"
#include "fsi_structureale.H"
#include "fsi_fluid_ale.H"
#include "fsi_fluid_xfem.H"
#include "fsi_utils.H"

#include "fs_monolithic.H"

#include "../drt_fluid/xfluid.H"
#include "../drt_fluid/xfluidfluid.H"

#include "../drt_scatra/scatra_utils.H"

#include "../drt_inpar/inpar_fsi.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_utils_createdis.H"

#ifdef PARALLEL
#include <mpi.h>
#endif

#ifdef PARALLEL
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_condition_selector.H"

#include "../drt_lib/drt_dofset_fixed_size.H"

/*----------------------------------------------------------------------*
 |                                                       m.gee 06/01    |
 | general problem data                                                 |
 | global variable GENPROB genprob is defined in global_control.c       |
 *----------------------------------------------------------------------*/
extern struct _GENPROB     genprob;


/*----------------------------------------------------------------------*/
// entry point for Fluid on Ale in DRT
/*----------------------------------------------------------------------*/
void fluid_ale_drt()
{
#ifdef PARALLEL
  const Epetra_MpiComm& epetrampicomm = dynamic_cast<const Epetra_MpiComm&>(DRT::Problem::Instance()->Dis(genprob.numff,0)->Comm());
  if (!(&epetrampicomm))
    dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  Epetra_MpiComm& comm = const_cast<Epetra_MpiComm&>(epetrampicomm);
#else
  Epetra_SerialComm comm;
#endif

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);
  RCP<DRT::Discretization> aledis   = DRT::Problem::Instance()->Dis(genprob.numaf,0);
  fluiddis->FillComplete();
  aledis->FillComplete();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
  {
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );

      alecreator->CreateMatchingDiscretization(fluiddis,aledis,-1);
    }

    if (comm.MyPID()==0)
    {
      cout << "\n\nCreating ALE discretisation ....\n\n";
    }
  }

  Teuchos::RCP<FSI::FluidAleAlgorithm> fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));
  if (genprob.restart)
  {
    // read the restart information, set vectors and variables
    fluid->ReadRestart(genprob.restart);
  }
  fluid->Timeloop();

  DRT::Problem::Instance()->AddFieldTest(fluid->MBFluidField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);
}


/*----------------------------------------------------------------------*/
// entry point for Fluid on XFEM in DRT
/*----------------------------------------------------------------------*/
void fluid_xfem_drt()
{
#ifdef PARALLEL
  const Epetra_MpiComm& epetrampicomm = dynamic_cast<const Epetra_MpiComm&>(DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm());
  if (!(&epetrampicomm))
    dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  Epetra_MpiComm& comm = const_cast<Epetra_MpiComm&>(epetrampicomm);
#else
  Epetra_SerialComm comm;
#endif


  // make sure the solid dis is filled
  RCP<DRT::Problem> problem = DRT::Problem::Instance();
  problem->Dis(genprob.numsf,0)->FillComplete();


  Teuchos::RCP<FSI::FluidXFEMAlgorithm> xfluid = Teuchos::rcp(new FSI::FluidXFEMAlgorithm(comm));
  if (genprob.restart)
  {
    // read the restart information, set vectors and variables
    xfluid->ReadRestart(genprob.restart);
  }
  xfluid->Timeloop();

  DRT::Problem::Instance()->AddFieldTest(xfluid->MBFluidField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);
}

/*----------------------------------------------------------------------*/
// entry point for Fluid on XFEM in DRT
/*----------------------------------------------------------------------*/
void fluid_xfem2_drt()
{
#ifdef PARALLEL
  const Epetra_MpiComm& epetrampicomm = dynamic_cast<const Epetra_MpiComm&>(DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm());
  if (!(&epetrampicomm))
    dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  Epetra_MpiComm& comm = const_cast<Epetra_MpiComm&>(epetrampicomm);
#else
  Epetra_SerialComm comm;
#endif

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> soliddis = problem->Dis(genprob.numsf,0);
  soliddis->FillComplete();

  RCP<DRT::Discretization> actdis = problem->Dis(genprob.numff,0);
  actdis->FillComplete();

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();

  // now we can reserve dofs for background fluid
  int numglobalnodes = actdis->NumGlobalNodes();
  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  actdis->ReplaceDofSet(maxdofset,true);
  actdis->FillComplete();

  actdis->GetDofSetProxy()->PrintAllDofsets(actdis->Comm());

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
//  Teuchos::RCP<LINALG::Solver> solver =
//    Teuchos::rcp(new LINALG::Solver(problem->FluidSolverParams(),
//                                    actdis->Comm(),
//                                    problem->ErrorFile()->Handle()));
//  //actdis->ComputeNullSpaceIfNecessary(solver->Params());
//
//  FLD::XFluid fluid( actdis,soliddis,*solver,problem->FluidDynamicParams(), xdyn);
//  fluid.Integrate();

   INPAR::XFEM::MovingBoundary moving_boundary = DRT::INPUT::IntegralValue<INPAR::XFEM::MovingBoundary>(xdyn,"XFLUID_BOUNDARY");


   if(moving_boundary == INPAR::XFEM::XFluidStationaryBoundary)
   {
     // no restart required, no moving interface

     // create instance of fluid basis algorithm
     const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();
     Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo = rcp(new ADAPTER::FluidBaseAlgorithm(fdyn,false));

     // run the simulation (timeloop() calls the xfluid-"integrate()" routine)
     fluidalgo->FluidField().TimeLoop();

     // perform result tests if required
     problem->AddFieldTest(fluidalgo->FluidField().CreateFieldTest());
     problem->TestAll(comm);
   }
   else if(moving_boundary == INPAR::XFEM::XFluidMovingBoundary)
   {
     // create instance of fluid xfem algorithm, for moving interfaces
     Teuchos::RCP<FSI::FluidXFEMAlgorithm> fluidalgo = Teuchos::rcp(new FSI::FluidXFEMAlgorithm(comm));
     if (genprob.restart)
     {
       // read the restart information, set vectors and variables
       fluidalgo->ReadRestart(genprob.restart);
     }

     // run the simulation
     fluidalgo->Timeloop();

     // perform result tests if required
     problem->AddFieldTest(fluidalgo->MBFluidField().CreateFieldTest());
     problem->TestAll(comm);
   }
   else if(moving_boundary == INPAR::XFEM::XFSIMovingBoundary) dserror("do not use XFSIMovingBoundary with prb fluid_xfem2");
   else dserror("not a valid XFLUID_BOUNDARY value");


}


/*----------------------------------------------------------------------*/
// entry point for Fluid-Fluid based  on XFEM in DRT
/*----------------------------------------------------------------------*/
void fluid_fluid_ale_drt()
{
  // create a communicator
  #ifdef PARALLEL
   RCP<Epetra_Comm> comm = rcp(DRT::Problem::Instance()->Dis(genprob.numff,0)->Comm().Clone());
  #else
    Epetra_SerialComm comm;
  #endif

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> bgfluiddis = problem->Dis(genprob.numff,0);
  bgfluiddis->FillComplete();

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();


  RCP<DRT::Discretization> embfluiddis = problem->Dis(genprob.numff,1);
  embfluiddis->FillComplete();

  // copy  bgfluid to embfluid
  const int numcolele = bgfluiddis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* ele = bgfluiddis->lColElement(i);
    const DRT::Node*const* elenodes = ele->Nodes();
    RCP<DRT::Element> newele = rcp(ele->Clone());
    embfluiddis->AddElement(newele);
    for (int inode=0; inode < ele->NumNode(); ++inode)
    {
      RCP<DRT::Node> newnode = rcp(elenodes[inode]->Clone());
      embfluiddis->AddNode(newnode);
    }
  }

  embfluiddis->FillComplete();

  //find MovingFluid's elements and nodes
  map<int, DRT::Node*> MovingFluidNodemap;
  map<int, RCP< DRT::Element> > MovingFluidelemap;
  DRT::UTILS::FindConditionObjects(*bgfluiddis, MovingFluidNodemap, MovingFluidelemap, "MovingFluid");

  // local vectors of nodes and elements of moving dis
  vector<int> MovingFluidNodeGIDs;
  vector<int> MovingFluideleGIDs;

  for( map<int, DRT::Node*>::iterator it = MovingFluidNodemap.begin(); it != MovingFluidNodemap.end(); ++it )
    MovingFluidNodeGIDs.push_back( it->first);

  for( map<int, RCP< DRT::Element> >::iterator it = MovingFluidelemap.begin(); it != MovingFluidelemap.end(); ++it )
    MovingFluideleGIDs.push_back( it->first);

  vector<int> NonMovingFluideleGIDs;
  vector<int> NonMovingFluidNodeGIDs;
  for (int iele=0; iele< bgfluiddis->NumMyColElements(); ++iele)
  {
    DRT::Element* bgele = bgfluiddis->lColElement(iele);
    vector<int>::iterator eleiter = find(MovingFluideleGIDs.begin(), MovingFluideleGIDs.end(),bgele->Id() );
    if (eleiter == MovingFluideleGIDs.end())
    {
      NonMovingFluideleGIDs.push_back(bgele->Id());
      int numnode = bgele->NumNode();
      for (int inode=0; inode <numnode; ++inode)
        NonMovingFluidNodeGIDs.push_back(bgele->Nodes()[inode]->Id());
    }
  }

  // copy the conditions to the embedded fluid discretization
  vector<string>          conditions_to_copy;
  conditions_to_copy.push_back("Dirichlet");
  conditions_to_copy.push_back("XFEMCoupling");

  // copy selected conditions to the new discretization
  for (vector<string>::const_iterator conditername = conditions_to_copy.begin();
       conditername != conditions_to_copy.end(); ++conditername)
 {
     vector<DRT::Condition*> conds;
     bgfluiddis->GetCondition(*conditername, conds);
     for (unsigned i=0; i<conds.size(); ++i)
     {
       // We use the same nodal ids and therefore we can just copy the conditions.
       embfluiddis->SetCondition(*conditername, rcp(new DRT::Condition(*conds[i])));
     }
   }

  //information how many processors work at all
  vector<int> allproc(embfluiddis->Comm().NumProc());

  // Gather all informations from all processors
  vector<int> MovingFluideleGIDsall;
  vector<int> MovingFluidNodeGIDsall;

  //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embfluiddis->Comm().NumProc(); ++i) allproc[i] = i;

  //gathers information of MovingFluideleGIDs of all processors
  LINALG::Gather<int>(MovingFluideleGIDs,MovingFluideleGIDsall,(int)embfluiddis->Comm().NumProc(),&allproc[0],embfluiddis->Comm());

  //gathers information of MovingFluidNodeGIDs of all processors
  LINALG::Gather<int>(MovingFluidNodeGIDs,MovingFluidNodeGIDsall,(int)embfluiddis->Comm().NumProc(),&allproc[0],embfluiddis->Comm());

  //information how many processors work at all
  vector<int> allprocbg(bgfluiddis->Comm().NumProc());

  vector<int> NonMovingFluideleGIDsall;
  vector<int> NonMovingFluidNodeGIDsall;

  //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<bgfluiddis->Comm().NumProc(); ++i) allprocbg[i] = i;

  //gathers information of NonMovingFluideleGIDs of all processors
  LINALG::Gather<int>(NonMovingFluideleGIDs,NonMovingFluideleGIDsall,(int)bgfluiddis->Comm().NumProc(),&allprocbg[0],bgfluiddis->Comm());

  //gathers information of NonMovingFluidNodeGIDs of all processors
  LINALG::Gather<int>(NonMovingFluidNodeGIDs,NonMovingFluidNodeGIDsall,(int)bgfluiddis->Comm().NumProc(),&allprocbg[0],bgfluiddis->Comm());

  // ----------------------------------------------------------------------------
  // preparing the background fluid discretization...

  // delete elements and nodes
  for(size_t mv=0; mv<MovingFluideleGIDsall.size(); ++mv)
    bgfluiddis->DeleteElement(MovingFluideleGIDsall.at(mv));

  for(size_t mv=0; mv<MovingFluidNodeGIDsall.size(); ++mv)
    bgfluiddis->DeleteNode(MovingFluidNodeGIDsall.at(mv));

  bgfluiddis->FillComplete();

  // now we can reserve dofs for background fluid
  int numglobalnodes = bgfluiddis->NumGlobalNodes();
  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  bgfluiddis->ReplaceDofSet(maxdofset,true);
  bgfluiddis->FillComplete();

#if defined(PARALLEL)
  vector<int> bgeleids;          // ele ids
  for (int i=0; i<bgfluiddis->NumMyRowElements(); ++i)
  {
    DRT::Element* bgele = bgfluiddis->lRowElement(i);
    int gid = bgele->Id();
    bgeleids.push_back(gid);
  }

  // Background discretization redistribution..
  RCP<Epetra_Map> bgroweles =  rcp(new Epetra_Map(-1,(int)bgeleids.size(),&bgeleids[0],0,bgfluiddis->Comm()));
  RCP<Epetra_Map> bgrownodes = Teuchos::null;
  RCP<Epetra_Map> bgcolnodes = Teuchos::null;

  DRT::UTILS::PartUsingParMetis(bgfluiddis,bgroweles,bgrownodes,bgcolnodes,comm,false);

  RCP<Epetra_Map> bgnewroweles  = Teuchos::null;
  RCP<Epetra_Map> bgnewcoleles  = Teuchos::null;
  bgfluiddis->BuildElementRowColumn(*bgrownodes,*bgcolnodes,bgnewroweles,bgnewcoleles);

  // export nodes and elements to the row map
  bgfluiddis->ExportRowNodes(*bgrownodes);
  bgfluiddis->ExportRowElements(*bgnewroweles);

  // export nodes and elements to the column map (create ghosting)
  bgfluiddis->ExportColumnNodes(*bgcolnodes);
  bgfluiddis->ExportColumnElements(*bgnewcoleles);

  bgfluiddis->FillComplete();
#endif
  //-------------------------------------------------------------------------

  // ----------------------------------------------------------------------------
  // preparing the embedded fluid discretization...

  // delete elements and nodes
  for(size_t nmv=0; nmv<NonMovingFluideleGIDsall.size(); ++nmv)
    embfluiddis->DeleteElement(NonMovingFluideleGIDsall.at(nmv));

  for(size_t nmv=0; nmv<NonMovingFluidNodeGIDsall.size(); ++nmv)
    embfluiddis->DeleteNode(NonMovingFluidNodeGIDsall.at(nmv));

  embfluiddis->FillComplete();

#if defined(PARALLEL)
  vector<int> eleids;          // ele ids
  for (int i=0; i<embfluiddis->NumMyRowElements(); ++i)
  {
    DRT::Element* ele = embfluiddis->lRowElement(i);
    int gid = ele->Id();
    eleids.push_back(gid);
  }

  // Embedded discretization redistribution..
  RCP<Epetra_Map> embroweles =  rcp(new Epetra_Map(-1,(int)eleids.size(),&eleids[0],0,embfluiddis->Comm()));
  RCP<Epetra_Map> embrownodes = Teuchos::null;
  RCP<Epetra_Map> embcolnodes = Teuchos::null;

  DRT::UTILS::PartUsingParMetis(embfluiddis,embroweles,embrownodes,embcolnodes,comm,false);

  RCP<Epetra_Map> embnewroweles  = Teuchos::null;
  RCP<Epetra_Map> embnewcoleles  = Teuchos::null;
  embfluiddis->BuildElementRowColumn(*embrownodes,*embcolnodes,embnewroweles,embnewcoleles);

  // export nodes and elements to the row map
  embfluiddis->ExportRowNodes(*embrownodes);
  embfluiddis->ExportRowElements(*embnewroweles);

  // export nodes and elements to the column map (create ghosting)
  embfluiddis->ExportColumnNodes(*embcolnodes);
  embfluiddis->ExportColumnElements(*embnewcoleles);

  embfluiddis->FillComplete();
#endif
  //------------------------------------------------------------------------------

  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
  {
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );

      alecreator->CreateMatchingDiscretization(embfluiddis,aledis,-1);
    }
    if (comm->MyPID()==0)
    {
      cout << "\n\nCreating ALE discretisation ....\n\n";
    }
  }

  aledis->FillComplete();

  Teuchos::RCP<FSI::FluidAleAlgorithm> alefluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(*comm));
  alefluid->Timeloop();

  problem->AddFieldTest(alefluid->MBFluidField().CreateFieldTest());
  problem->TestAll(*comm);

}

/*----------------------------------------------------------------------*/
// entry point for fluid-fluid-fsi based on XFEM
/*----------------------------------------------------------------------*/
void fluid_fluid_fsi_drt()
{
 // create a communicator
  #ifdef PARALLEL
   RCP<Epetra_Comm> comm = rcp(DRT::Problem::Instance()->Dis(genprob.numff,0)->Comm().Clone());
  #else
    Epetra_SerialComm comm;
  #endif

  /* |--str dofs--|--bgfluid dofs--|--embfluid dofs--|--ale dofs--|-> */

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> structdis = problem->Dis(genprob.numsf,0);
  structdis->FillComplete();

  RCP<DRT::Discretization> bgfluiddis = problem->Dis(genprob.numff,0);
  bgfluiddis->FillComplete();

  // reserve max size of dofs for the background fluid
  const Teuchos::ParameterList xdyn = DRT::Problem::Instance()->XFEMGeneralParams();

  RCP<DRT::Discretization> embfluiddis = problem->Dis(genprob.numff,1);
  embfluiddis->FillComplete();

  // copy  bgfluid to embfluid
  const int numcolele = bgfluiddis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* ele = bgfluiddis->lColElement(i);
    const DRT::Node*const* elenodes = ele->Nodes();
    RCP<DRT::Element> newele = rcp(ele->Clone());
    embfluiddis->AddElement(newele);
    for (int inode=0; inode < ele->NumNode(); ++inode)
    {
      RCP<DRT::Node> newnode = rcp(elenodes[inode]->Clone());
      embfluiddis->AddNode(newnode);
    }
  }

  embfluiddis->FillComplete();

  //find MovingFluid's elements and nodes
  map<int, DRT::Node*> MovingFluidNodemap;
  map<int, RCP< DRT::Element> > MovingFluidelemap;
  DRT::UTILS::FindConditionObjects(*bgfluiddis, MovingFluidNodemap, MovingFluidelemap, "MovingFluid");

  // local vectors of nodes and elements of moving dis
  vector<int> MovingFluidNodeGIDs;
  vector<int> MovingFluideleGIDs;

  for( map<int, DRT::Node*>::iterator it = MovingFluidNodemap.begin(); it != MovingFluidNodemap.end(); ++it )
    MovingFluidNodeGIDs.push_back( it->first);

  for( map<int, RCP< DRT::Element> >::iterator it = MovingFluidelemap.begin(); it != MovingFluidelemap.end(); ++it )
    MovingFluideleGIDs.push_back( it->first);

  // local vectors of nodes and elements of non-moving dis
  vector<int> NonMovingFluideleGIDs;
  vector<int> NonMovingFluidNodeGIDs;
  for (int iele=0; iele< bgfluiddis->NumMyColElements(); ++iele)
  {
    DRT::Element* bgele = bgfluiddis->lColElement(iele);
    vector<int>::iterator eleiter = find(MovingFluideleGIDs.begin(), MovingFluideleGIDs.end(),bgele->Id() );
    if (eleiter == MovingFluideleGIDs.end())
    {
      NonMovingFluideleGIDs.push_back(bgele->Id());
      int numnode = bgele->NumNode();
      for (int inode=0; inode <numnode; ++inode)
        NonMovingFluidNodeGIDs.push_back(bgele->Nodes()[inode]->Id());
    }
  }

  // copy the conditions to the embedded fluid discretization
  vector<string>          conditions_to_copy;
  conditions_to_copy.push_back("Dirichlet");
  conditions_to_copy.push_back("XFEMCoupling");
  conditions_to_copy.push_back("FluidFluidCoupling");

  // copy selected conditions to the new discretization
  for (vector<string>::const_iterator conditername = conditions_to_copy.begin();
       conditername != conditions_to_copy.end(); ++conditername)
  {
     vector<DRT::Condition*> conds;
     bgfluiddis->GetCondition(*conditername, conds);
     for (unsigned i=0; i<conds.size(); ++i)
     {
       // We use the same nodal ids and therefore we can just copy the conditions.
       embfluiddis->SetCondition(*conditername, rcp(new DRT::Condition(*conds[i])));
     }
   }

  //information how many processors work at all
  vector<int> allproc(embfluiddis->Comm().NumProc());

  // Gather all informations from all processors
  vector<int> MovingFluideleGIDsall;
  vector<int> MovingFluidNodeGIDsall;

  //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embfluiddis->Comm().NumProc(); ++i) allproc[i] = i;

  //gathers information of MovingFluideleGIDs of all processors
  LINALG::Gather<int>(MovingFluideleGIDs,MovingFluideleGIDsall,(int)embfluiddis->Comm().NumProc(),&allproc[0],embfluiddis->Comm());

  //gathers information of MovingFluidNodeGIDs of all processors
  LINALG::Gather<int>(MovingFluidNodeGIDs,MovingFluidNodeGIDsall,(int)embfluiddis->Comm().NumProc(),&allproc[0],embfluiddis->Comm());

  //information how many processors work at all
  vector<int> allprocbg(bgfluiddis->Comm().NumProc());

  vector<int> NonMovingFluideleGIDsall;
  vector<int> NonMovingFluidNodeGIDsall;

  //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<bgfluiddis->Comm().NumProc(); ++i) allprocbg[i] = i;

  //gathers information of NonMovingFluideleGIDs of all processors
  LINALG::Gather<int>(NonMovingFluideleGIDs,NonMovingFluideleGIDsall,(int)bgfluiddis->Comm().NumProc(),&allprocbg[0],bgfluiddis->Comm());

  //gathers information of NonMovingFluidNodeGIDs of all processors
  LINALG::Gather<int>(NonMovingFluidNodeGIDs,NonMovingFluidNodeGIDsall,(int)bgfluiddis->Comm().NumProc(),&allprocbg[0],bgfluiddis->Comm());

  // ----------------------------------------------------------------------------
  // preparing the background fluid discretization...

  // delete elements and nodes
  for(size_t mv=0; mv<MovingFluideleGIDsall.size(); ++mv)
    bgfluiddis->DeleteElement(MovingFluideleGIDsall.at(mv));

  for(size_t mv=0; mv<MovingFluidNodeGIDsall.size(); ++mv)
    bgfluiddis->DeleteNode(MovingFluidNodeGIDsall.at(mv));
  bgfluiddis->FillComplete();

  // now we can reserve dofs for background fluid
  int numglobalnodes = bgfluiddis->NumGlobalNodes();
  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  cout << " maxNumMyReservedDofs " << maxNumMyReservedDofs << endl;
  bgfluiddis->ReplaceDofSet(maxdofset,true);
  bgfluiddis->FillComplete();


#if defined(PARALLEL)
  vector<int> bgeleids;          // ele ids
  for (int i=0; i<bgfluiddis->NumMyRowElements(); ++i)
  {
    DRT::Element* bgele = bgfluiddis->lRowElement(i);
    int gid = bgele->Id();
    bgeleids.push_back(gid);
  }

  // Background discretization redistribution..
  RCP<Epetra_Map> bgroweles =  rcp(new Epetra_Map(-1,(int)bgeleids.size(),&bgeleids[0],0,bgfluiddis->Comm()));
  RCP<Epetra_Map> bgrownodes = Teuchos::null;
  RCP<Epetra_Map> bgcolnodes = Teuchos::null;

  DRT::UTILS::PartUsingParMetis(bgfluiddis,bgroweles,bgrownodes,bgcolnodes,comm,false);

  RCP<Epetra_Map> bgnewroweles  = Teuchos::null;
  RCP<Epetra_Map> bgnewcoleles  = Teuchos::null;
  bgfluiddis->BuildElementRowColumn(*bgrownodes,*bgcolnodes,bgnewroweles,bgnewcoleles);

  // export nodes and elements to the row map
  bgfluiddis->ExportRowNodes(*bgrownodes);
  bgfluiddis->ExportRowElements(*bgnewroweles);

  // export nodes and elements to the column map (create ghosting)
  bgfluiddis->ExportColumnNodes(*bgcolnodes);
  bgfluiddis->ExportColumnElements(*bgnewcoleles);

#endif
   bgfluiddis->FillComplete();

  //-------------------------------------------------------------------------

  // ----------------------------------------------------------------------------
  // preparing the embedded fluid discretization...

  // delete elements and nodes
  for(size_t nmv=0; nmv<NonMovingFluideleGIDsall.size(); ++nmv)
    embfluiddis->DeleteElement(NonMovingFluideleGIDsall.at(nmv));

  for(size_t nmv=0; nmv<NonMovingFluidNodeGIDsall.size(); ++nmv)
    embfluiddis->DeleteNode(NonMovingFluidNodeGIDsall.at(nmv));

  // new dofset for embfluiddis which begins after bgfluiddis dofs
  Teuchos::RCP<DRT::DofSet> newdofset = Teuchos::rcp(new DRT::DofSet());
  embfluiddis->ReplaceDofSet(newdofset,true);
  embfluiddis->FillComplete();

#if defined(PARALLEL)
  vector<int> eleids;          // ele ids
  for (int i=0; i<embfluiddis->NumMyRowElements(); ++i)
  {
    DRT::Element* ele = embfluiddis->lRowElement(i);
    int gid = ele->Id();
    eleids.push_back(gid);
  }

  // Embedded discretization redistribution..
  RCP<Epetra_Map> embroweles =  rcp(new Epetra_Map(-1,(int)eleids.size(),&eleids[0],0,embfluiddis->Comm()));
  RCP<Epetra_Map> embrownodes = Teuchos::null;
  RCP<Epetra_Map> embcolnodes = Teuchos::null;

  DRT::UTILS::PartUsingParMetis(embfluiddis,embroweles,embrownodes,embcolnodes,comm,false);

  RCP<Epetra_Map> embnewroweles  = Teuchos::null;
  RCP<Epetra_Map> embnewcoleles  = Teuchos::null;
  embfluiddis->BuildElementRowColumn(*embrownodes,*embcolnodes,embnewroweles,embnewcoleles);

  // export nodes and elements to the row map
  embfluiddis->ExportRowNodes(*embrownodes);
  embfluiddis->ExportRowElements(*embnewroweles);

  // export nodes and elements to the column map (create ghosting)
  embfluiddis->ExportColumnNodes(*embcolnodes);
  embfluiddis->ExportColumnElements(*embnewcoleles);
  embfluiddis->FillComplete();
#endif

  embfluiddis->FillComplete();

  //------------------------------------------------------------------------------

  // create ale elements if the ale discretization is empty
  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
  if (aledis->NumGlobalNodes()==0)
  {
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );

      alecreator->CreateMatchingDiscretization(embfluiddis,aledis,-1);
    }
    if (comm->MyPID()==0)
    {
      cout << "\n\nCreating ALE discretisation ....\n\n";
    }
  }
  aledis->FillComplete();

  // print all dofsets
  structdis->GetDofSetProxy()->PrintAllDofsets(*comm);

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();
  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_fluidfluid_monolithicstructuresplit:
  {
    Teuchos::RCP<FSI::Monolithic> fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicStructureSplit(*comm));

     // now do the coupling setup an create the combined dofmap
    fsi->SetupSystem();

    // here we go...
    fsi->Timeloop();

    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(*comm);
  }
  break;
  default:
  {
    // Any partitioned algorithm
    Teuchos::RCP<FSI::Partitioned> fsi;

    INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");

    if (method==INPAR::FSI::DirichletNeumannSlideale)
    {
        fsi = Teuchos::rcp(new FSI::DirichletNeumannSlideale(*comm));
    }
    else if (method==INPAR::FSI::DirichletNeumann)
    {
      fsi = Teuchos::rcp(new FSI::DirichletNeumann(*comm));
    }
    else if (method==INPAR::FSI::RobinNeumann)
    {
      fsi = Teuchos::rcp(new FSI::RobinNeumann(*comm));
    }
    else
    {
      fsi = Teuchos::rcp(new FSI::Robin(*comm));
    }

//     if (genprob.restart)
//     {
//       // read the restart information, set vectors and variables
//       fsi->ReadRestart(genprob.restart);
//     }
    fsi->Timeloop(fsi);
    DRT::Problem::Instance()->AddFieldTest(fsi->MBFluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(*comm);
  }
  }
}
/*----------------------------------------------------------------------*/
// entry point for (pure) free surface in DRT
/*----------------------------------------------------------------------*/
void fluid_freesurf_drt()
{
#ifdef PARALLEL
  const Epetra_MpiComm& epetrampicomm = dynamic_cast<const Epetra_MpiComm&>(DRT::Problem::Instance()->Dis(genprob.numff,0)->Comm());
  if (!(&epetrampicomm))
    dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  Epetra_MpiComm& comm = const_cast<Epetra_MpiComm&>(epetrampicomm);
#else
  Epetra_SerialComm comm;
#endif


  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  problem->Dis(genprob.numff,0)->FillComplete();
  problem->Dis(genprob.numaf,0)->FillComplete();

  // create ale elements if the ale discretization is empty
  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
  if (aledis->NumGlobalNodes()==0)
  {
    RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);

    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );

    alecreator->CreateMatchingDiscretization(fluiddis,aledis,-1);
  }
    //FSI::UTILS::CreateAleDiscretization();

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
  case fsi_iter_monolithiclagrange:
  {

    INPAR::FSI::LinearBlockSolver linearsolverstrategy = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

    if (linearsolverstrategy==INPAR::FSI::PartitionedAitken or
        linearsolverstrategy==INPAR::FSI::PartitionedVectorExtrapolation or
        linearsolverstrategy==INPAR::FSI::PartitionedJacobianFreeNewtonKrylov or
        linearsolverstrategy==INPAR::FSI::FSIAMG)
      dserror("No partitioned linear solver strategy or FSIAMG supported in Monolithic Free Surface Algorithm. Use PreconditionedKrylov");

    Teuchos::RCP<FSI::MonolithicMainFS> fsi;

    // Monolithic Free Surface Algorithm

    fsi = Teuchos::rcp(new FSI::MonolithicFS(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
    }

    fsi->Timeloop(fsi);

    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
    break;
  }
  default:
  {
    Teuchos::RCP<FSI::FluidAleAlgorithm> fluid;

    // Partitioned FS Algorithm
    fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));

    fluid->Timeloop();

    DRT::Problem::Instance()->AddFieldTest(fluid->MBFluidField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
    break;
  }
  }
}


/*----------------------------------------------------------------------*/
// entry point for FSI using ALE in DRT
/*----------------------------------------------------------------------*/
void fsi_ale_drt()
{
#ifdef PARALLEL
  const Epetra_MpiComm& epetrampicomm = dynamic_cast<const Epetra_MpiComm&>(DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm());
  if (!(&epetrampicomm))
    dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  Epetra_MpiComm& comm = const_cast<Epetra_MpiComm&>(epetrampicomm);
#else
  Epetra_SerialComm comm;
#endif

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       structure dof < fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  problem->Dis(genprob.numsf,0)->FillComplete();
  problem->Dis(genprob.numff,0)->FillComplete();
  problem->Dis(genprob.numaf,0)->FillComplete();

  // create ale elements if the ale discretization is empty
  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
  if (aledis->NumGlobalNodes()==0)
  {
    RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);

    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );

    alecreator->CreateMatchingDiscretization(fluiddis,aledis,-1);
  }
  //FSI::UTILS::CreateAleDiscretization();

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_pseudo_structureale:
  {
    // pseudo FSI problem used to find starting configuration

    Teuchos::RCP<FSI::StructureALE> fsi = Teuchos::rcp(new FSI::StructureALE(comm));

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
    }

    fsi->Timeloop();

    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
    break;
  }
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
  case fsi_iter_monolithiclagrange:
  case fsi_iter_lung_monolithicstructuresplit:
  case fsi_iter_lung_monolithicfluidsplit:
  case fsi_iter_constr_monolithicfluidsplit:
  case fsi_iter_constr_monolithicstructuresplit:
  case fsi_iter_mortar_monolithicstructuresplit:
  case fsi_iter_mortar_monolithicfluidsplit:
  {
    Teuchos::RCP<FSI::MonolithicNOX> fsi;

    INPAR::FSI::LinearBlockSolver linearsolverstrategy = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

    // call constructor to initialise the base class
    if (linearsolverstrategy==INPAR::FSI::PartitionedAitken or
        linearsolverstrategy==INPAR::FSI::PartitionedVectorExtrapolation or
        linearsolverstrategy==INPAR::FSI::PartitionedJacobianFreeNewtonKrylov)
    {
      fsi = Teuchos::rcp(new FSI::PartitionedMonolithic(comm));
    }
    else if (coupling==fsi_iter_monolithicfluidsplit)
    {
      fsi = Teuchos::rcp(new FSI::MonolithicFluidSplit(comm));
    }
    else if (coupling==fsi_iter_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm));
    }
    else if (coupling==fsi_iter_monolithiclagrange)
    {
      fsi = Teuchos::rcp(new FSI::MonolithicLagrange(comm));
    }
    else if (coupling==fsi_iter_lung_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::LungMonolithicStructureSplit(comm));
    }
    else if (coupling==fsi_iter_lung_monolithicfluidsplit)
    {
      fsi = Teuchos::rcp(new FSI::LungMonolithicFluidSplit(comm));
    }
    else if (coupling==fsi_iter_constr_monolithicfluidsplit)
    {
      fsi = Teuchos::rcp(new FSI::ConstrMonolithicFluidSplit(comm));
    }
    else if (coupling==fsi_iter_constr_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::ConstrMonolithicStructureSplit(comm));
    }
    else if (coupling==fsi_iter_mortar_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::MortarMonolithicStructureSplit(comm));
    }
    else if (coupling==fsi_iter_mortar_monolithicfluidsplit)
    {
      fsi = Teuchos::rcp(new FSI::MortarMonolithicFluidSplit(comm));
    }
    else
    {
      dserror("Cannot find appropriate monolithic solver for coupling %d and linear strategy %d",coupling,linearsolverstrategy);
    }

    // read the restart information, set vectors and variables ---
    // be careful, dofmaps might be changed here in a Redistribute call
    if (genprob.restart)
    {
      fsi->ReadRestart(genprob.restart);
    }

    // now do the coupling setup an create the combined dofmap
    fsi->SetupSystem();

    // here we go...
    fsi->Timeloop(fsi);

    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
    break;
  }
  default:
  {
    // Any partitioned algorithm. Stable of working horses.

    Teuchos::RCP<FSI::Partitioned> fsi;

    INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");

    if (method==INPAR::FSI::DirichletNeumannSlideale)
    {
        fsi = Teuchos::rcp(new FSI::DirichletNeumannSlideale(comm));
    }
    else if (method==INPAR::FSI::DirichletNeumann)
    {
      fsi = Teuchos::rcp(new FSI::DirichletNeumann(comm));
    }
    else if (method==INPAR::FSI::RobinNeumann)
    {
      fsi = Teuchos::rcp(new FSI::RobinNeumann(comm));
    }
    else
    {
      fsi = Teuchos::rcp(new FSI::Robin(comm));
    }

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
    }

    fsi->Timeloop(fsi);
    DRT::Problem::Instance()->AddFieldTest(fsi->MBFluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
  }
  }

  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
}




/*----------------------------------------------------------------------*/
// entry point for FSI using XFEM in DRT
/*----------------------------------------------------------------------*/
void xfsi_drt()
{
#ifdef PARALLEL
  const Epetra_MpiComm& epetrampicomm = dynamic_cast<const Epetra_MpiComm&>(DRT::Problem::Instance()->Dis(genprob.numsf,0)->Comm());
  if (!(&epetrampicomm))
    dserror("ERROR: casting Epetra_Comm -> Epetra_MpiComm failed");
  Epetra_MpiComm& comm = const_cast<Epetra_MpiComm&>(epetrampicomm);
#else
  Epetra_SerialComm comm;
#endif

  if (comm.MyPID() == 0)
  {
    cout << endl;
    cout << YELLOW_LIGHT << "       @..@    " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "      (----)      " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "     ( >__< )   " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "     ^^ ~~ ^^  " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "     _     _ _______ _______ _____" << END_COLOR << endl;
    cout << YELLOW_LIGHT << "      \\\\__/  |______ |______   |  " << END_COLOR << endl;
    cout << YELLOW_LIGHT << "     _/  \\\\_ |       ______| __|__" << END_COLOR << endl;
    cout <<  endl << endl;
  }

  RCP<DRT::Problem> problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

#if 0

  // create ale elements if the ale discretization is empty
  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
  if (aledis->NumGlobalNodes()==0)
  {
    RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,1);

    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
      Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );

    alecreator->CreateMatchingDiscretization(fluiddis,aledis,-1);
  }

#endif

  RCP<DRT::Discretization> soliddis = problem->Dis(genprob.numsf,0);
  soliddis->FillComplete();

  RCP<DRT::Discretization> actdis = problem->Dis(genprob.numff,0);
  actdis->FillComplete();

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();

  // now we can reserve dofs for background fluid
  int numglobalnodes = actdis->NumGlobalNodes();
  cout << "numglobalnodes" << numglobalnodes << endl;
  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  actdis->ReplaceDofSet(maxdofset,true);
  actdis->FillComplete();

  actdis->GetDofSetProxy()->PrintAllDofsets(actdis->Comm());

  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_monolithicxfem:
  {
    dserror("fsi_iter_monolithicxfem not implemented yet");

    INPAR::FSI::LinearBlockSolver linearsolverstrategy = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

    if (linearsolverstrategy!=INPAR::FSI::PreconditionedKrylov)
      dserror("Only Newton-Krylov scheme with XFEM fluid");

    Teuchos::RCP<FSI::MonolithicXFEM> fsi;
    fsi = Teuchos::rcp(new FSI::MonolithicXFEM(comm));

    // read the restart information, set vectors and variables ---
    // be careful, dofmaps might be changed here in a Redistribute call
    if (genprob.restart)
    {
      fsi->ReadRestart(genprob.restart);
    }

    // here we go...
    fsi->Timeloop();

    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);

    break;
  }
  case fsi_pseudo_structureale:
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
  case fsi_iter_monolithiclagrange:
    dserror("Unreasonable choice");
  default:
  {
    // Any partitioned algorithm. Stable of working horses.

    Teuchos::RCP<FSI::Partitioned> fsi;

    INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");

    if (method==INPAR::FSI::DirichletNeumann)
    {
      fsi = rcp(new FSI::DirichletNeumann(comm));
    }
    else
      dserror("only Dirichlet-Neumann partitioned schemes with XFEM");

    if (genprob.restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(genprob.restart);
    }

    fsi->Timeloop(fsi);

    DRT::Problem::Instance()->AddFieldTest(fsi->MBFluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
  }
  }

  Teuchos::TimeMonitor::summarize();
}


///*----------------------------------------------------------------------*/
//// entry point for FSI using XFEM in DRT
///*----------------------------------------------------------------------*/
//void xfsi_drt()
//{
//#ifdef PARALLEL
//  Epetra_MpiComm comm(MPI_COMM_WORLD);
//#else
//  Epetra_SerialComm comm;
//#endif
//
//  if (comm.MyPID() == 0)
//  {
//    cout << endl;
//    cout << YELLOW_LIGHT << "       @..@    " << END_COLOR << endl;
//    cout << YELLOW_LIGHT << "      (----)      " << END_COLOR << endl;
//    cout << YELLOW_LIGHT << "     ( >__< )   " << END_COLOR << endl;
//    cout << YELLOW_LIGHT << "     ^^ ~~ ^^  " << END_COLOR << endl;
//    cout << YELLOW_LIGHT << "     _     _ _______ _______ _____" << END_COLOR << endl;
//    cout << YELLOW_LIGHT << "      \\\\__/  |______ |______   |  " << END_COLOR << endl;
//    cout << YELLOW_LIGHT << "     _/  \\\\_ |       ______| __|__" << END_COLOR << endl;
//    cout <<  endl << endl;
//  }
//
//  RCP<DRT::Problem> problem = DRT::Problem::Instance();
//  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();
//
//#if 0
//
//  // create ale elements if the ale discretization is empty
//  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
//  if (aledis->NumGlobalNodes()==0)
//  {
//    RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,1);
//
//    Teuchos::RCP<DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy> > alecreator =
//      Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<FSI::UTILS::AleFluidCloneStrategy>() );
//
//    alecreator->CreateMatchingDiscretization(fluiddis,aledis,-1);
//  }
//
//#endif
//
//  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
//  switch (coupling)
//  {
//  case fsi_iter_monolithicxfem:
//  {
//    dserror("fsi_iter_monolithicxfem not implemented yet");
//
//    INPAR::FSI::LinearBlockSolver linearsolverstrategy = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");
//
//    if (linearsolverstrategy!=INPAR::FSI::PreconditionedKrylov)
//      dserror("Only Newton-Krylov scheme with XFEM fluid");
//
//    Teuchos::RCP<FSI::MonolithicXFEM> fsi;
//    fsi = Teuchos::rcp(new FSI::MonolithicXFEM(comm));
//
//    // read the restart information, set vectors and variables ---
//    // be careful, dofmaps might be changed here in a Redistribute call
//    if (genprob.restart)
//    {
//      fsi->ReadRestart(genprob.restart);
//    }
//
//    // here we go...
//    fsi->Timeloop();
//
//    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
//    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
//    DRT::Problem::Instance()->TestAll(comm);
//
//    break;
//  }
//  case fsi_pseudo_structureale:
//  case fsi_iter_monolithicfluidsplit:
//  case fsi_iter_monolithicstructuresplit:
//  case fsi_iter_monolithiclagrange:
//    dserror("Unreasonable choice");
//  default:
//  {
//    // Any partitioned algorithm. Stable of working horses.
//
//    Teuchos::RCP<FSI::Partitioned> fsi;
//
//    INPAR::FSI::PartitionedCouplingMethod method =
//      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");
//
//    if (method==INPAR::FSI::DirichletNeumann)
//    {
//      fsi = rcp(new FSI::DirichletNeumann(comm));
//    }
//    else
//      dserror("only Dirichlet-Neumann partitioned schemes with XFEM");
//
//    if (genprob.restart)
//    {
//      // read the restart information, set vectors and variables
//      fsi->ReadRestart(genprob.restart);
//    }
//
//    fsi->Timeloop(fsi);
//
//    DRT::Problem::Instance()->AddFieldTest(fsi->MBFluidField().CreateFieldTest());
//    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
//    DRT::Problem::Instance()->TestAll(comm);
//  }
//  }
//
//  Teuchos::TimeMonitor::summarize();
//}


#endif
