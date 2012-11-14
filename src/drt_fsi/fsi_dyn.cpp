


#include <string>
#include <vector>
#include <set>
#include <functional>

#include <Teuchos_TimeMonitor.hpp>

#include "fsi_dyn.H"
#include "fsi_dirichletneumann.H"
#include "fsi_dirichletneumannslideale.H"
#include "fsi_monolithicfluidsplit.H"
#include "fsi_monolithicstructuresplit.H"
#include "fsi_lungmonolithic.H"
#include "fsi_lungmonolithic_structuresplit.H"
#include "fsi_lungmonolithic_fluidsplit.H"
#include "fsi_constrmonolithic_fluidsplit.H"
#include "fsi_constrmonolithic_structuresplit.H"
#include "fsi_mortarmonolithic_structuresplit.H"
#include "fsi_fluidfluidmonolithic_structuresplit_nonox.H"
#include "fsi_mortarmonolithic_fluidsplit.H"
#include "fsi_structureale.H"
#include "fsi_fluid_ale.H"
#include "fsi_fluid_xfem.H"
#include "fsi_utils.H"
#include "fsi_resulttest.H"

#include "fs_monolithic.H"

#include "../drt_fluid/xfluid.H"
#include "../drt_fluid/xfluidfluid.H"

#include "../drt_scatra/scatra_utils.H"

#include "../drt_ale/ale_utils_clonestrategy.H"

#include "../drt_inpar/inpar_fsi.H"
#include "../drt_lib/drt_resulttest.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../linalg/linalg_utils.H"

#include "../drt_lib/drt_globalproblem.H"

#include "../drt_lib/drt_condition_utils.H"
#include "../drt_lib/drt_condition_selector.H"

#include "../drt_lib/drt_dofset_fixed_size.H"

#include "../drt_adapter/ad_str_structure.H"
#include "../drt_adapter/adapter_coupling.H"
#include "../drt_adapter/adapter_coupling_mortar.H"
#include "../drt_adapter/ad_fld_fluid.H"
#include "../drt_adapter/ad_fld_moving_boundary.H"

#include "../drt_lib/drt_discret_xfem.H"


/*----------------------------------------------------------------------*/
// entry point for Fluid on Ale in DRT
/*----------------------------------------------------------------------*/
void fluid_ale_drt()
{
  DRT::Problem* problem = DRT::Problem::Instance();

  const Epetra_Comm& comm = problem->GetDis("fluid")->Comm();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  RCP<DRT::Discretization> aledis   = problem->GetDis("ale");
  fluiddis->FillComplete();
  aledis->FillComplete();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
  }
  else  // filled ale discretization
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis,aledis))
      dserror("Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
  }

  Teuchos::RCP<FSI::FluidAleAlgorithm> fluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(comm));
  const int restart = problem->Restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    fluid->ReadRestart(restart);
  }
  fluid->Timeloop();

  DRT::Problem::Instance()->AddFieldTest(fluid->MBFluidField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(comm);
}


/*----------------------------------------------------------------------*/
// entry point for Fluid on XFEM in DRT
/*----------------------------------------------------------------------*/
void fluid_xfem2_drt()
{
#ifdef PARALLEL
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("structure")->Comm();
#else
  Epetra_SerialComm comm;
#endif

  DRT::Problem* problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> soliddis = problem->GetDis("structure");
  soliddis->FillComplete();

  RCP<DRT::Discretization> actdis = problem->GetDis("fluid");
  actdis->FillComplete();


  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();

  // now we can reserve dofs for background fluid
  int numglobalnodes = actdis->NumGlobalNodes();
  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  actdis->ReplaceDofSet(maxdofset,true);
  actdis->FillComplete();

  actdis->GetDofSetProxy()->PrintAllDofsets(actdis->Comm());



  const Teuchos::ParameterList& xfdyn     = DRT::Problem::Instance()->XFluidDynamicParams();

  INPAR::XFEM::MovingBoundary moving_boundary = DRT::INPUT::IntegralValue<INPAR::XFEM::MovingBoundary>(xfdyn.sublist("GENERAL"),"XFLUID_BOUNDARY");

  if(moving_boundary == INPAR::XFEM::XFluidStationaryBoundary)
  {
    // no restart required, no moving interface

    // create instance of fluid basis algorithm
    const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn,false));

    // run the simulation, calls the xfluid-"integrate()" routine
    fluidalgo->FluidField().Integrate();

    // perform result tests if required
    problem->AddFieldTest(fluidalgo->FluidField().CreateFieldTest());
    problem->TestAll(comm);
  }
  else if(moving_boundary == INPAR::XFEM::XFluidMovingBoundary)
  {
    // create instance of fluid xfem algorithm, for moving interfaces
    Teuchos::RCP<FSI::FluidXFEMAlgorithm> fluidalgo = Teuchos::rcp(new FSI::FluidXFEMAlgorithm(comm));

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fluidalgo->ReadRestart(restart);
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
  DRT::Problem* problem = DRT::Problem::Instance();

  // create a communicator
  #ifdef PARALLEL
   RCP<Epetra_Comm> comm = Teuchos::rcp(problem->GetDis("fluid")->Comm().Clone());
  #else
    Epetra_SerialComm comm;
  #endif

  RCP<DRT::Discretization> bgfluiddis = problem->GetDis("fluid");
  bgfluiddis->FillComplete();

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();


  RCP<DRT::Discretization> embfluiddis = problem->GetDis("xfluid");
  embfluiddis->FillComplete();

  // -------------------------------------------------------------------
  // --------------- copy  bgfluid to embfluid
  const int numcolele = bgfluiddis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* ele = bgfluiddis->lColElement(i);
    const DRT::Node*const* elenodes = ele->Nodes();
    RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
    embfluiddis->AddElement(newele);
    for (int inode=0; inode < ele->NumNode(); ++inode)
    {
      RCP<DRT::Node> newnode = Teuchos::rcp(elenodes[inode]->Clone());
      embfluiddis->AddNode(newnode);
    }
  }

  embfluiddis->FillComplete();

  // -------------------------------------------------------------------
  // ---------------- find MovingFluid's elements and nodes
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

  // ----------------------------------------------------------------
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
       embfluiddis->SetCondition(*conditername, Teuchos::rcp(new DRT::Condition(*conds[i])));
     }
   }

  // --------------------------------------------------------------------------
  // ------------------ gather information for moving fluid  -------------------

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

  // -------------------------------------------------------------------------------
  // -------------- now build the nonmoving vectors from the gathered moving vectors
  vector<int> NonMovingFluideleGIDs;
  vector<int> NonMovingFluidNodeGIDs;
  for (int iele=0; iele< bgfluiddis->NumMyColElements(); ++iele)
  {
    DRT::Element* bgele = bgfluiddis->lColElement(iele);
    vector<int>::iterator eleiter = find(MovingFluideleGIDsall.begin(), MovingFluideleGIDsall.end(),bgele->Id() );
    if (eleiter == MovingFluideleGIDsall.end())
    {
      NonMovingFluideleGIDs.push_back(bgele->Id());
      int numnode = bgele->NumNode();
      for (int inode=0; inode <numnode; ++inode)
        NonMovingFluidNodeGIDs.push_back(bgele->Nodes()[inode]->Id());
    }
  }

  // --------------------------------------------------------------------------
  // ------------------ gather information for non moving fluid ---------------
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
  RCP<Epetra_Map> bgroweles =  Teuchos::rcp(new Epetra_Map(-1,(int)bgeleids.size(),&bgeleids[0],0,bgfluiddis->Comm()));
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

  embfluiddis->CheckFilledGlobally();

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
  RCP<Epetra_Map> embroweles =  Teuchos::rcp(new Epetra_Map(-1,(int)eleids.size(),&eleids[0],0,embfluiddis->Comm()));
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

  RCP<DRT::Discretization> aledis = problem->GetDis("ale");

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(embfluiddis,aledis);
  }
  else  // ale discretization in input file
    dserror("Providing an ALE mesh is not supported for this problemtype.");

  aledis->FillComplete();

  Teuchos::RCP<FSI::FluidAleAlgorithm> alefluid = Teuchos::rcp(new FSI::FluidAleAlgorithm(*comm));
  const int restart = DRT::Problem::Instance()->Restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    alefluid->ReadRestart(restart);
  }
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
   RCP<Epetra_Comm> comm = Teuchos::rcp(DRT::Problem::Instance()->GetDis("fluid")->Comm().Clone());
  #else
    Epetra_SerialComm comm;
  #endif

  /* |--str dofs--|--bgfluid dofs--|--embfluid dofs--|--ale dofs--|-> */

  DRT::Problem* problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  structdis->FillComplete();

  RCP<DRT::Discretization> bgfluiddis = problem->GetDis("fluid");
  bgfluiddis->FillComplete();

  // reserve max size of dofs for the background fluid
  const Teuchos::ParameterList xdyn = DRT::Problem::Instance()->XFEMGeneralParams();

  RCP<DRT::Discretization> embfluiddis = problem->GetDis("xfluid");
  embfluiddis->FillComplete();

  // -------------------------------------------------------------------
  // --------------- copy  bgfluid to embfluid
  const int numcolele = bgfluiddis->NumMyColElements();
  for (int i=0; i<numcolele; ++i)
  {
    DRT::Element* ele = bgfluiddis->lColElement(i);
    const DRT::Node*const* elenodes = ele->Nodes();
    RCP<DRT::Element> newele = Teuchos::rcp(ele->Clone());
    embfluiddis->AddElement(newele);
    for (int inode=0; inode < ele->NumNode(); ++inode)
    {
      RCP<DRT::Node> newnode = Teuchos::rcp(elenodes[inode]->Clone());
      embfluiddis->AddNode(newnode);
    }
  }

  embfluiddis->FillComplete();

  // -------------------------------------------------------------------
  // ---------------- find MovingFluid's elements and nodes
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


  // ----------------------------------------------------------------
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
       embfluiddis->SetCondition(*conditername, Teuchos::rcp(new DRT::Condition(*conds[i])));
     }
   }

  // --------------------------------------------------------------------------
  // ------------------ gather information for moving fluid -------------------
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

  // -------------------------------------------------------------------------------
  // -------------- now build the nonmoving vectors from the gathered moving vectors
  vector<int> NonMovingFluideleGIDs;
  vector<int> NonMovingFluidNodeGIDs;
  for (int iele=0; iele< bgfluiddis->NumMyColElements(); ++iele)
  {
    DRT::Element* bgele = bgfluiddis->lColElement(iele);
    vector<int>::iterator eleiter = find(MovingFluideleGIDsall.begin(), MovingFluideleGIDsall.end(),bgele->Id() );
    if (eleiter == MovingFluideleGIDsall.end())
    {
      NonMovingFluideleGIDs.push_back(bgele->Id());
      int numnode = bgele->NumNode();
      for (int inode=0; inode <numnode; ++inode)
        NonMovingFluidNodeGIDs.push_back(bgele->Nodes()[inode]->Id());
    }
  }

  // --------------------------------------------------------------------------
  // ------------------ gather information for non moving fluid ---------------
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
  RCP<Epetra_Map> bgroweles =  Teuchos::rcp(new Epetra_Map(-1,(int)bgeleids.size(),&bgeleids[0],0,bgfluiddis->Comm()));
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

  embfluiddis->CheckFilledGlobally();

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
  RCP<Epetra_Map> embroweles =  Teuchos::rcp(new Epetra_Map(-1,(int)eleids.size(),&eleids[0],0,embfluiddis->Comm()));
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
  RCP<DRT::Discretization> aledis = problem->GetDis("ale");
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(embfluiddis,aledis);
  }
  else  // ale discretization in input file
    dserror("Providing an ALE mesh is not supported for problemtype Fluid_Fluid_FSI.");

  aledis->FillComplete();

  // print all dofsets
  structdis->GetDofSetProxy()->PrintAllDofsets(*comm);

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();
  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_fluidfluid_monolithicstructuresplit:
  {
    Teuchos::RCP<FSI::MonolithicNoNOX> fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicStructureSplitNoNOX(*comm,fsidyn));

     // now do the coupling setup an create the combined dofmap
    fsi->SetupSystem();

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(restart);
    }

    // here we go...
    fsi->Timeloop();

    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());
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
        fsi = Teuchos::rcp(new FSI::DirichletNeumannSlideale(*comm));
    else if (method==INPAR::FSI::DirichletNeumann)
      fsi = Teuchos::rcp(new FSI::DirichletNeumann(*comm));
    else
      dserror("unsupported partitioned FSI scheme");

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(restart);
    }


    fsi->Timeloop(fsi);
    DRT::Problem::Instance()->AddFieldTest(fsi->MBFluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());
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
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();
#else
  Epetra_SerialComm comm;
#endif


  DRT::Problem* problem = DRT::Problem::Instance();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  problem->GetDis("fluid")->FillComplete();
  problem->GetDis("ale")->FillComplete();

  // get discretizations
  RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  RCP<DRT::Discretization> aledis = problem->GetDis("ale");

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
  }
  else  // filled ale discretization
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis,aledis))
      dserror("Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
  }

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
  {

    INPAR::FSI::LinearBlockSolver linearsolverstrategy = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

    if (linearsolverstrategy==INPAR::FSI::FSIAMG)
      dserror("No FSIAMG supported in Monolithic Free Surface Algorithm. Use PreconditionedKrylov");

    Teuchos::RCP<FSI::MonolithicMainFS> fsi;

    // Monolithic Free Surface Algorithm

    fsi = Teuchos::rcp(new FSI::MonolithicFS(comm,fsidyn));

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(restart);
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
  DRT::Problem* problem = DRT::Problem::Instance();

  #ifdef PARALLEL
  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();
#else
  Epetra_SerialComm comm;
#endif

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       structure dof < fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  problem->GetDis("structure")->FillComplete();
  problem->GetDis("fluid")->FillComplete();
  problem->GetDis("ale")->FillComplete();

  // get discretizations
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0) // empty ale discretization
  {
    DRT::UTILS::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(fluiddis,aledis);
  }
  else  // filled ale discretization
  {
    if (!FSI::UTILS::FluidAleNodesDisjoint(fluiddis,aledis))
      dserror("Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
  }

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_pseudo_structureale:
  {
    // pseudo FSI problem used to find starting configuration

    Teuchos::RCP<FSI::StructureALE> fsi = Teuchos::rcp(new FSI::StructureALE(comm));

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(restart);
    }

    fsi->Timeloop();

    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
    break;
  }
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
  case fsi_iter_lung_monolithicstructuresplit:
  case fsi_iter_lung_monolithicfluidsplit:
  case fsi_iter_constr_monolithicfluidsplit:
  case fsi_iter_constr_monolithicstructuresplit:
  case fsi_iter_mortar_monolithicstructuresplit:
  case fsi_iter_mortar_monolithicfluidsplit:
  {
    Teuchos::RCP<FSI::Monolithic> fsi;

    INPAR::FSI::LinearBlockSolver linearsolverstrategy = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsidyn,"LINEARBLOCKSOLVER");

    // call constructor to initialise the base class
    if (coupling==fsi_iter_monolithicfluidsplit)
    {
      fsi = Teuchos::rcp(new FSI::MonolithicFluidSplit(comm,fsidyn));
    }
    else if (coupling==fsi_iter_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::MonolithicStructureSplit(comm,fsidyn));
    }
    else if (coupling==fsi_iter_lung_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::LungMonolithicStructureSplit(comm,fsidyn));
    }
    else if (coupling==fsi_iter_lung_monolithicfluidsplit)
    {
      fsi = Teuchos::rcp(new FSI::LungMonolithicFluidSplit(comm,fsidyn));
    }
    else if (coupling==fsi_iter_constr_monolithicfluidsplit)
    {
      fsi = Teuchos::rcp(new FSI::ConstrMonolithicFluidSplit(comm,fsidyn));
    }
    else if (coupling==fsi_iter_constr_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::ConstrMonolithicStructureSplit(comm,fsidyn));
    }
    else if (coupling==fsi_iter_mortar_monolithicstructuresplit)
    {
      fsi = Teuchos::rcp(new FSI::MortarMonolithicStructureSplit(comm,fsidyn));
    }
    else if (coupling==fsi_iter_mortar_monolithicfluidsplit)
    {
      fsi = Teuchos::rcp(new FSI::MortarMonolithicFluidSplit(comm,fsidyn));
    }
    else
    {
      dserror("Cannot find appropriate monolithic solver for coupling %d and linear strategy %d",coupling,linearsolverstrategy);
    }

    // read the restart information, set vectors and variables ---
    // be careful, dofmaps might be changed here in a Redistribute call
    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      fsi->ReadRestart(restart);
    }

    // now do the coupling setup an create the combined dofmap
    fsi->SetupSystem();

    // here we go...
    fsi->Timeloop(fsi);

    // create result tests for single fields
    DRT::Problem::Instance()->AddFieldTest(fsi->AleField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());

    // create fsi specific result test
    Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new FSI::FSIResultTest(fsi,fsidyn));
    DRT::Problem::Instance()->AddFieldTest(fsitest);

    // do the actual testing
    DRT::Problem::Instance()->TestAll(comm);

    break;
  }
  default:
  {
    // Any partitioned algorithm.

    Teuchos::RCP<FSI::Partitioned> fsi;

    INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");

    if (method==INPAR::FSI::DirichletNeumannSlideale)
      fsi = Teuchos::rcp(new FSI::DirichletNeumannSlideale(comm));
    else if (method==INPAR::FSI::DirichletNeumann)
      fsi = Teuchos::rcp(new FSI::DirichletNeumann(comm));
    else
      dserror("unsupported partitioned FSI scheme");

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(restart);
    }

    fsi->Timeloop(fsi);

    // create result tests for single fields
    DRT::Problem::Instance()->AddFieldTest(fsi->MBFluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());

    // do the actual testing
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
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("structure")->Comm();
#else
  Epetra_SerialComm comm;
#endif

  if (comm.MyPID() == 0)
  {
    cout << endl;
    cout << "       @..@    " << endl;
    cout << "      (----)      " << endl;
    cout << "     ( >__< )   " << endl;
    cout << "     ^^ ~~ ^^  " << endl;
    cout << "     _     _ _______ _______ _____" << endl;
    cout << "      \\\\__/  |______ |______   |  " << endl;
    cout << "     _/  \\\\_ |       ______| __|__" << endl;
    cout <<  endl << endl;
  }

  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  RCP<DRT::Discretization> soliddis = problem->GetDis("structure");
  soliddis->FillComplete();

  RCP<DRT::Discretization> actdis = problem->GetDis("fluid");
  actdis->FillComplete();

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();

  // now we can reserve dofs for background fluid
  int numglobalnodes = actdis->NumGlobalNodes();
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
//    const int restart = DRT::Problem::Instance()->Restart();
//    if (restart)
//    {
//      fsi->ReadRestart(restart);
//    }
//
//    // here we go...
//    fsi->Timeloop();
//
//    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
//    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());
//    DRT::Problem::Instance()->TestAll(comm);

    break;
  }
  case fsi_pseudo_structureale:
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
    dserror("Unreasonable choice");
  default:
  {
    // Any partitioned algorithm. Stable of working horses.

    Teuchos::RCP<FSI::Partitioned> fsi;

    INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn,"PARTITIONED");

    if (method==INPAR::FSI::DirichletNeumann)
    {
      fsi = Teuchos::rcp(new FSI::DirichletNeumann(comm));
    }
    else
      dserror("only Dirichlet-Neumann partitioned schemes with XFEM");

    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fsi->ReadRestart(restart);
    }

    fsi->Timeloop(fsi);

    DRT::Problem::Instance()->AddFieldTest(fsi->MBFluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
  }
  }

  Teuchos::TimeMonitor::summarize();
}


