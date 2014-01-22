/*----------------------------------------------------------------------*/
/*!
\file fsi_dyn.cpp

\brief Entry routines for FSI problems and some other problem types as well

<pre>
Maintainer: Matthias Mayr
            mayr@lnm.mw.tum.de
            http://www.mhpc.mw.tum.de
            089 - 289-10362
</pre>
*/

/*----------------------------------------------------------------------*/

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
#include "fsi_lungmonolithic_fluidsplit.H"
#include "fsi_lungmonolithic_structuresplit.H"
#include "fsi_constrmonolithic_fluidsplit.H"
#include "fsi_constrmonolithic_structuresplit.H"
#include "fsi_mortarmonolithic_fluidsplit.H"
#include "fsi_mortarmonolithic_structuresplit.H"
#include "fsi_fluidfluidmonolithic_structuresplit_nonox.H"
#include "fsi_fluidfluidmonolithic_fluidsplit_nonox.H"
#include "fsi_fluidfluidmonolithic_structuresplit.H"
#include "fsi_structureale.H"
#include "fsi_fluid_ale.H"
#include "fsi_utils.H"
#include "fsi_resulttest.H"

#include "../drt_fsi_xfem/fsi_xfem_fluid.H"
#include "../drt_fsi_xfem/fsi_xfem_monolithic.H"

#include "fs_monolithic.H"

#include "../drt_fluid_xfluid/xfluid.H"
#include "../drt_fluid_xfluid/xfluidfluid.H"

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
#include "../drt_adapter/ad_str_fsiwrapper.H"
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

  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis   = problem->GetDis("ale");
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
void fluid_xfem_drt()
{
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("structure")->Comm();

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> soliddis = problem->GetDis("structure");
  soliddis->FillComplete();

  Teuchos::RCP<DRT::Discretization> actdis = problem->GetDis("fluid");
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
    const Teuchos::ParameterList& fdyn = DRT::Problem::Instance()->FluidDynamicParams();
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn,fdyn,"fluid",false));

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
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(problem->GetDis("fluid")->Comm().Clone());

  Teuchos::RCP<DRT::Discretization> bgfluiddis = problem->GetDis("fluid");
  bgfluiddis->FillComplete();

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();


  Teuchos::RCP<DRT::Discretization> embfluiddis = problem->GetDis("xfluid");
  embfluiddis->FillComplete();

  // -------------------------------------------------------------------
  // ---------------- find MovingFluid's elements and nodes
  std::map<int, DRT::Node*> MovingFluidNodemap;
  std::map<int, Teuchos::RCP< DRT::Element> > MovingFluidelemap;
  DRT::UTILS::FindConditionObjects(*bgfluiddis, MovingFluidNodemap, MovingFluidelemap, "MovingFluid");

  // local vectors of nodes and elements of moving dis
  std::vector<int> MovingFluidNodeGIDs;
  std::vector<int> MovingFluideleGIDs;

  for( std::map<int, DRT::Node*>::iterator it = MovingFluidNodemap.begin(); it != MovingFluidNodemap.end(); ++it )
    MovingFluidNodeGIDs.push_back( it->first);

  for( std::map<int, Teuchos::RCP< DRT::Element> >::iterator it = MovingFluidelemap.begin(); it != MovingFluidelemap.end(); ++it )
    MovingFluideleGIDs.push_back( it->first);

  // --------------------------------------------------------------------------
  // ------------------ gather information for moving fluid  -------------------

  //information how many processors work at all
  std::vector<int> allproc(embfluiddis->Comm().NumProc());

  // Gather all informations from all processors
  std::vector<int> MovingFluideleGIDsall;
  std::vector<int> MovingFluidNodeGIDsall;

  //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embfluiddis->Comm().NumProc(); ++i) allproc[i] = i;

  //gathers information of MovingFluideleGIDs of all processors
  LINALG::Gather<int>(MovingFluideleGIDs,MovingFluideleGIDsall,(int)embfluiddis->Comm().NumProc(),&allproc[0],embfluiddis->Comm());

  //gathers information of MovingFluidNodeGIDs of all processors
  LINALG::Gather<int>(MovingFluidNodeGIDs,MovingFluidNodeGIDsall,(int)embfluiddis->Comm().NumProc(),&allproc[0],embfluiddis->Comm());

  //information how many processors work at all
  std::vector<int> allprocbg(bgfluiddis->Comm().NumProc());

  // -------------------------------------------------------------------------------
  // -------------- now build the nonmoving vectors from the gathered moving vectors
  std::vector<int> NonMovingFluideleGIDs;
  std::vector<int> NonMovingFluidNodeGIDs;
  for (int iele=0; iele< bgfluiddis->NumMyColElements(); ++iele)
  {
    DRT::Element* bgele = bgfluiddis->lColElement(iele);
    std::vector<int>::iterator eleiter = find(MovingFluideleGIDsall.begin(), MovingFluideleGIDsall.end(),bgele->Id() );
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
  std::vector<int> NonMovingFluideleGIDsall;
  std::vector<int> NonMovingFluidNodeGIDsall;

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

  std::vector<int> bgeleids;          // ele ids
  for (int i=0; i<bgfluiddis->NumMyRowElements(); ++i)
  {
    DRT::Element* bgele = bgfluiddis->lRowElement(i);
    int gid = bgele->Id();
    bgeleids.push_back(gid);
  }

  // Background discretization redistribution..
  Teuchos::RCP<Epetra_Map> bgroweles =  Teuchos::rcp(new Epetra_Map(-1,(int)bgeleids.size(),&bgeleids[0],0,bgfluiddis->Comm()));
  Teuchos::RCP<Epetra_Map> bgrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> bgcolnodes = Teuchos::null;

  DRT::UTILS::PartUsingParMetis(bgfluiddis,bgroweles,bgrownodes,bgcolnodes,comm,false);

  Teuchos::RCP<Epetra_Map> bgnewroweles  = Teuchos::null;
  Teuchos::RCP<Epetra_Map> bgnewcoleles  = Teuchos::null;
  bgfluiddis->BuildElementRowColumn(*bgrownodes,*bgcolnodes,bgnewroweles,bgnewcoleles);

  // export nodes and elements to the row map
  bgfluiddis->ExportRowNodes(*bgrownodes);
  bgfluiddis->ExportRowElements(*bgnewroweles);

  // export nodes and elements to the column map (create ghosting)
  bgfluiddis->ExportColumnNodes(*bgcolnodes);
  bgfluiddis->ExportColumnElements(*bgnewcoleles);

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


  std::vector<int> eleids;          // ele ids
  for (int i=0; i<embfluiddis->NumMyRowElements(); ++i)
  {
    DRT::Element* ele = embfluiddis->lRowElement(i);
    int gid = ele->Id();
    eleids.push_back(gid);
  }

  // Embedded discretization redistribution..
  Teuchos::RCP<Epetra_Map> embroweles =  Teuchos::rcp(new Epetra_Map(-1,(int)eleids.size(),&eleids[0],0,embfluiddis->Comm()));
  Teuchos::RCP<Epetra_Map> embrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> embcolnodes = Teuchos::null;

  DRT::UTILS::PartUsingParMetis(embfluiddis,embroweles,embrownodes,embcolnodes,comm,false);

  Teuchos::RCP<Epetra_Map> embnewroweles  = Teuchos::null;
  Teuchos::RCP<Epetra_Map> embnewcoleles  = Teuchos::null;
  embfluiddis->BuildElementRowColumn(*embrownodes,*embcolnodes,embnewroweles,embnewcoleles);

  // export nodes and elements to the row map
  embfluiddis->ExportRowNodes(*embrownodes);
  embfluiddis->ExportRowElements(*embnewroweles);

  // export nodes and elements to the column map (create ghosting)
  embfluiddis->ExportColumnNodes(*embcolnodes);
  embfluiddis->ExportColumnElements(*embnewcoleles);

  embfluiddis->FillComplete();
  //------------------------------------------------------------------------------

  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

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
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(DRT::Problem::Instance()->GetDis("fluid")->Comm().Clone());

  /* |--str dofs--|--bgfluid dofs--|--embfluid dofs--|--ale dofs--|-> */

  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> structdis = problem->GetDis("structure");
  structdis->FillComplete();

  Teuchos::RCP<DRT::Discretization> bgfluiddis = problem->GetDis("fluid");
  bgfluiddis->FillComplete();

  // reserve max size of dofs for the background fluid
  const Teuchos::ParameterList xdyn = DRT::Problem::Instance()->XFEMGeneralParams();

  Teuchos::RCP<DRT::Discretization> embfluiddis = problem->GetDis("xfluid");
  embfluiddis->FillComplete();

  // -------------------------------------------------------------------
  // ---------------- find MovingFluid's elements and nodes
  std::map<int, DRT::Node*> MovingFluidNodemap;
  std::map<int, Teuchos::RCP< DRT::Element> > MovingFluidelemap;
  DRT::UTILS::FindConditionObjects(*bgfluiddis, MovingFluidNodemap, MovingFluidelemap, "MovingFluid");

  // local vectors of nodes and elements of moving dis
  std::vector<int> MovingFluidNodeGIDs;
  std::vector<int> MovingFluideleGIDs;

  for( std::map<int, DRT::Node*>::iterator it = MovingFluidNodemap.begin(); it != MovingFluidNodemap.end(); ++it )
    MovingFluidNodeGIDs.push_back( it->first);

  for( std::map<int, Teuchos::RCP< DRT::Element> >::iterator it = MovingFluidelemap.begin(); it != MovingFluidelemap.end(); ++it )
    MovingFluideleGIDs.push_back( it->first);

  // --------------------------------------------------------------------------
  // ------------------ gather information for moving fluid -------------------
  //information how many processors work at all
  std::vector<int> allproc(embfluiddis->Comm().NumProc());

  // Gather all informations from all processors
  std::vector<int> MovingFluideleGIDsall;
  std::vector<int> MovingFluidNodeGIDsall;

  //in case of n processors allproc becomes a vector with entries (0,1,...,n-1)
  for (int i=0; i<embfluiddis->Comm().NumProc(); ++i) allproc[i] = i;

  //gathers information of MovingFluideleGIDs of all processors
  LINALG::Gather<int>(MovingFluideleGIDs,MovingFluideleGIDsall,(int)embfluiddis->Comm().NumProc(),&allproc[0],embfluiddis->Comm());

  //gathers information of MovingFluidNodeGIDs of all processors
  LINALG::Gather<int>(MovingFluidNodeGIDs,MovingFluidNodeGIDsall,(int)embfluiddis->Comm().NumProc(),&allproc[0],embfluiddis->Comm());

  // -------------------------------------------------------------------------------
  // -------------- now build the nonmoving vectors from the gathered moving vectors
  std::vector<int> NonMovingFluideleGIDs;
  std::vector<int> NonMovingFluidNodeGIDs;
  for (int iele=0; iele< bgfluiddis->NumMyColElements(); ++iele)
  {
    DRT::Element* bgele = bgfluiddis->lColElement(iele);
    std::vector<int>::iterator eleiter = find(MovingFluideleGIDsall.begin(), MovingFluideleGIDsall.end(),bgele->Id() );
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
  std::vector<int> allprocbg(bgfluiddis->Comm().NumProc());

  std::vector<int> NonMovingFluideleGIDsall;
  std::vector<int> NonMovingFluidNodeGIDsall;

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
  IO::cout << " maxNumMyReservedDofs " << maxNumMyReservedDofs << IO::endl;
  bgfluiddis->ReplaceDofSet(maxdofset,true);
  bgfluiddis->FillComplete();

  std::vector<int> bgeleids;          // ele ids
  for (int i=0; i<bgfluiddis->NumMyRowElements(); ++i)
  {
    DRT::Element* bgele = bgfluiddis->lRowElement(i);
    int gid = bgele->Id();
    bgeleids.push_back(gid);
  }

  // Background discretization redistribution..
  Teuchos::RCP<Epetra_Map> bgroweles =  Teuchos::rcp(new Epetra_Map(-1,(int)bgeleids.size(),&bgeleids[0],0,bgfluiddis->Comm()));
  Teuchos::RCP<Epetra_Map> bgrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> bgcolnodes = Teuchos::null;

  DRT::UTILS::PartUsingParMetis(bgfluiddis,bgroweles,bgrownodes,bgcolnodes,comm,false);

  Teuchos::RCP<Epetra_Map> bgnewroweles  = Teuchos::null;
  Teuchos::RCP<Epetra_Map> bgnewcoleles  = Teuchos::null;
  bgfluiddis->BuildElementRowColumn(*bgrownodes,*bgcolnodes,bgnewroweles,bgnewcoleles);

  // export nodes and elements to the row map
  bgfluiddis->ExportRowNodes(*bgrownodes);
  bgfluiddis->ExportRowElements(*bgnewroweles);

  // export nodes and elements to the column map (create ghosting)
  bgfluiddis->ExportColumnNodes(*bgcolnodes);
  bgfluiddis->ExportColumnElements(*bgnewcoleles);

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

  std::vector<int> eleids;          // ele ids
  for (int i=0; i<embfluiddis->NumMyRowElements(); ++i)
  {
    DRT::Element* ele = embfluiddis->lRowElement(i);
    int gid = ele->Id();
    eleids.push_back(gid);
  }

  // Embedded discretization redistribution..
  Teuchos::RCP<Epetra_Map> embroweles =  Teuchos::rcp(new Epetra_Map(-1,(int)eleids.size(),&eleids[0],0,embfluiddis->Comm()));
  Teuchos::RCP<Epetra_Map> embrownodes = Teuchos::null;
  Teuchos::RCP<Epetra_Map> embcolnodes = Teuchos::null;

  DRT::UTILS::PartUsingParMetis(embfluiddis,embroweles,embrownodes,embcolnodes,comm,false);

  Teuchos::RCP<Epetra_Map> embnewroweles  = Teuchos::null;
  Teuchos::RCP<Epetra_Map> embnewcoleles  = Teuchos::null;
  embfluiddis->BuildElementRowColumn(*embrownodes,*embcolnodes,embnewroweles,embnewcoleles);

  // export nodes and elements to the row map
  embfluiddis->ExportRowNodes(*embrownodes);
  embfluiddis->ExportRowElements(*embnewroweles);

  // export nodes and elements to the column map (create ghosting)
  embfluiddis->ExportColumnNodes(*embcolnodes);
  embfluiddis->ExportColumnElements(*embnewcoleles);
  embfluiddis->FillComplete();

  embfluiddis->FillComplete();

  //------------------------------------------------------------------------------

  // create ale elements if the ale discretization is empty
  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");
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

      // create fsi specific result test
      Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new FSI::FSIResultTest(fsi,fsidyn));
      DRT::Problem::Instance()->AddFieldTest(fsitest);

      // do the actual testing
      DRT::Problem::Instance()->TestAll(*comm);

      break;
    }
    case fsi_iter_fluidfluid_monolithicfluidsplit:
    {
      Teuchos::RCP<FSI::MonolithicNoNOX> fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicFluidSplitNoNOX(*comm,fsidyn));

       //Setup the coupling, create the combined dofmap
      fsi->SetupSystem();

      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        // Read the restart information, set vectors and variables
        fsi->ReadRestart(restart);
      }

      fsi->Timeloop();

      DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());

      // create fsi specific result test
      Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new FSI::FSIResultTest(fsi,fsidyn));
      DRT::Problem::Instance()->AddFieldTest(fsitest);

      // do the actual testing
      DRT::Problem::Instance()->TestAll(*comm);

      break;
    }
    case fsi_iter_fluidfluid_monolithicstructuresplit_nox:
    {
      Teuchos::RCP<FSI::Monolithic> fsi = Teuchos::rcp(new FSI::FluidFluidMonolithicStructureSplit(*comm,fsidyn));

      // now do the coupling setup an create the combined dofmap
      fsi->SetupSystem();

      const int restart = DRT::Problem::Instance()->Restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        fsi->ReadRestart(restart);
      }

      // here we go...
      fsi->Timeloop(fsi);

      DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
      DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());
      // create fsi specific result test
      Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new FSI::FSIResultTest(fsi,fsidyn));
      DRT::Problem::Instance()->AddFieldTest(fsitest);
      // do the actual testing
      DRT::Problem::Instance()->TestAll(*comm);

      break;
    }
    default:
    {
      // Any partitioned algorithm
      Teuchos::RCP<FSI::Partitioned> fsi;

      INPAR::FSI::PartitionedCouplingMethod method =
        DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn.sublist("PARTITIONED SOLVER"),"PARTITIONED");

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

      // do the actual testing
      DRT::Problem::Instance()->TestAll(*comm);

      break;
    }
  }
}

/*----------------------------------------------------------------------*/
// entry point for (pure) free surface in DRT
/*----------------------------------------------------------------------*/
void fluid_freesurf_drt()
{
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();

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
  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  Teuchos::RCP<DRT::Discretization> aledis = problem->GetDis("ale");

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

  const Teuchos::ParameterList& fsidyn = problem->FSIDynamicParams();
  const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
  {

    INPAR::FSI::LinearBlockSolver linearsolverstrategy
      = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono,"LINEARBLOCKSOLVER");

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

  const Epetra_Comm& comm = problem->GetDis("structure")->Comm();

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
    // monolithic solver settings
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

    Teuchos::RCP<FSI::Monolithic> fsi;

    INPAR::FSI::LinearBlockSolver linearsolverstrategy
      = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono,"LINEARBLOCKSOLVER");

    // call constructor to initialize the base class
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

    // calculate errors in comparison to analytical solution
    fsi->FluidField().CalculateError();

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
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn.sublist("PARTITIONED SOLVER"),"PARTITIONED");

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

    break;
  }
  }

  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
}




/*----------------------------------------------------------------------*/
// entry point for FSI using XFEM in DRT
/*----------------------------------------------------------------------*/
void xfsi_drt()
{
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("structure")->Comm();

  if (comm.MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "       @..@    " << std::endl;
    std::cout << "      (----)      " << std::endl;
    std::cout << "     ( >__< )   " << std::endl;
    std::cout << "     ^^ ~~ ^^  " << std::endl;
    std::cout << "     _     _ _______ _______ _____" << std::endl;
    std::cout << "      \\\\__/  |______ |______   |  " << std::endl;
    std::cout << "     _/  \\\\_ |       ______| __|__" << std::endl;
    std::cout <<  std::endl << std::endl;
  }

  DRT::Problem* problem = DRT::Problem::Instance();
  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  Teuchos::RCP<DRT::Discretization> soliddis = problem->GetDis("structure");
  soliddis->FillComplete();

  Teuchos::RCP<DRT::Discretization> fluiddis = problem->GetDis("fluid");
  fluiddis->FillComplete();

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();

  // now we can reserve dofs for background fluid
  int numglobalnodes = fluiddis->NumGlobalNodes();
  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  fluiddis->ReplaceDofSet(maxdofset,true);
  fluiddis->FillComplete();

  // print all dofsets
  fluiddis->GetDofSetProxy()->PrintAllDofsets(fluiddis->Comm());

  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_xfem_monolithic:
  {
    // monolithic solver settings
    const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

    INPAR::FSI::LinearBlockSolver linearsolverstrategy
      = DRT::INPUT::IntegralValue<INPAR::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

    if (linearsolverstrategy!=INPAR::FSI::PreconditionedKrylov)
      dserror("Only Newton-Krylov scheme with XFEM fluid");

    // create the MonolithicXFEM object that does the whole work
    Teuchos::RCP<FSI::AlgorithmXFEM> fsi = Teuchos::rcp(new FSI::MonolithicXFEM(comm, fsidyn));

    // setup the system (block-DOF-row maps, systemmatrix etc.) for the monolithic XFEM system
    fsi->SetupSystem();

    // read the restart information, set vectors and variables ---
    // be careful, dofmaps might be changed here in a Redistribute call
    const int restart = DRT::Problem::Instance()->Restart();
    if (restart)
    {
      fsi->ReadRestart(restart);
    }

    // here we go...
    fsi->Timeloop();

    DRT::Problem::Instance()->AddFieldTest(fsi->FluidField()->CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(fsi->StructureField()->CreateFieldTest());

//    // create FSI specific result test
//    Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new FSI::FSIResultTest(fsi,fsidyn));
//    DRT::Problem::Instance()->AddFieldTest(fsitest);

    // do the actual testing
    DRT::Problem::Instance()->TestAll(comm);

    break;
  }
  case fsi_pseudo_structureale:
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
    dserror("Unreasonable choice");
    break;
  default:
  {
    // Any partitioned algorithm. Stable of working horses.

    Teuchos::RCP<FSI::Partitioned> fsi;

    INPAR::FSI::PartitionedCouplingMethod method =
      DRT::INPUT::IntegralValue<INPAR::FSI::PartitionedCouplingMethod>(fsidyn.sublist("PARTITIONED SOLVER"),"PARTITIONED");

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

    break;
  }
  }

  Teuchos::TimeMonitor::summarize();
}


