/*!----------------------------------------------------------------------
\file fluid_dyn_nln_drt.cpp
\brief Main control routine for all fluid (in)stationary solvers,

     including instationary solvers based on

     o one-step-theta time-integration scheme

     o two-step BDF2 time-integration scheme
       (with potential one-step-theta start algorithm)

     o generalized-alpha time-integration scheme

     and stationary solver.

<pre>
Maintainer: Ursula Rasthofer
            rasthofer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15236
</pre>

 *----------------------------------------------------------------------*/
/*
#include <ctime>
#include <cstdlib>
#include <iostream>
*/

#include "fluid_dyn_nln_drt.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/ad_fld_base_algorithm.H"
#include "../drt_fluid_turbulence/turbulent_flow_algorithm.H"
#include "../drt_lib/drt_condition_utils.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_dofset_fixed_size.H"
#include "../drt_lib/drt_utils_parmetis.H"
#include "../drt_lib/drt_discret_xfem.H"

/*----------------------------------------------------------------------*
 * Main control routine for fluid including various solvers:
 *
 *        o instationary one-step-theta
 *        o instationary BDF2
 *        o instationary generalized-alpha (two versions)
 *        o stationary
 *
 *----------------------------------------------------------------------*/
void dyn_fluid_drt(const int restart)
{
  // create a communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis("fluid")->Comm();

  // access to some parameter lists
  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  // prepares a turbulent flow simulation with generation of turbulent inflow during the
  // actual simulation
  // this is done in two steps
  // 1. computation of inflow until it reaches a fully turbulent state
  // 2. computation of the main problem after restart
  // Remark: we restart the simulation to save procs!
  if ((DRT::INPUT::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"),"TURBULENTINFLOW")==true) and
     (restart<fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
  {
    if (comm.MyPID()==0)
    {
      std::cout << "#-----------------------------------------------#" << std::endl;
      std::cout << "#      ENTER TURBULENT INFLOW COMPUTATION       #" << std::endl;
      std::cout << "#-----------------------------------------------#" << std::endl;
    }

    // create instance of fluid turbulent flow algorithm
    Teuchos::RCP<FLD::TurbulentFlowAlgorithm> turbfluidalgo = Teuchos::rcp(new FLD::TurbulentFlowAlgorithm(comm,fdyn));

    // read the restart information, set vectors and variables
    if (restart) turbfluidalgo->ReadRestart(restart);

    // run simulation for a separate part of the complete domain to get turbulent flow in it
    // after restart a turbulent inflow profile is computed in the separate inflow section and
    // transferred as dirichlet boundary condition to the problem domain of interest
    // this finally allows to get high quality turbulent inflow conditions during simulation of the
    // actual flow
    turbfluidalgo->TimeLoop();

    // perform result tests if required
    DRT::Problem::Instance()->AddFieldTest(turbfluidalgo->DoResultCheck());
    DRT::Problem::Instance()->TestAll(comm);
  }
  // solve a simple fluid problem
  else
  {
    // create instance of fluid basis algorithm
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn,fdyn,"fluid",false));

    // read the restart information, set vectors and variables
    if (restart) fluidalgo->FluidField()->ReadRestart(restart);

    // run the simulation
//    fluidalgo->FluidField()->TimeLoop();
    fluidalgo->FluidField()->Integrate();

    // perform result tests if required
    DRT::Problem::Instance()->AddFieldTest(fluidalgo->FluidField()->CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
  }

  // have fun with your results!
  return;

} // end of dyn_fluid_drt()

//----------------------------------------------------------------------------
// main routine for fluid_fluid problems
//-------------------------------------------------------------------------
void fluid_fluid_drt(const int restart)
{
  DRT::Problem* problem = DRT::Problem::Instance();

  Teuchos::RCP<DRT::Discretization> embfluiddis = problem->GetDis("fluid");
  embfluiddis->FillComplete();

  // create a communicator
  Teuchos::RCP<Epetra_Comm> comm = Teuchos::rcp(DRT::Problem::Instance()->GetDis("xfluid")->Comm().Clone());

  Teuchos::RCP<DRT::DiscretizationXFEM> bgfluiddis =  Teuchos::rcp_dynamic_cast<DRT::DiscretizationXFEM>(problem->GetDis("xfluid"));
  bgfluiddis->FillComplete();

  const Teuchos::ParameterList xdyn = DRT::Problem::Instance()->XFEMGeneralParams();

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

  // Gather all informations from all processors
  //information how many processors work at all
  std::vector<int> allproc(embfluiddis->Comm().NumProc());

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

#if defined(PARALLEL) && defined(PARMETIS)
  DRT::UTILS::PartUsingParMetis(embfluiddis,embroweles,embrownodes,embcolnodes,comm,false);
#endif

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
  int maxNumMyReservedDofsperNode = (xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofsperNode,numglobalnodes));
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

#if defined(PARALLEL) && defined(PARMETIS)
  DRT::UTILS::PartUsingParMetis(bgfluiddis,bgroweles,bgrownodes,bgcolnodes,comm,false);
#endif

  Teuchos::RCP<Epetra_Map> bgnewroweles  = Teuchos::null;
  Teuchos::RCP<Epetra_Map> bgnewcoleles  = Teuchos::null;
  bgfluiddis->BuildElementRowColumn(*bgrownodes,*bgcolnodes,bgnewroweles,bgnewcoleles);

  // export nodes and elements to the row map
  bgfluiddis->ExportRowNodes(*bgrownodes);
  bgfluiddis->ExportRowElements(*bgnewroweles);

  // export nodes and elements to the column map (create ghosting)
  bgfluiddis->ExportColumnNodes(*bgcolnodes);
  bgfluiddis->ExportColumnElements(*bgnewcoleles);

  bgfluiddis->InitialFillComplete();

  //------------------------------------------------------------------------------

  // access to some parameter lists
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  // create instance of fluid basis algorithm
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo = Teuchos::rcp(new ADAPTER::FluidBaseAlgorithm(fdyn,fdyn,"fluid",false));

  // read the restart information, set vectors and variables
  if (restart)
    fluidalgo->FluidField()->ReadRestart(restart);


  // run the simulation
//  fluidalgo->FluidField()->TimeLoop();
  fluidalgo->FluidField()->Integrate();

  // perform result tests if required
  DRT::Problem::Instance()->AddFieldTest(fluidalgo->FluidField()->CreateFieldTest());
  DRT::Problem::Instance()->TestAll(*comm);
  return;


} // end of fluid_fluid_drt()

