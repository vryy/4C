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
Maintainer: Peter Gamnitzer
            gamnitzer@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
</pre>

 *----------------------------------------------------------------------*/
#ifdef CCADISCRET
/*
#include <ctime>
#include <cstdlib>
#include <iostream>
*/

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "fluid_dyn_nln_drt.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_fluid_base_algorithm.H"
#include "turbulent_flow_algorithm.H"
#include "../drt_lib/drt_utils_createdis.H"
#include "../linalg/linalg_utils.H"
#include "../drt_lib/drt_dofset_fixed_size.H"

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
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

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
    Teuchos::RCP<FLD::TurbulentFlowAlgorithm> turbfluidalgo = rcp(new FLD::TurbulentFlowAlgorithm(comm,fdyn));

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
    Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo = rcp(new ADAPTER::FluidBaseAlgorithm(fdyn,false));

    // read the restart information, set vectors and variables
    if (restart) fluidalgo->FluidField().ReadRestart(restart);

    // run the simulation
    fluidalgo->FluidField().TimeLoop();

    // perform result tests if required
    DRT::Problem::Instance()->AddFieldTest(fluidalgo->FluidField().CreateFieldTest());
    DRT::Problem::Instance()->TestAll(comm);
  }

  // have fun with your results!
  return;

} // end of dyn_fluid_drt()
//----------------------------------------------------------------------------
// main routine for fluid_fluid problems
//-------------------------------------------------------------------------
void fluid_fluid_drt()
{
  // create a communicator
  #ifdef PARALLEL
   RCP<Epetra_Comm> comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    Epetra_SerialComm comm;
  #endif

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> bgfluiddis = problem->Dis(genprob.numff,0);

  const Teuchos::ParameterList xdyn = DRT::Problem::Instance()->XFEMGeneralParams();

  // compute numnode
  int numglobalnodes = 0;
  int numlocalnodes = bgfluiddis->NumMyColNodes();
  (bgfluiddis->Comm()).SumAll(&numlocalnodes,&numglobalnodes,1);

  // reserve max size of dofs for the background fluid
  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  bgfluiddis->ReplaceDofSet(maxdofset);

  bgfluiddis->FillComplete();

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

  vector<string> conditions_to_copy_mf;
  conditions_to_copy_mf.push_back("MovingFluid");
  Teuchos::RCP<DRT::Discretization> MovingFluiddis = DRT::UTILS::CreateDiscretizationFromCondition(bgfluiddis, "MovingFluid",
           "MovingFluid", "VELE3", conditions_to_copy_mf);

  // delete embedded fluid's node and elements from the background fluid
  vector<int> MovingFluideleGIDs;
  for (int iele=0; iele< MovingFluiddis->NumMyColElements(); ++iele)
  {
    DRT::Element* MovingFluidele = MovingFluiddis->lColElement(iele);
    MovingFluideleGIDs.push_back(MovingFluidele->Id());
  }

  vector<int> MovingFluidNodeGIDs;
  for (int node=0; node<MovingFluiddis->NumMyColNodes(); node++)
  {
    DRT::Node*  MovingFluidnode = MovingFluiddis->lColNode(node);
    MovingFluidNodeGIDs.push_back(MovingFluidnode->Id());
  }

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

  // access to some parameter lists
  //const Teuchos::ParameterList& probtype = DRT::Problem::Instance()->ProblemTypeParams();
  const Teuchos::ParameterList& fdyn     = DRT::Problem::Instance()->FluidDynamicParams();

  // create instance of fluid basis algorithm
  Teuchos::RCP<ADAPTER::FluidBaseAlgorithm> fluidalgo = rcp(new ADAPTER::FluidBaseAlgorithm(fdyn,false));

  // read the restart information, set vectors and variables
  //if (restart) fluidalgo->FluidField().ReadRestart(restart);

  // run the simulation
  fluidalgo->FluidField().TimeLoop();

  // perform result tests if required
  DRT::Problem::Instance()->AddFieldTest(fluidalgo->FluidField().CreateFieldTest());
  DRT::Problem::Instance()->TestAll(*comm);
  return;

} // end of fluid_fluid_drt()

#endif  // #ifdef CCADISCRET
