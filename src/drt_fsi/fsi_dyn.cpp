
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

#include "../drt_fs3i/fs3i.H"
#include "../drt_fs3i/fs3i_1wc.H"
#include "../drt_fs3i/fs3i_biofilm_growth.H"

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
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  RefCountPtr<DRT::Discretization> aledis = DRT::Problem::Instance()->Dis(genprob.numaf,0);
  if (!aledis->Filled()) aledis->FillComplete();

  // create ale elements if the ale discretization is empty
  if (aledis->NumGlobalNodes()==0)
  {
    {
      RCP<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(genprob.numff,0);

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
  Epetra_MpiComm comm(MPI_COMM_WORLD);
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
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> soliddis = problem->Dis(genprob.numsf,0);
  soliddis->FillComplete();

  RCP<DRT::Discretization> actdis = problem->Dis(genprob.numff,0);

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();

  // compute number of nodes
  int numglobalnodes = 0.0;
  int numlocalnodes = actdis->NumMyColNodes();
  (actdis->Comm()).SumAll(&numlocalnodes,&numglobalnodes,1);

  // reserve max size of dofs for the fluid dis
  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  actdis->ReplaceDofSet(maxdofset);

  actdis->FillComplete();

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
  Teuchos::RCP<LINALG::Solver> solver =
    Teuchos::rcp(new LINALG::Solver(problem->FluidSolverParams(),
                                    actdis->Comm(),
                                    problem->ErrorFile()->Handle()));
  //actdis->ComputeNullSpaceIfNecessary(solver->Params());

  FLD::XFluid fluid( actdis,soliddis,*solver,problem->FluidDynamicParams(), xdyn);
  fluid.Integrate();

//   Teuchos::RCP<FSI::FluidXFEMAlgorithm> xfluid = Teuchos::rcp(new FSI::FluidXFEMAlgorithm(comm));
//   if (genprob.restart)
//   {
//     // read the restart information, set vectors and variables
//     xfluid->ReadRestart(genprob.restart);
//   }
//   xfluid->Timeloop();

  DRT::Problem::Instance()->AddFieldTest(Teuchos::rcp(new FLD::XFluidResultTest2(&fluid)));
  DRT::Problem::Instance()->TestAll(comm);
}

/*----------------------------------------------------------------------*/
// entry point for Fluid-Fluid based  on XFEM in DRT
/*----------------------------------------------------------------------*/
void fluid_fluid_ale_drt()
{
  // create a communicator
  #ifdef PARALLEL
   RCP<Epetra_Comm> comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    Epetra_SerialComm comm;
  #endif

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> bgfluiddis = problem->Dis(genprob.numff,0);

  const Teuchos::ParameterList xdyn = problem->XFEMGeneralParams();

  // compute number of nodes
  int numglobalnodes = 0.0;
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

  // -------------------------------------------------------------------
  // create a solver
  // -------------------------------------------------------------------
//   Teuchos::RCP<LINALG::Solver> solver =
//     Teuchos::rcp(new LINALG::Solver(problem->FluidSolverParams(),
//                                     bgfluiddis->Comm(),
//                                     problem->ErrorFile()->Handle()));
//   FLD::XFluidFluid fluid(bgfluiddis,embfluiddis,*solver,problem->FluidDynamicParams());
//   fluid.IntegrateFluidFluid();
//   DRT::Problem::Instance()->AddFieldTest(Teuchos::rcp(new FLD::XFluidFluidResultTest(&fluid)));
//   DRT::Problem::Instance()->TestAll(comm);
}

/*----------------------------------------------------------------------*/
// entry point for fluid-fluid-fsi based on XFEM
/*----------------------------------------------------------------------*/
void fluid_fluid_fsi_drt()
{
 // create a communicator
  #ifdef PARALLEL
   RCP<Epetra_Comm> comm = rcp(new Epetra_MpiComm(MPI_COMM_WORLD));
  #else
    Epetra_SerialComm comm;
  #endif

  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> structdis = problem->Dis(genprob.numsf,0);
  structdis->FillComplete();

  RCP<DRT::Discretization> bgfluiddis = problem->Dis(genprob.numff,0);

  // reserve max size of dofs for the background fluid
  const Teuchos::ParameterList xdyn = DRT::Problem::Instance()->XFEMGeneralParams();

  // compute number of nodes
  int numglobalnodes = 0.0;
  int numlocalnodes = bgfluiddis->NumMyColNodes();
  (bgfluiddis->Comm()).SumAll(&numlocalnodes,&numglobalnodes,1);

  int maxNumMyReservedDofs = numglobalnodes*(xdyn.get<int>("MAX_NUM_DOFSETS"))*4;
  Teuchos::RCP<DRT::FixedSizeDofSet> maxdofset = Teuchos::rcp(new DRT::FixedSizeDofSet(maxNumMyReservedDofs));
  bgfluiddis->ReplaceDofSet(maxdofset);

  bgfluiddis->FillComplete();

  RCP<DRT::Discretization> embfluiddis = problem->Dis(genprob.numff,1);
  embfluiddis->FillComplete();

  RCP<DRT::Discretization> aledis = problem->Dis(genprob.numaf,0);
  aledis->FillComplete();

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

  structdis->FillComplete();
  bgfluiddis->FillComplete();
  embfluiddis->FillComplete();
  aledis->FillComplete();

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
    fsi->Timeloop(fsi);

//     DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
//     DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
//     DRT::Problem::Instance()->TestAll(comm);
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
  Epetra_MpiComm comm(MPI_COMM_WORLD);
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
  Epetra_MpiComm comm(MPI_COMM_WORLD);
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
    Teuchos::RCP<FSI::Monolithic> fsi;

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
  Epetra_MpiComm comm(MPI_COMM_WORLD);
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

  int coupling = DRT::INPUT::IntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_monolithicxfem:
  {
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


/*----------------------------------------------------------------------*/
// entry point for gas exchange lung fsi
/*----------------------------------------------------------------------*/
void fsi_lung_gas()
{
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
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
  problem->Dis(genprob.numscatra,0)->FillComplete();
  problem->Dis(genprob.numscatra,1)->FillComplete();

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

  // access the fluid discretization
  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(1,0);
  // access the structure discretization
  RefCountPtr<DRT::Discretization> structdis = DRT::Problem::Instance()->Dis(0,0);
  // access the fluid scatra discretization
  RefCountPtr<DRT::Discretization> fluidscatradis = DRT::Problem::Instance()->Dis(3,0);
  // access the fluid structure discretization
  RefCountPtr<DRT::Discretization> structscatradis = DRT::Problem::Instance()->Dis(3,1);

  // get material map for the transport elements
  std::map<std::pair<string,string>,std::map<int,int> > clonefieldmatmap = DRT::Problem::Instance()->ClonedMaterialMap();
  if (clonefieldmatmap.size() < 2)
    dserror("at least 2 matlists needed for lung gas exchange");

  // FLUID SCATRA
  // we use the fluid discretization as layout for the scalar transport discretization
  if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

  // create fluid scatra elements if the fluid scatra discretization is empty
  if (fluidscatradis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);

    // create the fluid scatra discretization
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );

      std::pair<string,string> key("fluid","scatra1");
      std::map<int,int> fluidmatmap = clonefieldmatmap[key];

      clonewizard->CreateMatchingDiscretization(fluiddis,fluidscatradis,fluidmatmap);
    }
    if (comm.MyPID()==0)
      cout<<"Created scalar transport discretization from fluid field in...."
          <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("Fluid AND ScaTra discretization present. This is not supported.");

  // STRUCTURE SCATRA
  // we use the structure discretization as layout for the scalar transport discretization
  if (structdis->NumGlobalNodes()==0) dserror("Structure discretization is empty!");

  // create structure scatra elements if the structure scatra discretization is empty
  if (structscatradis->NumGlobalNodes()==0)
  {
    Epetra_Time time(comm);

    // create the structure scatra discretization
    {
      Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard =
        Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );

      std::pair<string,string> key("structure","scatra2");
      std::map<int,int> structmatmap = clonefieldmatmap[key];

      clonewizard->CreateMatchingDiscretization(structdis,structscatradis,structmatmap);
    }
    if (comm.MyPID()==0)
      cout<<"Created scalar transport discretization from structure field in...."
          <<time.ElapsedTime() << " secs\n\n";
  }
  else
    dserror("Structure AND ScaTra discretization present. This is not supported.");

  const Teuchos::ParameterList& fsidyn   = problem->FSIDynamicParams();

  int coupling = Teuchos::getIntegralValue<int>(fsidyn,"COUPALGO");
  switch (coupling)
  {
  case fsi_iter_monolithicfluidsplit:
  case fsi_iter_monolithicstructuresplit:
  {
    Teuchos::RCP<FSI::Monolithic> fsi;

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
    else
    {
      dserror("Cannot find appropriate monolithic solver for coupling %d and linear strategy %d",coupling,linearsolverstrategy);
    }

    // access the biofilm-specific parameter list
        const Teuchos::ParameterList& biofilmcontrol = DRT::Problem::Instance()->BIOFILMControlParams();
    	const int surfgrowth = DRT::INPUT::IntegralValue<int>(biofilmcontrol,"SURFACEGROWTH");

    	if (surfgrowth==1)
    		{
    		cout << "\n\n---------> Surface growth algorithm <---------"<< endl;

    		// create an FSI::BiofilmGrowth instance
    		Teuchos::RCP<FS3I::BiofilmGrowth> biofilmgrowth = Teuchos::rcp(new FS3I::BiofilmGrowth(fsi, comm, fsidyn));
    		biofilmgrowth->ReadRestart();
    		//biofilmgrowth->SetupFSISystem();
    		biofilmgrowth->OutTimeloop();
    		}
    		else
    		{

			Teuchos::RCP<FS3I::FS3I_1WC> lungscatra = Teuchos::rcp(new FS3I::FS3I_1WC(fsi));

			lungscatra->ReadRestart();
			lungscatra->SetupFSISystem();
			lungscatra->Timeloop();
    		}


		DRT::Problem::Instance()->AddFieldTest(fsi->FluidField().CreateFieldTest());
		DRT::Problem::Instance()->AddFieldTest(fsi->StructureField().CreateFieldTest());
		DRT::Problem::Instance()->TestAll(comm);
		break;
  }
  default:
  {
    dserror("coupling algorithm not yet implemented for lung fsi gas exchange");
    break;
  }
  }

  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
}

#endif
