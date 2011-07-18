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
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif
  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  RCP<DRT::Discretization> bgfluiddis = problem->Dis(genprob.numff,0);
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
      const DRT::Node*const* bgelenodes = bgele->Nodes();
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

  // delete elements and nodes
  for(int mv=0; mv<MovingFluideleGIDs.size(); ++mv)
    bgfluiddis->DeleteElement(MovingFluideleGIDs.at(mv));

  for(int mv=0; mv<MovingFluidNodeGIDs.size(); ++mv)
    bgfluiddis->DeleteNode(MovingFluidNodeGIDs.at(mv));

  bgfluiddis->FillComplete();

  for(int nmv=0; nmv<NonMovingFluideleGIDs.size(); ++nmv)
    embfluiddis->DeleteElement(NonMovingFluideleGIDs.at(nmv));

  for(int nmv=0; nmv<NonMovingFluidNodeGIDs.size(); ++nmv)
    embfluiddis->DeleteNode(NonMovingFluidNodeGIDs.at(nmv));

  embfluiddis->FillComplete();

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
  DRT::Problem::Instance()->TestAll(comm);

  return;

} // end of fluid_fluid_drt()

#endif  // #ifdef CCADISCRET
