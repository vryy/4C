/*----------------------------------------------------------------------*/
/*!
\file scatra_dyn.cpp
\brief entry point for (passive) scalar transport problems

<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15245
</pre>
*/
/*----------------------------------------------------------------------*/
#ifdef CCADISCRET

#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_StandardParameterEntryValidators.hpp>

#ifdef PARALLEL
#include <mpi.h>
#endif

#include "scatra_dyn.H"
#include "passive_scatra_algorithm.H"
#include "scatra_utils.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"
#include "scatra_resulttest.H"


/*----------------------------------------------------------------------*
 * Main control routine for scalar transport problems, icl. various solvers
 *
 *        o Laplace-/ Poisson equation (zero velocity field)
 *          (with linear and nonlinear boundary conditons)
 *        o transport of passive scalar in velocity field given by spatial function
 *        o transport of passive scalar in velocity field given by Navier-Stokes
 *          (one-way coupling)
 *
 *----------------------------------------------------------------------*/
void scatra_dyn(int disnumff, int disnumscatra, int restart)
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // access the problem-specific parameter list
  const Teuchos::ParameterList& scatradyn     = DRT::Problem::Instance()->ScalarTransportDynamicParams();

  // access the fluid discretization
  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(disnumff,0);
  if (!fluiddis->Filled()) fluiddis->FillComplete();
  
  // access the scatra discretization
  RefCountPtr<DRT::Discretization> scatradis = DRT::Problem::Instance()->Dis(disnumscatra,0);
  if (!scatradis->Filled()) scatradis->FillComplete();

  // set velocity field
  int veltype = Teuchos::getIntegralValue<int>(scatradyn,"VELOCITYFIELD");
  switch (veltype)
  {
    case 0:  // zero  (see case 1)
    case 1:  // function
    {
      // we directly use the elements from the scalar transport elements section
      if (scatradis->NumGlobalNodes()==0)
        dserror("No elements in the ---TRANSPORT ELEMENTS section");

      // create instance of scalar transport basis algorithm (empty fluid discretization)
      Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = rcp(new ADAPTER::ScaTraBaseAlgorithm(scatradyn));

      // read the restart information, set vectors and variables
      if (restart) scatraonly->ScaTraField().ReadRestart(restart);

      // set velocity field 
      //(this is done only once. Time-dependent velocity fields are not supported)
      (scatraonly->ScaTraField()).SetVelocityField(veltype,scatradyn.get<int>("VELFUNCNO"));

      // enter time loop to solve problem with given convective velocity
      (scatraonly->ScaTraField()).TimeLoop();

      // perform the result test if required
      DRT::ResultTestManager testmanager(comm);
      testmanager.AddFieldTest(scatraonly->CreateScaTraFieldTest());
      testmanager.TestAll();

      break;
    }
    case 2:  // Navier_Stokes
    {
      // we use the fluid discretization as layout for the scalar transport discretization
      if (fluiddis->NumGlobalNodes()==0)
        dserror("No fluid discretization found!");
 
      // create scatra elements if the scatra discretization is empty
      if (scatradis->NumGlobalNodes()==0)
      {
        Epetra_Time time(comm);
        std::map<string,string> conditions_to_copy;
        conditions_to_copy.insert(pair<string,string>("TransportDirichlet","Dirichlet"));
        conditions_to_copy.insert(pair<string,string>("TransportPointNeumann","PointNeumann"));
        conditions_to_copy.insert(pair<string,string>("TransportLineNeumann","LineNeumann"));
        conditions_to_copy.insert(pair<string,string>("TransportSurfaceNeumann","SurfaceNeumann"));
        conditions_to_copy.insert(pair<string,string>("TransportVolumeNeumann","VolumeNeumann"));
        conditions_to_copy.insert(pair<string,string>("SurfacePeriodic","SurfacePeriodic"));
        conditions_to_copy.insert(pair<string,string>("FluidStressCalc","FluxCalculation")); // a hack
        SCATRA::CreateScaTraDiscretization(fluiddis,scatradis,conditions_to_copy,false);
        if (comm.MyPID()==0)
        cout<<"Created scalar transport discretization from fluid field in...."
        <<time.ElapsedTime() << " secs\n\n";
      }
      else
        dserror("Fluid AND ScaTra discretization present. This is not supported.");

      // create an one-way coupling algorithm instance
      Teuchos::RCP<SCATRA::PassiveScaTraAlgorithm> algo = Teuchos::rcp(new SCATRA::PassiveScaTraAlgorithm(comm,scatradyn));

      if (restart)
      {
        // read the restart information, set vectors and variables
        algo->ReadRestart(restart);
      }

      // solve the whole (one-way-coupled) problem
      algo->TimeLoop();

      // summarize the performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform the result test
      DRT::ResultTestManager testmanager(comm);
      testmanager.AddFieldTest(algo->FluidField().CreateFieldTest());
      testmanager.AddFieldTest(algo->CreateScaTraFieldTest());
      testmanager.TestAll();

      break;
    } // case 2
    default:
      dserror("unknown velocity field type for transport of passive scalar");
  }

  return;

} // end of scatra_dyn()

#endif  // #ifdef CCADISCRET
