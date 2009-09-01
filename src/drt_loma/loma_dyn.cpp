/*!----------------------------------------------------------------------
\file loma_dyn.H
\brief Control routine for low-Mach-number flow module.


<pre>
Maintainer: Volker Gravemeier
            vgravem@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089/28915245
</pre>

*----------------------------------------------------------------------*/

#ifdef CCADISCRET

#include <string>
#include <iostream>

#ifdef PARALLEL
#include <mpi.h>
#include <Epetra_MpiComm.h>
#else
#include <Epetra_SerialComm.h>
#endif

#include "loma_dyn.H"
#include "loma_algorithm.H"
#include "../drt_scatra/scatra_utils.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*/
// entry point for LOMA in DRT
/*----------------------------------------------------------------------*/
void loma_dyn(int disnumff,int disnumscatra, int restart)
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // print warning to screen
  if (comm.MyPID()==0) cout<<"You are now about to enter the module for low-Mach-number flow!"<<endl;

  // access the fluid discretization
  RefCountPtr<DRT::Discretization> fluiddis = DRT::Problem::Instance()->Dis(disnumff,0);
  if (!fluiddis->Filled()) fluiddis->FillComplete();

  // access the (typically empty) scatra discretization
  RefCountPtr<DRT::Discretization> scatradis = DRT::Problem::Instance()->Dis(disnumscatra,0);
  if (!scatradis->Filled()) scatradis->FillComplete();

  // access the problem-specific parameter list
  const Teuchos::ParameterList& lomacontrol = DRT::Problem::Instance()->LOMAControlParams();

  // access the scalar transport parameter list
  const Teuchos::ParameterList& scatradyn = DRT::Problem::Instance()->ScalarTransportDynamicParams();
  int veltype = Teuchos::getIntegralValue<int>(scatradyn,"VELOCITYFIELD");

  // choose algorithm depending on velocity field type
  switch (veltype)
  {
  case 0:  // zero  (see case 1)
  case 1:  // function
  {
    // we directly use the elements from the scalar transport elements section
    if (scatradis->NumGlobalNodes()==0)
      dserror("No elements in the ---TRANSPORT ELEMENTS section");

    // create instance of scalar transport basis algorithm (empty fluid discretization)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = rcp(new ADAPTER::ScaTraBaseAlgorithm(lomacontrol,false));

    // read the restart information, set vectors and variables
    if (restart) scatraonly->ScaTraField().ReadRestart(restart);

    // set velocity field
    //(this is done only once. Time-dependent velocity fields are not supported)
    (scatraonly->ScaTraField()).SetVelocityField();

    // enter time loop to solve problem with given convective velocity
    (scatraonly->ScaTraField()).TimeLoop();

    // perform the result test if required
    DRT::Problem::Instance()->AddFieldTest(scatraonly->CreateScaTraFieldTest());
    DRT::Problem::Instance()->TestAll(comm);

    break;
  }
  case 2:  // Navier_Stokes
  {
    // we use the fluid discretization as layout for the scalar transport discretization
    if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

    // create scatra elements if the scatra discretization is empty (typical case)
    if (scatradis->NumGlobalNodes()==0)
    {
      Epetra_Time time(comm);
      std::map<string,string> conditions_to_copy;
      conditions_to_copy.insert(pair<string,string>("TransportDirichlet","Dirichlet"));
      conditions_to_copy.insert(pair<string,string>("TransportPointNeumann","PointNeumann"));
      conditions_to_copy.insert(pair<string,string>("TransportLineNeumann","LineNeumann"));
      conditions_to_copy.insert(pair<string,string>("TransportSurfaceNeumann","SurfaceNeumann"));
      conditions_to_copy.insert(pair<string,string>("TransportVolumeNeumann","VolumeNeumann"));
      conditions_to_copy.insert(pair<string,string>("LinePeriodic","LinePeriodic"));
      conditions_to_copy.insert(pair<string,string>("SurfacePeriodic","SurfacePeriodic"));
      conditions_to_copy.insert(pair<string,string>("TransportNeumannInflow","TransportNeumannInflow"));
      conditions_to_copy.insert(pair<string,string>("FluidStressCalc","FluxCalculation"));

      // fetch the desired material id for the transport elements
      const int matid = scatradyn.get<int>("MATID");

      // create the scatra discretization
      SCATRA::CreateScaTraDiscretization(fluiddis,scatradis,conditions_to_copy,matid,false);
      if (comm.MyPID()==0)
        cout<<"Created scalar transport discretization from fluid field in...."
        <<time.ElapsedTime() << " secs\n\n";
    }
    else dserror("Fluid AND Scatra discretization present. This is not supported.");

    // create a LOMA::Algorithm instance
    Teuchos::RCP<LOMA::Algorithm> loma = Teuchos::rcp(new LOMA::Algorithm(comm,lomacontrol));

    // read the restart information, set vectors and variables
    if (restart) loma->ReadRestart(restart);

    // enter LOMA algorithm
    loma->TimeLoop();

    // summarize the performance measurements
    Teuchos::TimeMonitor::summarize();

    // perform the result test
    DRT::Problem::Instance()->AddFieldTest(loma->FluidField().CreateFieldTest());
    DRT::Problem::Instance()->AddFieldTest(loma->CreateScaTraFieldTest());
    DRT::Problem::Instance()->TestAll(comm);

    break;
  } // case 2
  default:
    dserror("Unknown velocity field type for low-Mach-number flow: %d",veltype);
  }

  return;

} // loma_dyn()


#endif  // CCADISCRET
