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
#include "../drt_inpar/drt_validparameters.H"
#include "../drt_scatra/scatra_utils.H"
#include <Teuchos_TimeMonitor.hpp>
#include <Teuchos_Time.hpp>
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_utils_createdis.H"
#include <Epetra_Time.h>


/*----------------------------------------------------------------------*/
// entry point for LOMA in DRT
/*----------------------------------------------------------------------*/
void loma_dyn(int disnumff, int disnumscatra, int restart)
{
  // create a communicator
#ifdef PARALLEL
  Epetra_MpiComm comm(MPI_COMM_WORLD);
#else
  Epetra_SerialComm comm;
#endif

  // print warning to screen
  if (comm.MyPID()==0)
    cout << "You are now about to enter the module for low-Mach-number flow!" <<endl;

  // define abbreviation
  RCP<DRT::Problem> problem = DRT::Problem::Instance();

  // access fluid and (typically empty) scatra discretization
  RefCountPtr<DRT::Discretization> fluiddis  = problem->Dis(disnumff,0);
  RefCountPtr<DRT::Discretization> scatradis = problem->Dis(disnumscatra,0);

  // ensure that all dofs are assigned in the right order such that
  // dof numbers are created with fluid dof < scatra/elch dof
  fluiddis->FillComplete();
  scatradis->FillComplete();

  // access problem-specific parameter list for LOMA
  const Teuchos::ParameterList& lomacontrol = problem->LOMAControlParams();

  // access parameter list for scatra
  const Teuchos::ParameterList& scatradyn = problem->ScalarTransportDynamicParams();

  // identify type of velocity field
  const INPAR::SCATRA::VelocityField veltype = DRT::INPUT::IntegralValue<INPAR::SCATRA::VelocityField>(scatradyn,"VELOCITYFIELD");

  // choose algorithm depending on type of velocity field
  switch (veltype)
  {
  case INPAR::SCATRA::velocity_zero:  // zero velocity field (see case 1)
  case INPAR::SCATRA::velocity_function:  // velocity field prescribed by function
  {
    // directly use elements from input section 'transport elements'
    if (scatradis->NumGlobalNodes()==0)
      dserror("No elements in input section ---TRANSPORT ELEMENTS!");

    // create instance of scalar transport basis algorithm (no fluid discretization)
    Teuchos::RCP<ADAPTER::ScaTraBaseAlgorithm> scatraonly = rcp(new ADAPTER::ScaTraBaseAlgorithm(lomacontrol,false));

    // read restart information
    if (restart) scatraonly->ScaTraField().ReadRestart(restart);

    // set velocity field (done only once, time-dependent fields not supported)
    (scatraonly->ScaTraField()).SetVelocityField();

    // enter time loop to solve problem with given convective velocity field
    (scatraonly->ScaTraField()).TimeLoop();

    // perform result test if required
    problem->AddFieldTest(scatraonly->CreateScaTraFieldTest());
    problem->TestAll(comm);

    break;
  }
  case INPAR::SCATRA::velocity_Navier_Stokes:  // Navier_Stokes
  {
    // use fluid discretization as layout for scatra discretization
    if (fluiddis->NumGlobalNodes()==0) dserror("Fluid discretization is empty!");

    // create scatra elements if scatra discretization is empty (typical case)
    if (scatradis->NumGlobalNodes()==0)
    {
      Epetra_Time time(comm);

      // get desired material id for scatra elements
      const int matid = SCATRA::GetScaTraMatID(scatradyn);

      // create scatra discretization
      {
        Teuchos::RCP<DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy> > clonewizard = Teuchos::rcp(new DRT::UTILS::DiscretizationCreator<SCATRA::ScatraFluidCloneStrategy>() );

        clonewizard->CreateMatchingDiscretization(fluiddis,scatradis,matid);
      }
      if (comm.MyPID()==0)
        cout << "Created scalar transport discretization from fluid discretization in ... "
        << time.ElapsedTime() << " secs\n\n";
    }
    else dserror("Fluid AND ScaTra discretization present. This is not supported.");

    // create a LOMA::Algorithm instance
    Teuchos::RCP<LOMA::Algorithm> loma = Teuchos::rcp(new LOMA::Algorithm(comm,lomacontrol));

    // read restart information
    if (restart) loma->ReadRestart(restart);

    // enter LOMA algorithm
    loma->TimeLoop();

    // summarize performance measurements
    Teuchos::TimeMonitor::summarize();

    // perform result test if required
    problem->AddFieldTest(loma->FluidField().CreateFieldTest());
    problem->AddFieldTest(loma->CreateScaTraFieldTest());
    problem->TestAll(comm);

    break;
  }
  default:
    dserror("Unknown velocity field type for low-Mach-number flow: %d",veltype);
  }

  return;

} // loma_dyn()


#endif  // CCADISCRET
