/*----------------------------------------------------------------------*/
/*!
 \file porofluidmultiphase_dyn.cpp

 \brief entry point (global control routine) for porous multiphase flow

   \level 3

   \maintainer  Anh-Tu Vuong
                vuong@lnm.mw.tum.de
                http://www.lnm.mw.tum.de
                089 - 289-15251
 *----------------------------------------------------------------------*/



#include "porofluidmultiphase_dyn.H"

#include "porofluidmultiphase_utils.H"

#include "porofluidmultiphase_timint_implicit.H"
#include "porofluidmultiphase_timint_ost.H"

#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_discret.H"

#include "../drt_io/io.H"
#include "../drt_io/io_control.H"


#include <Teuchos_TimeMonitor.hpp>

/*-------------------------------------------------------------------------------*
 | Main control routine for poro fluid multiphase problems           vuong 08/16 |
 *-------------------------------------------------------------------------------*/
void porofluidmultiphase_dyn(int restart)
{
  // define the discretization names
  const std::string fluid_disname = "porofluid";
  const std::string stuct_disname = "structure";

  // access the communicator
  const Epetra_Comm& comm = DRT::Problem::Instance()->GetDis(fluid_disname)->Comm();

  // access the problem
  DRT::Problem* problem = DRT::Problem::Instance();

  // print problem type
  if (comm.MyPID() == 0)
  {
    std::cout << "###################################################"
        << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: "
        << DRT::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################"
        << std::endl;
  }

  //Parameter reading
  // access structural dynamic params list which will be possibly modified while creating the time integrator
  const Teuchos::ParameterList& porodyn  = problem->PoroFluidMultiPhaseDynamicParams();

  // get the solver number used for poro fluid solver
  const int linsolvernumber = porodyn.get<int>("LINEAR_SOLVER");

  // -------------------------------------------------------------------
  // access the discretization
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = DRT::Problem::Instance()->GetDis(fluid_disname);

  // -------------------------------------------------------------------
  // assign dof set for solid pressures
  // -------------------------------------------------------------------
  const int nds_solidpressure = actdis->BuildDofSetAuxProxy(1,0,0,false);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0,0.0);

  // Creation of Coupled Problem algortihm.
  Teuchos::RCP<POROFLUIDMULTIPHASE::TimIntImpl> algo = Teuchos::null;

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme timintscheme =
    DRT::INPUT::IntegralValue<INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme>(porodyn,"TIMEINTEGR");

  switch(timintscheme)
  {
  case INPAR::POROFLUIDMULTIPHASE::timeint_one_step_theta:
  {
    // create algorithm
    algo = Teuchos::rcp(new POROFLUIDMULTIPHASE::TimIntOneStepTheta(
            actdis,
            linsolvernumber,
            porodyn,
            porodyn,
            DRT::Problem::Instance()->ErrorFile()->Handle(),
            output));
    break;
  }
  default:
    dserror("Unknown time-integration scheme for multiphase poro fluid problem");
    break;
  }

  // initialize
  algo->Init(
      false, // eulerian formulation
      -1,    //  no displacements
      -1,    // no velocities
      nds_solidpressure, // dof set for post processing solid pressure
      -1   // no scalar field
      );

  // read the restart information, set vectors and variables
  if (restart)
    algo->ReadRestart(restart);

  // assign poro material for evaluation of porosity
  // note: to be done after potential restart, as in ReadRestart()
  //       the secondary material is destroyed
  POROFLUIDMULTIPHASE::UTILS::SetupMaterial(comm,stuct_disname,fluid_disname);

  //4.- Run of the actual problem.
  algo->TimeLoop();

  // 4.3.- Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test if required
  problem->AddFieldTest(algo->CreateFieldTest());
  problem->TestAll(comm);

  return;

} // poromultiphase_dyn
