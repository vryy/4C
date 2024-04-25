/*----------------------------------------------------------------------*/
/*! \file
 \brief entry point (global control routine) for porous multiphase flow

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_porofluidmultiphase_dyn.hpp"

#include "4C_global_data.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_lib_discret.hpp"
#include "4C_lib_dofset_predefineddofnumber.hpp"
#include "4C_porofluidmultiphase_timint_implicit.hpp"
#include "4C_porofluidmultiphase_timint_ost.hpp"
#include "4C_porofluidmultiphase_utils.hpp"

#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*-------------------------------------------------------------------------------*
 | Main control routine for poro fluid multiphase problems           vuong 08/16 |
 *-------------------------------------------------------------------------------*/
void porofluidmultiphase_dyn(int restart)
{
  // define the discretization names
  const std::string fluid_disname = "porofluid";
  const std::string stuct_disname = "structure";
  const std::string artery_disname = "artery";

  // access the communicator
  const Epetra_Comm& comm = GLOBAL::Problem::Instance()->GetDis(fluid_disname)->Comm();

  // access the problem
  GLOBAL::Problem* problem = GLOBAL::Problem::Instance();

  // print problem type and logo
  if (comm.MyPID() == 0)
  {
    POROFLUIDMULTIPHASE::PrintLogo();
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << GLOBAL::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // Parameter reading
  // access structural dynamic params list which will be possibly modified while creating the time
  // integrator
  const Teuchos::ParameterList& porodyn = problem->PoroFluidMultiPhaseDynamicParams();

  // get the solver number used for poro fluid solver
  const int linsolvernumber = porodyn.get<int>("LINEAR_SOLVER");

  // -------------------------------------------------------------------
  // access the discretization(s)
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::Discretization> actdis = Teuchos::null;
  actdis = GLOBAL::Problem::Instance()->GetDis(fluid_disname);

  // possible interaction partners as seen from the artery elements
  // [artelegid; contelegid_1, ...contelegid_n]
  std::map<int, std::set<int>> nearbyelepairs;

  if (GLOBAL::Problem::Instance()->DoesExistDis(artery_disname))
  {
    Teuchos::RCP<DRT::Discretization> arterydis = Teuchos::null;
    arterydis = GLOBAL::Problem::Instance()->GetDis(artery_disname);
    // get the coupling method
    INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod arterycoupl =
        CORE::UTILS::IntegralValue<INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod>(
            porodyn.sublist("ARTERY COUPLING"), "ARTERY_COUPLING_METHOD");

    // lateral surface coupling active?
    const bool evaluate_on_lateral_surface = CORE::UTILS::IntegralValue<int>(
        porodyn.sublist("ARTERY COUPLING"), "LATERAL_SURFACE_COUPLING");

    switch (arterycoupl)
    {
      case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::gpts:
      case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::mp:
      case INPAR::ARTNET::ArteryPoroMultiphaseScatraCouplingMethod::ntp:
      {
        actdis->FillComplete();
        nearbyelepairs = POROFLUIDMULTIPHASE::UTILS::ExtendedGhostingArteryDiscretization(
            actdis, arterydis, evaluate_on_lateral_surface, arterycoupl);
        break;
      }
      default:
      {
        break;
      }
    }
  }

  // -------------------------------------------------------------------
  // assign dof set for solid pressures
  // -------------------------------------------------------------------
  Teuchos::RCP<DRT::DofSetInterface> dofsetaux =
      Teuchos::rcp(new DRT::DofSetPredefinedDoFNumber(1, 0, 0, false));
  const int nds_solidpressure = actdis->AddDofSet(dofsetaux);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  actdis->FillComplete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<IO::DiscretizationWriter> output = actdis->Writer();
  output->WriteMesh(0, 0.0);

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme timintscheme =
      CORE::UTILS::IntegralValue<INPAR::POROFLUIDMULTIPHASE::TimeIntegrationScheme>(
          porodyn, "TIMEINTEGR");

  // build poro fluid time integrator
  Teuchos::RCP<ADAPTER::PoroFluidMultiphase> algo = POROFLUIDMULTIPHASE::UTILS::CreateAlgorithm(
      timintscheme, actdis, linsolvernumber, porodyn, porodyn, output);

  // initialize
  algo->Init(false,       // eulerian formulation
      -1,                 //  no displacements
      -1,                 // no velocities
      nds_solidpressure,  // dof set for post processing solid pressure
      -1,                 // no scalar field
      &nearbyelepairs);   // possible interaction pairs

  // read the restart information, set vectors and variables
  if (restart) algo->ReadRestart(restart);

  // assign poro material for evaluation of porosity
  // note: to be done after potential restart, as in ReadRestart()
  //       the secondary material is destroyed
  POROFLUIDMULTIPHASE::UTILS::SetupMaterial(comm, stuct_disname, fluid_disname);

  // 4.- Run of the actual problem.
  algo->TimeLoop();

  // 4.3.- Summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test if required
  problem->AddFieldTest(algo->CreateFieldTest());
  problem->TestAll(comm);

  return;

}  // poromultiphase_dyn

FOUR_C_NAMESPACE_CLOSE
