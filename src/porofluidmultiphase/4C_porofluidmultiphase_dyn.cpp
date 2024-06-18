/*----------------------------------------------------------------------*/
/*! \file
 \brief entry point (global control routine) for porous multiphase flow

   \level 3

 *----------------------------------------------------------------------*/



#include "4C_porofluidmultiphase_dyn.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_bio.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
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
  const Epetra_Comm& comm = Global::Problem::Instance()->GetDis(fluid_disname)->Comm();

  // access the problem
  Global::Problem* problem = Global::Problem::Instance();

  // print problem type and logo
  if (comm.MyPID() == 0)
  {
    POROFLUIDMULTIPHASE::PrintLogo();
    std::cout << "###################################################" << std::endl;
    std::cout << "# YOUR PROBLEM TYPE: " << Global::Problem::Instance()->ProblemName() << std::endl;
    std::cout << "###################################################" << std::endl;
  }

  // Parameter reading
  // access structural dynamic params list which will be possibly modified while creating the time
  // integrator
  const Teuchos::ParameterList& porodyn = problem->poro_fluid_multi_phase_dynamic_params();

  // get the solver number used for poro fluid solver
  const int linsolvernumber = porodyn.get<int>("LINEAR_SOLVER");

  // -------------------------------------------------------------------
  // access the discretization(s)
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::FE::Discretization> actdis = Teuchos::null;
  actdis = Global::Problem::Instance()->GetDis(fluid_disname);

  // possible interaction partners as seen from the artery elements
  // [artelegid; contelegid_1, ...contelegid_n]
  std::map<int, std::set<int>> nearbyelepairs;

  if (Global::Problem::Instance()->DoesExistDis(artery_disname))
  {
    Teuchos::RCP<Core::FE::Discretization> arterydis = Teuchos::null;
    arterydis = Global::Problem::Instance()->GetDis(artery_disname);
    // get the coupling method
    Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod arterycoupl =
        Core::UTILS::IntegralValue<Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod>(
            porodyn.sublist("ARTERY COUPLING"), "ARTERY_COUPLING_METHOD");

    // lateral surface coupling active?
    const bool evaluate_on_lateral_surface = Core::UTILS::IntegralValue<int>(
        porodyn.sublist("ARTERY COUPLING"), "LATERAL_SURFACE_COUPLING");

    switch (arterycoupl)
    {
      case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::gpts:
      case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::mp:
      case Inpar::ArteryNetwork::ArteryPoroMultiphaseScatraCouplingMethod::ntp:
      {
        actdis->fill_complete();
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
  Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux =
      Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(1, 0, 0, false));
  const int nds_solidpressure = actdis->AddDofSet(dofsetaux);

  // -------------------------------------------------------------------
  // set degrees of freedom in the discretization
  // -------------------------------------------------------------------
  actdis->fill_complete();

  // -------------------------------------------------------------------
  // context for output and restart
  // -------------------------------------------------------------------
  Teuchos::RCP<Core::IO::DiscretizationWriter> output = actdis->Writer();
  output->write_mesh(0, 0.0);

  // -------------------------------------------------------------------
  // algorithm construction depending on
  // time-integration (or stationary) scheme
  // -------------------------------------------------------------------
  Inpar::POROFLUIDMULTIPHASE::TimeIntegrationScheme timintscheme =
      Core::UTILS::IntegralValue<Inpar::POROFLUIDMULTIPHASE::TimeIntegrationScheme>(
          porodyn, "TIMEINTEGR");

  // build poro fluid time integrator
  Teuchos::RCP<Adapter::PoroFluidMultiphase> algo = POROFLUIDMULTIPHASE::UTILS::CreateAlgorithm(
      timintscheme, actdis, linsolvernumber, porodyn, porodyn, output);

  // initialize
  algo->init(false,       // eulerian formulation
      -1,                 //  no displacements
      -1,                 // no velocities
      nds_solidpressure,  // dof set for post processing solid pressure
      -1,                 // no scalar field
      &nearbyelepairs);   // possible interaction pairs

  // read the restart information, set vectors and variables
  if (restart) algo->read_restart(restart);

  // assign poro material for evaluation of porosity
  // note: to be done after potential restart, as in read_restart()
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
