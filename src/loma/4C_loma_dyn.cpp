/*---------------------------------------------------------------------*/
/*! \file

\brief Control routine for low-Mach-number flow module.

\level 2


*/
/*---------------------------------------------------------------------*/


#include "4C_loma_dyn.hpp"

#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_validparameters.hpp"
#include "4C_loma_algorithm.hpp"
#include "4C_scatra_ele.hpp"
#include "4C_scatra_timint_implicit.hpp"
#include "4C_scatra_utils_clonestrategy.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <iostream>
#include <string>

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*/
// entry point for LOMA in DRT
/*----------------------------------------------------------------------*/
void loma_dyn(int restart)
{
  // create a communicator
  const Epetra_Comm& comm = Global::Problem::instance()->get_dis("fluid")->get_comm();

  // print warning to screen
  if (comm.MyPID() == 0)
    std::cout << "You are now about to enter the module for low-Mach-number flow!" << std::endl;

  // define abbreviation
  Global::Problem* problem = Global::Problem::instance();

  // access fluid and (typically empty) scatra discretization
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->get_dis("fluid");
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->get_dis("scatra");

  // ensure that all dofs are assigned in the right order such that
  // dof numbers are created with fluid dof < scatra/elch dof
  fluiddis->fill_complete();
  scatradis->fill_complete();

  // access problem-specific parameter list for LOMA
  const Teuchos::ParameterList& lomacontrol = problem->loma_control_params();

  // access parameter list for scatra
  const Teuchos::ParameterList& scatradyn = problem->scalar_transport_dynamic_params();

  // access parameter list for fluid
  const Teuchos::ParameterList& fdyn = problem->fluid_dynamic_params();

  // identify type of velocity field
  const Inpar::ScaTra::VelocityField veltype =
      Core::UTILS::IntegralValue<Inpar::ScaTra::VelocityField>(scatradyn, "VELOCITYFIELD");

  // choose algorithm depending on type of velocity field
  switch (veltype)
  {
    case Inpar::ScaTra::velocity_zero:      // zero velocity field (see case 1)
    case Inpar::ScaTra::velocity_function:  // velocity field prescribed by function
    {
      // directly use elements from input section 'transport elements'
      if (scatradis->num_global_nodes() == 0)
        FOUR_C_THROW("No elements in input section ---TRANSPORT ELEMENTS!");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        FOUR_C_THROW(
            "no linear solver defined for LOMA problem. Please set LINEAR_SOLVER in SCALAR "
            "TRANSPORT DYNAMIC to a valid number!");

      // create instance of scalar transport basis algorithm (no fluid discretization)
      Teuchos::RCP<Adapter::ScaTraBaseAlgorithm> scatraonly =
          Teuchos::rcp(new Adapter::ScaTraBaseAlgorithm(
              lomacontrol, scatradyn, Global::Problem::instance()->solver_params(linsolvernumber)));

      // add proxy of velocity related degrees of freedom to scatra discretization
      Teuchos::RCP<Core::DOFSets::DofSetInterface> dofsetaux =
          Teuchos::rcp(new Core::DOFSets::DofSetPredefinedDoFNumber(
              Global::Problem::instance()->n_dim() + 1, 0, 0, true));
      if (scatradis->add_dof_set(dofsetaux) != 1)
        FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
      scatraonly->sca_tra_field()->set_number_of_dof_set_velocity(1);

      // now we can call init() on base algo
      scatraonly->init();

      // only now we must call setup() on the scatra time integrator.
      // all objects relying on the parallel distribution are
      // created and pointers are set.
      // calls setup() on the time integrator inside
      scatraonly->setup();

      // read restart information
      if (restart) (scatraonly->sca_tra_field())->read_restart(restart);

      // set initial velocity field
      // note: The order read_restart() before set_velocity_field() is important here!!
      // for time-dependent velocity fields, set_velocity_field() is additionally called in each
      // prepare_time_step()-call
      (scatraonly->sca_tra_field())->set_velocity_field();

      // enter time loop to solve problem with given convective velocity field
      (scatraonly->sca_tra_field())->time_loop();

      // perform result test if required
      problem->add_field_test(scatraonly->create_sca_tra_field_test());
      problem->test_all(comm);

      break;
    }
    case Inpar::ScaTra::velocity_Navier_Stokes:  // Navier_Stokes
    {
      // use fluid discretization as layout for scatra discretization
      if (fluiddis->num_global_nodes() == 0) FOUR_C_THROW("Fluid discretization is empty!");

      // to generate turbulent flow in the inflow section only, it is not necessary to
      // solve the transport equation for the temperature
      // therefore, use problem type fluid
      if ((Core::UTILS::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") ==
              true) and
          (restart < fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
        FOUR_C_THROW("Choose problem type fluid to generate turbulent flow in the inflow section!");

      // create scatra elements if scatra discretization is empty (typical case)
      if (scatradis->num_global_nodes() == 0)
      {
        // fill scatra discretization by cloning fluid discretization
        Core::FE::CloneDiscretization<ScaTra::ScatraFluidCloneStrategy>(
            fluiddis, scatradis, Global::Problem::instance()->cloning_material_map());

        // set implementation type of cloned scatra elements to loma
        for (int i = 0; i < scatradis->num_my_col_elements(); ++i)
        {
          Discret::ELEMENTS::Transport* element =
              dynamic_cast<Discret::ELEMENTS::Transport*>(scatradis->l_col_element(i));
          if (element == nullptr)
            FOUR_C_THROW("Invalid element type!");
          else
            element->set_impl_type(Inpar::ScaTra::impltype_loma);
        }
      }
      else
        FOUR_C_THROW("Fluid AND ScaTra discretization present. This is not supported.");

      // get linear solver id from SCALAR TRANSPORT DYNAMIC
      const int linsolvernumber = scatradyn.get<int>("LINEAR_SOLVER");
      if (linsolvernumber == (-1))
        FOUR_C_THROW(
            "no linear solver defined for LOMA problem. Please set LINEAR_SOLVER in SCALAR "
            "TRANSPORT DYNAMIC to a valid number!");

      // create a LowMach::Algorithm instance
      Teuchos::RCP<LowMach::Algorithm> loma = Teuchos::rcp(new LowMach::Algorithm(
          comm, lomacontrol, Global::Problem::instance()->solver_params(linsolvernumber)));

      // add proxy of fluid transport degrees of freedom to scatra discretization
      if (scatradis->add_dof_set(fluiddis->get_dof_set_proxy()) != 1)
        FOUR_C_THROW("Scatra discretization has illegal number of dofsets!");
      loma->sca_tra_field()->set_number_of_dof_set_velocity(1);

      loma->init();

      loma->setup();

      // read restart information
      // in case a inflow generation in the inflow section has been performed, there are not any
      // scatra results available and the initial field is used
      if (restart)
      {
        if ((Core::UTILS::IntegralValue<int>(fdyn.sublist("TURBULENT INFLOW"), "TURBULENTINFLOW") ==
                true) and
            (restart == fdyn.sublist("TURBULENT INFLOW").get<int>("NUMINFLOWSTEP")))
          loma->read_inflow_restart(restart);
        else
          loma->read_restart(restart);
      }

      // enter LOMA algorithm
      loma->time_loop();

      // summarize performance measurements
      Teuchos::TimeMonitor::summarize();

      // perform result test if required
      problem->add_field_test(loma->fluid_field()->create_field_test());
      problem->add_field_test(loma->create_sca_tra_field_test());
      problem->test_all(comm);

      break;
    }
    default:
      FOUR_C_THROW("Unknown velocity field type for low-Mach-number flow: %d", veltype);
      break;
  }

  return;

}  // loma_dyn()

FOUR_C_NAMESPACE_CLOSE
