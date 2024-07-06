/*----------------------------------------------------------------------*/
/*! \file

\brief Basis of all algorithms that perform a coupling between Navier-Stokes
       and (active or passive) scalar transport equations

\level 1


*/
/*----------------------------------------------------------------------*/


#include "4C_adapter_scatra_fluid_coupling_algorithm.hpp"

#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_turbulence_statistic_manager.hpp"
#include "4C_fluid_xfluid.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_levelset_algorithm.hpp"
#include "4C_xfem_discretization.hpp"

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::ScaTraFluidCouplingAlgorithm::ScaTraFluidCouplingAlgorithm(const Epetra_Comm& comm,
    const Teuchos::ParameterList& prbdyn, bool isale, const std::string scatra_disname,
    const Teuchos::ParameterList& solverparams)
    : AlgorithmBase(comm, prbdyn),
      FluidBaseAlgorithm(prbdyn, Global::Problem::instance()->fluid_dynamic_params(), "fluid",
          isale,
          false),  // false -> no immediate initialization of fluid time integration
      ScaTraBaseAlgorithm(prbdyn, Global::Problem::instance()->scalar_transport_dynamic_params(),
          solverparams, scatra_disname, isale),
      fieldcoupling_(Core::UTILS::IntegralValue<Inpar::ScaTra::FieldCoupling>(
          Global::Problem::instance()->scalar_transport_dynamic_params(), "FIELDCOUPLING")),
      volcoupl_fluidscatra_(Teuchos::null),
      params_(prbdyn),
      scatra_disname_(scatra_disname),
      issetup_(false),
      isinit_(false)
{
  // keep constructor empty
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidCouplingAlgorithm::init()
{
  set_is_setup(false);

  Adapter::ScaTraBaseAlgorithm::init();

  // perform algorithm specific initialization stuff
  do_algorithm_specific_init();

  // do potential volmortar business
  setup_field_coupling("fluid", scatra_disname_);

  set_is_init(true);
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidCouplingAlgorithm::setup()
{
  check_is_init();

  // initialize scatra time integration scheme
  Adapter::ScaTraBaseAlgorithm::setup();

  // initialize fluid time integration scheme
  fluid_field()->init();

  // setup coupling adapter
  if (not volcoupl_fluidscatra_.is_null())
    volcoupl_fluidscatra_->setup(Global::Problem::instance()->volmortar_params(),
        Global::Problem::instance()->cut_general_params());

  // set also initial field
  set_initial_flow_field(Global::Problem::instance()->fluid_dynamic_params());

  if (Global::Problem::instance()->get_problem_type() != Core::ProblemType::fluid_xfem_ls)
  {
    // transfer the initial convective velocity from initial fluid field to scalar transport field
    // subgrid scales not transferred since they are zero at time t=0.0
    if (volcoupl_fluidscatra_.is_null())
      sca_tra_field()->set_velocity_field(
          fluid_field()->convective_vel(), Teuchos::null, Teuchos::null, Teuchos::null);
    else
      sca_tra_field()->set_velocity_field(
          volcoupl_fluidscatra_->apply_vector_mapping21(fluid_field()->convective_vel()),
          Teuchos::null, Teuchos::null, Teuchos::null);
  }

  // ensure that both single field solvers use the same
  // time integration scheme
  switch (sca_tra_field()->method_name())
  {
    case Inpar::ScaTra::timeint_stationary:
    {
      if (fluid_field()->tim_int_scheme() != Inpar::FLUID::timeint_stationary)
        if (get_comm().MyPID() == 0)
          FOUR_C_THROW("Fluid and scatra time integration schemes do not match!");
      break;
    }
    case Inpar::ScaTra::timeint_one_step_theta:
    {
      if (fluid_field()->tim_int_scheme() != Inpar::FLUID::timeint_one_step_theta)
        if (get_comm().MyPID() == 0)
          std::cout << "WARNING: Fluid and scatra time integration schemes do not match!"
                    << std::endl;
      break;
    }
    case Inpar::ScaTra::timeint_bdf2:
    {
      if (fluid_field()->tim_int_scheme() != Inpar::FLUID::timeint_bdf2)
        if (get_comm().MyPID() == 0)
          std::cout << "WARNING: Fluid and scatra time integration schemes do not match!"
                    << std::endl;
      break;
    }
    case Inpar::ScaTra::timeint_gen_alpha:
    {
      if (fluid_field()->tim_int_scheme() != Inpar::FLUID::timeint_npgenalpha and
          fluid_field()->tim_int_scheme() != Inpar::FLUID::timeint_afgenalpha)
        if (get_comm().MyPID() == 0)
          std::cout << "WARNING: Fluid and scatra time integration schemes do not match!"
                    << std::endl;
      break;
    }
    default:
    {
      FOUR_C_THROW("Time integration scheme for scalar transport not recognized!");
      break;
    }
  }

  // if applicable, provide scatra data to the turbulence statistics
  if (fluid_field()->turbulence_statistic_manager() != Teuchos::null and
      sca_tra_field()->method_name() != Inpar::ScaTra::timeint_stationary)
  {
    // Now, the statistics manager has access to the scatra time integration
    fluid_field()->turbulence_statistic_manager()->add_sca_tra_field(sca_tra_field());
  }

  // if available, allow scatra field to access dynamic Smagorinsky filter
  if (fluid_field()->dyn_smag_filter() != Teuchos::null)
    sca_tra_field()->access_dyn_smag_filter(fluid_field()->dyn_smag_filter());

  // if available, allow scatra field to access dynamic Vreman
  if (fluid_field()->vreman() != Teuchos::null)
    sca_tra_field()->access_vreman(fluid_field()->vreman());

  // safety check:
  if (volcoupl_fluidscatra_ == Teuchos::null and
      fieldcoupling_ == Inpar::ScaTra::coupling_volmortar)
    FOUR_C_THROW("Something went terrible wrong. Sorry about this!");

  set_is_setup(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidCouplingAlgorithm::setup_field_coupling(
    const std::string fluid_disname, const std::string scatra_disname)
{
  Global::Problem* problem = Global::Problem::instance();
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->get_dis(fluid_disname);
  Teuchos::RCP<Core::FE::Discretization> scatradis = problem->get_dis(scatra_disname);

  if (fieldcoupling_ == Inpar::ScaTra::coupling_volmortar)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_fluidscatra_ = Teuchos::rcp(new Core::Adapter::MortarVolCoupl());

    // setup projection matrices (use default material strategy)
    volcoupl_fluidscatra_->init(problem->n_dim(), fluiddis, scatradis, nullptr, nullptr, nullptr,
        nullptr, Teuchos::null, true);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::ScaTraFluidCouplingAlgorithm::fluid_to_scatra(
    const Teuchos::RCP<const Epetra_Vector> fluidvector) const
{
  switch (fieldcoupling_)
  {
    case Inpar::ScaTra::coupling_match:
      return fluidvector;
      break;
    case Inpar::ScaTra::coupling_volmortar:
      return volcoupl_fluidscatra_->apply_vector_mapping21(fluidvector);
      break;
    default:
      FOUR_C_THROW("unknown field coupling type");
      return Teuchos::null;
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<const Epetra_Vector> Adapter::ScaTraFluidCouplingAlgorithm::scatra_to_fluid(
    const Teuchos::RCP<const Epetra_Vector> scatravector) const
{
  switch (fieldcoupling_)
  {
    case Inpar::ScaTra::coupling_match:
      return scatravector;
      break;
    case Inpar::ScaTra::coupling_volmortar:
      return volcoupl_fluidscatra_->apply_vector_mapping12(scatravector);
      break;
    default:
      FOUR_C_THROW("unknown field coupling type");
      return Teuchos::null;
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidCouplingAlgorithm::read_restart(int step)
{
  fluid_field()->read_restart(step);
  sca_tra_field()->read_restart(step);
  set_time_step(fluid_field()->time(), step);

  // read scatra-specific restart data for turbulence statistics
  if (fluid_field()->turbulence_statistic_manager() != Teuchos::null)
  {
    Core::IO::DiscretizationReader reader(
        sca_tra_field()->discretization(), Global::Problem::instance()->input_control_file(), step);
    fluid_field()->turbulence_statistic_manager()->read_restart_sca_tra(reader, step);
  }
}

FOUR_C_NAMESPACE_CLOSE
