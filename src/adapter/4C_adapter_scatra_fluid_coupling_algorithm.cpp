// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_adapter_scatra_fluid_coupling_algorithm.hpp"

#include "4C_comm_mpi_utils.hpp"
#include "4C_coupling_adapter_volmortar.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fluid_turbulence_statistic_manager.hpp"
#include "4C_fluid_xfluid.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_levelset_algorithm.hpp"
#include "4C_xfem_discretization.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Adapter::ScaTraFluidCouplingAlgorithm::ScaTraFluidCouplingAlgorithm(Global::Problem& problem,
    MPI_Comm comm, const Teuchos::ParameterList& prbdyn, bool isale,
    const std::string scatra_disname, const Teuchos::ParameterList& solverparams)
    : AlgorithmBase(problem, comm, prbdyn),
      FluidBaseAlgorithm(problem, prbdyn, problem.fluid_dynamic_params(), "fluid", isale,
          false),  // false -> no immediate initialization of fluid time integration
      ScaTraBaseAlgorithm(problem, prbdyn, problem.scalar_transport_dynamic_params(), solverparams,
          scatra_disname, isale),
      fieldcoupling_(Teuchos::getIntegralValue<ScaTra::FieldCoupling>(
          problem.scalar_transport_dynamic_params(), "FIELDCOUPLING")),
      volcoupl_fluidscatra_(nullptr),
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
  auto& problem = AlgorithmBase::problem();

  check_is_init();

  // initialize scatra time integration scheme
  Adapter::ScaTraBaseAlgorithm::setup();

  // initialize fluid time integration scheme
  fluid_field()->init();

  // setup coupling adapter
  if (volcoupl_fluidscatra_)
    volcoupl_fluidscatra_->setup(problem.volmortar_params(), problem.cut_general_params());

  // set also initial field
  set_initial_flow_field(problem.fluid_dynamic_params());

  // transfer the initial convective velocity from initial fluid field to scalar transport field
  // subgrid scales not transferred since they are zero at time t=0.0
  if (!volcoupl_fluidscatra_)
  {
    scatra_field()->set_convective_velocity(*fluid_field()->convective_vel());
  }
  else
  {
    scatra_field()->set_convective_velocity(
        *volcoupl_fluidscatra_->apply_vector_mapping21(*fluid_field()->convective_vel()));
  }

  // ensure that both single field solvers use the same
  // time integration scheme
  switch (scatra_field()->method_name())
  {
    case ScaTra::timeint_stationary:
    {
      if (fluid_field()->tim_int_scheme() != FLUID::timeint_stationary)
        if (Core::Communication::my_mpi_rank(get_comm()) == 0)
          FOUR_C_THROW("Fluid and scatra time integration schemes do not match!");
      break;
    }
    case ScaTra::timeint_one_step_theta:
    {
      if (fluid_field()->tim_int_scheme() != FLUID::timeint_one_step_theta)
        if (Core::Communication::my_mpi_rank(get_comm()) == 0)
          std::cout << "WARNING: Fluid and scatra time integration schemes do not match!"
                    << std::endl;
      break;
    }
    case ScaTra::timeint_bdf2:
    {
      if (fluid_field()->tim_int_scheme() != FLUID::timeint_bdf2)
        if (Core::Communication::my_mpi_rank(get_comm()) == 0)
          std::cout << "WARNING: Fluid and scatra time integration schemes do not match!"
                    << std::endl;
      break;
    }
    case ScaTra::timeint_gen_alpha:
    {
      if (fluid_field()->tim_int_scheme() != FLUID::timeint_npgenalpha and
          fluid_field()->tim_int_scheme() != FLUID::timeint_afgenalpha)
        if (Core::Communication::my_mpi_rank(get_comm()) == 0)
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
  if (fluid_field()->turbulence_statistic_manager() != nullptr and
      scatra_field()->method_name() != ScaTra::timeint_stationary)
  {
    // Now, the statistics manager has access to the scatra time integration
    fluid_field()->turbulence_statistic_manager()->add_scatra_field(scatra_field());
  }

  // if available, allow scatra field to access dynamic Smagorinsky filter
  if (fluid_field()->dyn_smag_filter() != nullptr)
    scatra_field()->access_dyn_smag_filter(fluid_field()->dyn_smag_filter());

  // if available, allow scatra field to access dynamic Vreman
  if (fluid_field()->vreman() != nullptr) scatra_field()->access_vreman(fluid_field()->vreman());

  // safety check:
  if (volcoupl_fluidscatra_ == nullptr and fieldcoupling_ == ScaTra::coupling_volmortar)
    FOUR_C_THROW("Something went terrible wrong. Sorry about this!");

  set_is_setup(true);
  return;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidCouplingAlgorithm::setup_field_coupling(
    const std::string fluid_disname, const std::string scatra_disname)
{
  auto& problem = AlgorithmBase::problem();
  std::shared_ptr<Core::FE::Discretization> fluiddis = problem.get_dis(fluid_disname);
  std::shared_ptr<Core::FE::Discretization> scatradis = problem.get_dis(scatra_disname);

  if (fieldcoupling_ == ScaTra::coupling_volmortar)
  {
    // Scheme: non matching meshes --> volumetric mortar coupling...
    volcoupl_fluidscatra_ = std::make_shared<Coupling::Adapter::MortarVolCoupl>();

    // setup projection matrices (use default material strategy)
    volcoupl_fluidscatra_->init(
        problem.n_dim(), fluiddis, scatradis, nullptr, nullptr, nullptr, nullptr, nullptr, true);
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Adapter::ScaTraFluidCouplingAlgorithm::fluid_to_scatra(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> fluidvector) const
{
  switch (fieldcoupling_)
  {
    case ScaTra::coupling_match:
      return fluidvector;
      break;
    case ScaTra::coupling_volmortar:
      return volcoupl_fluidscatra_->apply_vector_mapping21(*fluidvector);
      break;
    default:
      FOUR_C_THROW("unknown field coupling type");
      return nullptr;
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
std::shared_ptr<const Core::LinAlg::Vector<double>>
Adapter::ScaTraFluidCouplingAlgorithm::scatra_to_fluid(
    const std::shared_ptr<const Core::LinAlg::Vector<double>> scatravector) const
{
  switch (fieldcoupling_)
  {
    case ScaTra::coupling_match:
      return scatravector;
      break;
    case ScaTra::coupling_volmortar:
      return volcoupl_fluidscatra_->apply_vector_mapping12(*scatravector);
      break;
    default:
      FOUR_C_THROW("unknown field coupling type");
      return nullptr;
      break;
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void Adapter::ScaTraFluidCouplingAlgorithm::read_restart(int step)
{
  auto& problem = AlgorithmBase::problem();

  fluid_field()->read_restart(step);
  scatra_field()->read_restart(step);
  set_time_step(fluid_field()->time(), step);

  // read scatra-specific restart data for turbulence statistics
  if (fluid_field()->turbulence_statistic_manager() != nullptr)
  {
    Core::IO::DiscretizationReader reader(
        *scatra_field()->discretization(), problem.input_control_file(), step);
    fluid_field()->turbulence_statistic_manager()->read_restart_scatra(reader, step);
  }
}

FOUR_C_NAMESPACE_CLOSE
