// This file is part of 4C multiphysics licensed under the
// GNU Lesser General Public License v3.0 or later.
//
// See the LICENSE.md file in the top-level for license information.
//
// SPDX-License-Identifier: LGPL-3.0-or-later

#include "4C_porofluid_pressure_based_dyn.hpp"

#include "4C_art_net_input.hpp"
#include "4C_fem_discretization.hpp"
#include "4C_fem_dofset_predefineddofnumber.hpp"
#include "4C_global_data.hpp"
#include "4C_io.hpp"
#include "4C_io_control.hpp"
#include "4C_porofluid_pressure_based_algorithm.hpp"
#include "4C_porofluid_pressure_based_algorithm_dependencies.hpp"
#include "4C_porofluid_pressure_based_utils.hpp"

#include <Teuchos_StandardParameterEntryValidators.hpp>

#include <functional>

FOUR_C_NAMESPACE_OPEN

namespace
{
  struct PorofluidPressureBasedMainContext
  {
    std::shared_ptr<Core::FE::Discretization> fluid_discretization;
    std::shared_ptr<Core::FE::Discretization> artery_discretization;
    bool has_artery_discretization;
    std::map<std::pair<std::string, std::string>, std::map<int, int>> cloning_material_map;
    PoroPressureBased::PorofluidAlgorithmDeps algorithm_deps;
    MPI_Comm comm;
    std::string problem_name;
    const Teuchos::ParameterList& porofluid_dynamic_parameters;
    std::function<void(MPI_Comm)> run_result_tests;
  };

  PorofluidPressureBasedMainContext make_main_context_from_problem(Global::Problem& problem)
  {
    const std::string fluid_disname = "porofluid";
    const std::string artery_disname = "artery";

    std::shared_ptr<Core::FE::Discretization> fluid_discretization = problem.get_dis(fluid_disname);
    const bool has_artery_discretization = problem.does_exist_dis(artery_disname);
    const std::shared_ptr<Core::FE::Discretization> artery_discretization =
        has_artery_discretization ? problem.get_dis(artery_disname) : nullptr;

    return PorofluidPressureBasedMainContext{
        .fluid_discretization = fluid_discretization,
        .artery_discretization = artery_discretization,
        .has_artery_discretization = has_artery_discretization,
        .cloning_material_map = problem.cloning_material_map(),
        .algorithm_deps = PoroPressureBased::make_algorithm_deps_from_problem(problem),
        .comm = fluid_discretization->get_comm(),
        .problem_name = problem.problem_name(),
        .porofluid_dynamic_parameters = problem.porofluid_pressure_based_dynamic_params(),
        .run_result_tests = [&problem](MPI_Comm comm) { problem.test_all(comm); },
    };
  }

  void run_porofluid_pressure_based(
      const PorofluidPressureBasedMainContext& context, const int restart)
  {
    // define the discretization names
    const std::string fluid_disname = "porofluid";
    const std::string struct_disname = "structure";
    const std::string artery_disname = "artery";

    // print problem type
    if (Core::Communication::my_mpi_rank(context.comm) == 0)
    {
      std::cout << "###################################################" << std::endl;
      std::cout << "# YOUR PROBLEM TYPE: " << context.problem_name << std::endl;
      std::cout << "###################################################" << std::endl;
    }

    // Parameter reading
    // access structural dynamic params list which will be possibly modified while creating the
    // time integrator
    const Teuchos::ParameterList& porodyn = context.porofluid_dynamic_parameters;

    // get the solver number used for poro fluid solver
    const int linsolvernumber = porodyn.sublist("nonlinear_solver").get<int>("linear_solver_id");

    // -------------------------------------------------------------------
    // access the discretization(s)
    // -------------------------------------------------------------------
    std::shared_ptr<Core::FE::Discretization> actdis = context.fluid_discretization;

    // possible interaction partners as seen from the artery elements
    // [artelegid; contelegid_1, ...contelegid_n]
    std::map<int, std::set<int>> nearby_ele_pairs;

    if (context.has_artery_discretization)
    {
      const std::shared_ptr<Core::FE::Discretization> arterydis = context.artery_discretization;
      // get the coupling method
      auto arterycoupl =
          Teuchos::getIntegralValue<ArteryNetwork::ArteryPorofluidElastScatraCouplingMethod>(
              porodyn.sublist("artery_coupling"), "coupling_method");

      // lateral surface coupling active?
      const bool evaluate_on_lateral_surface =
          porodyn.sublist("artery_coupling").get<bool>("lateral_surface_coupling");

      switch (arterycoupl)
      {
        case ArteryNetwork::ArteryPorofluidElastScatraCouplingMethod::gauss_point_to_segment:
        case ArteryNetwork::ArteryPorofluidElastScatraCouplingMethod::mortar_penalty:
        case ArteryNetwork::ArteryPorofluidElastScatraCouplingMethod::node_to_point:
        {
          actdis->fill_complete();
          nearby_ele_pairs = PoroPressureBased::extended_ghosting_artery_discretization(
              *actdis, arterydis, evaluate_on_lateral_surface, arterycoupl);
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
    std::shared_ptr<Core::DOFSets::DofSetInterface> dofsetaux =
        std::make_shared<Core::DOFSets::DofSetPredefinedDoFNumber>(1, 0, 0, false);
    const int nds_solidpressure = actdis->add_dof_set(dofsetaux);

    // -------------------------------------------------------------------
    // set degrees of freedom in the discretization
    // -------------------------------------------------------------------
    actdis->fill_complete();

    // -------------------------------------------------------------------
    // context for output and restart
    // -------------------------------------------------------------------
    std::shared_ptr<Core::IO::DiscretizationWriter> output = actdis->writer();
    output->write_mesh(0, 0.0);

    // -------------------------------------------------------------------
    // algorithm construction depending on
    // time-integration (or stationary) scheme
    // -------------------------------------------------------------------
    auto timintscheme =
        porodyn.sublist("time_integration").get<PoroPressureBased::TimeIntegrationScheme>("scheme");

    // build poro fluid time integrator
    std::shared_ptr<Adapter::PoroFluidMultiphase> algo = PoroPressureBased::create_algorithm(
        timintscheme, actdis, linsolvernumber, porodyn, porodyn, output, context.algorithm_deps);

    // initialize
    algo->init(false,        // eulerian formulation
        -1,                  //  no displacements
        -1,                  // no velocities
        nds_solidpressure,   // dof set for post processing solid pressure
        -1,                  // no scalar field
        &nearby_ele_pairs);  // possible interaction pairs

    // read the restart information, set vectors and variables
    if (restart) algo->read_restart(restart);

    // assign poro material for evaluation of porosity
    // note: to be done after potential restart, as in read_restart()
    //       the secondary material is destroyed
    PoroPressureBased::setup_material(
        context.comm, *actdis, context.cloning_material_map, struct_disname, fluid_disname);

    // 4.- Run of the actual problem.
    algo->time_loop();

    // perform the result test if required
    context.algorithm_deps.add_field_test(algo->create_field_test());
    context.run_result_tests(context.comm);
  }
}  // namespace

/*-------------------------------------------------------------------------------*
 | Main control routine for poro fluid multiphase problems           vuong 08/16 |
 *-------------------------------------------------------------------------------*/
void porofluid_pressure_based_dyn(int restart)
{
  const PorofluidPressureBasedMainContext context =
      make_main_context_from_problem(*Global::Problem::instance());
  run_porofluid_pressure_based(context, restart);
}  // poromultiphase_dyn

FOUR_C_NAMESPACE_CLOSE
