#include "4C_fsi_dyn.hpp"

#include "4C_adapter_ale_fsi.hpp"
#include "4C_adapter_fld_fluid_fsi.hpp"
#include "4C_adapter_fld_fluid_xfsi.hpp"
#include "4C_adapter_fld_moving_boundary.hpp"
#include "4C_adapter_str_fpsiwrapper.hpp"
#include "4C_adapter_str_fsiwrapper.hpp"
#include "4C_adapter_str_poro_wrapper.hpp"
#include "4C_adapter_str_structure.hpp"
#include "4C_ale_ale3.hpp"
#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_beam3_base.hpp"
#include "4C_binstrategy.hpp"
#include "4C_coupling_adapter.hpp"
#include "4C_coupling_adapter_mortar.hpp"
#include "4C_fem_condition_selector.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_dofset_fixed_size.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fluid_ele.hpp"
#include "4C_fluid_xfluid.hpp"
#include "4C_fluid_xfluid_fluid.hpp"
#include "4C_fsi_dirichletneumann_disp.hpp"
#include "4C_fsi_dirichletneumann_factory.hpp"
#include "4C_fsi_dirichletneumann_vel.hpp"
#include "4C_fsi_dirichletneumann_volcoupl.hpp"
#include "4C_fsi_dirichletneumannslideale.hpp"
#include "4C_fsi_fluid_ale.hpp"
#include "4C_fsi_fluidfluidmonolithic_fluidsplit.hpp"
#include "4C_fsi_fluidfluidmonolithic_fluidsplit_nonox.hpp"
#include "4C_fsi_fluidfluidmonolithic_structuresplit.hpp"
#include "4C_fsi_fluidfluidmonolithic_structuresplit_nonox.hpp"
#include "4C_fsi_monolithicfluidsplit.hpp"
#include "4C_fsi_monolithicstructuresplit.hpp"
#include "4C_fsi_mortarmonolithic_fluidsplit.hpp"
#include "4C_fsi_mortarmonolithic_fluidsplit_sp.hpp"
#include "4C_fsi_mortarmonolithic_structuresplit.hpp"
#include "4C_fsi_resulttest.hpp"
#include "4C_fsi_slidingmonolithic_fluidsplit.hpp"
#include "4C_fsi_slidingmonolithic_structuresplit.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_fsi_xfem_fluid.hpp"
#include "4C_fsi_xfem_monolithic.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_fbi.hpp"
#include "4C_inpar_fsi.hpp"
#include "4C_linalg_utils_sparse_algebra_math.hpp"
#include "4C_poroelast_utils_clonestrategy.hpp"
#include "4C_poroelast_utils_setup.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_so3_base.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_utils_parameter_list.hpp"
#include "4C_utils_result_test.hpp"
#include "4C_xfem_discretization.hpp"

#include <Teuchos_TimeMonitor.hpp>

#include <functional>
#include <set>
#include <string>
#include <vector>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*/
// entry point for Fluid on Ale in DRT
/*----------------------------------------------------------------------*/
void fluid_ale_drt()
{
  Global::Problem* problem = Global::Problem::instance();

  const Epetra_Comm& comm = problem->get_dis("fluid")->get_comm();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!

  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->get_dis("fluid");
  // check for xfem discretization
  if (problem->x_fluid_dynamic_params().sublist("GENERAL").get<bool>("XFLUIDFLUID"))
  {
    FLD::XFluid::setup_fluid_discretization();
  }
  else
  {
    fluiddis->fill_complete();
  }

  Teuchos::RCP<Core::FE::Discretization> aledis = problem->get_dis("ale");
  aledis->fill_complete();

  // create ale elements if the ale discretization is empty
  if (aledis->num_global_nodes() == 0)
  {
    Core::FE::clone_discretization<ALE::Utils::AleCloneStrategy>(
        *fluiddis, *aledis, Global::Problem::instance()->cloning_material_map());
    aledis->fill_complete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->evaluate(params);
  }
  else  // filled ale discretization
  {
    if (!FSI::Utils::fluid_ale_nodes_disjoint(*fluiddis, *aledis))
      FOUR_C_THROW(
          "Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
  }

  Teuchos::RCP<FSI::FluidAleAlgorithm> fluid = Teuchos::make_rcp<FSI::FluidAleAlgorithm>(comm);
  const int restart = problem->restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    fluid->read_restart(restart);
  }
  fluid->timeloop();

  Global::Problem::instance()->add_field_test(fluid->mb_fluid_field()->create_field_test());
  Global::Problem::instance()->test_all(comm);
}

/*----------------------------------------------------------------------*/
// entry point for Fluid on XFEM in DRT
/*----------------------------------------------------------------------*/
void fluid_xfem_drt()
{
  const Epetra_Comm& comm = Global::Problem::instance()->get_dis("structure")->get_comm();

  Global::Problem* problem = Global::Problem::instance();

  Teuchos::RCP<Core::FE::Discretization> soliddis = problem->get_dis("structure");
  soliddis->fill_complete();

  FLD::XFluid::setup_fluid_discretization();

  const Teuchos::ParameterList xfluid = problem->x_fluid_dynamic_params();
  bool alefluid = xfluid.sublist("GENERAL").get<bool>("ALE_XFluid");

  if (alefluid)  // in ale case
  {
    Teuchos::RCP<Core::FE::Discretization> aledis = problem->get_dis("ale");
    aledis->fill_complete();

    // create ale elements if the ale discretization is empty
    if (aledis->num_global_nodes() == 0)
    {
      Core::FE::clone_discretization<ALE::Utils::AleCloneStrategy>(
          *problem->get_dis("fluid"), *aledis, Global::Problem::instance()->cloning_material_map());
      aledis->fill_complete();
      // setup material in every ALE element
      Teuchos::ParameterList params;
      params.set<std::string>("action", "setup_material");
      aledis->evaluate(params);
    }
    else  // filled ale discretization
    {
      if (!FSI::Utils::fluid_ale_nodes_disjoint(*problem->get_dis("fluid"), *aledis))
        FOUR_C_THROW(
            "Fluid and ALE nodes have the same node numbers. "
            "This it not allowed since it causes problems with Dirichlet BCs. "
            "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");
    }
  }

  if (alefluid)
  {
    // create instance of fluid xfem algorithm, for moving interfaces
    FSI::FluidXFEMAlgorithm fluidalgo(comm);

    const int restart = Global::Problem::instance()->restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fluidalgo.read_restart(restart);
    }

    // run the simulation
    fluidalgo.timeloop();

    // perform result tests if required
    problem->add_field_test(fluidalgo.mb_fluid_field()->create_field_test());
    problem->test_all(comm);
  }
  else
  {
    //--------------------------------------------------------------
    // create instance of fluid basis algorithm
    const Teuchos::ParameterList& fdyn = Global::Problem::instance()->fluid_dynamic_params();

    Adapter::FluidBaseAlgorithm fluidalgo(fdyn, fdyn, "fluid", false);

    //--------------------------------------------------------------
    // restart the simulation
    const int restart = Global::Problem::instance()->restart();
    if (restart)
    {
      // read the restart information, set vectors and variables
      fluidalgo.fluid_field()->read_restart(restart);
    }

    //--------------------------------------------------------------
    // run the simulation
    fluidalgo.fluid_field()->integrate();

    //--------------------------------------------------------------
    // perform result tests if required
    problem->add_field_test(fluidalgo.fluid_field()->create_field_test());
    problem->test_all(comm);
  }
}

/*----------------------------------------------------------------------*/
// entry point for FSI using multidimensional immersed method (FBI)
/*----------------------------------------------------------------------*/
void fsi_immersed_drt()
{
  Global::Problem* problem = Global::Problem::instance();

  Teuchos::RCP<Core::FE::Discretization> structdis = problem->get_dis("structure");
  const Epetra_Comm& comm = structdis->get_comm();

  auto correct_node = [](const Core::Nodes::Node& node) -> decltype(auto)
  {
    const Core::Elements::Element* element = node.elements()[0];
    const auto* beamelement = dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(element);
    if (beamelement != nullptr && !beamelement->is_centerline_node(node))
      return *element->nodes()[0];
    else
      return node;
  };

  auto determine_relevant_points = [correct_node](const Core::FE::Discretization& discret,
                                       const Core::Elements::Element& ele,
                                       Teuchos::RCP<const Core::LinAlg::Vector<double>> disnp)
      -> std::vector<std::array<double, 3>>
  {
    if (dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(&ele))
    {
      return Core::Binstrategy::DefaultRelevantPoints{
          .correct_node = correct_node,
      }(discret, ele, disnp);
    }
    else
      return Core::Binstrategy::DefaultRelevantPoints{}(discret, ele, disnp);
  };

  if (structdis->get_condition("PointCoupling") != nullptr)
  {
    structdis->fill_complete(false, false, false);
    Teuchos::ParameterList binning_params = Global::Problem::instance()->binning_strategy_params();
    Core::Utils::add_enum_class_to_parameter_list<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
        binning_params);
    Core::Rebalance::rebalance_discretizations_by_binning(binning_params,
        Global::Problem::instance()->output_control_file(), {structdis}, correct_node,
        determine_relevant_points, true);
  }
  else if (not structdis->filled() || not structdis->have_dofs())
  {
    structdis->fill_complete();
  }

  problem->get_dis("fluid")->fill_complete();

  // get discretizations
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->get_dis("fluid");

  // create vector of discr.
  std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
  dis.push_back(fluiddis);
  dis.push_back(structdis);

  // binning strategy is created
  Teuchos::ParameterList binning_params = Global::Problem::instance()->binning_strategy_params();
  Core::Utils::add_enum_class_to_parameter_list<Core::FE::ShapeFunctionType>(
      "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
      binning_params);


  auto binningstrategy = Teuchos::make_rcp<Core::Binstrategy::BinningStrategy>(binning_params,
      Global::Problem::instance()->output_control_file(), comm, comm.MyPID(), correct_node,
      determine_relevant_points, dis);

  const Teuchos::ParameterList& fbidyn = problem->fbi_params();

  Inpar::FBI::BeamToFluidPreSortStrategy presort_strategy =
      Teuchos::getIntegralValue<Inpar::FBI::BeamToFluidPreSortStrategy>(fbidyn, "PRESORT_STRATEGY");

  // redistribute discr. with help of binning strategy
  if (presort_strategy == Inpar::FBI::BeamToFluidPreSortStrategy::binning)
  {
    std::vector<Teuchos::RCP<Epetra_Map>> stdelecolmap;
    std::vector<Teuchos::RCP<Epetra_Map>> stdnodecolmap;
    Teuchos::RCP<Epetra_Map> rowbins =
        binningstrategy
            ->do_weighted_partitioning_of_bins_and_extend_ghosting_of_discret_to_one_bin_layer(
                dis, stdelecolmap, stdnodecolmap);
    binningstrategy->fill_bins_into_bin_discretization(*rowbins);
    binningstrategy->fill_bins_into_bin_discretization(*rowbins);
  }



  const Teuchos::ParameterList& fsidyn = problem->fsi_dynamic_params();

  // Any partitioned algorithm.
  Teuchos::RCP<FSI::Partitioned> fsi;

  auto method = Teuchos::getIntegralValue<Inpar::FSI::PartitionedCouplingMethod>(
      fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED");
  if (method == Inpar::FSI::DirichletNeumann)
  {
    fsi = FSI::DirichletNeumannFactory::create_algorithm(comm, fsidyn);
    Teuchos::rcp_dynamic_cast<FSI::DirichletNeumann>(fsi, true)->setup();
  }
  else
    FOUR_C_THROW("unsupported partitioned FSI scheme");

  if (presort_strategy == Inpar::FBI::BeamToFluidPreSortStrategy::binning)
  {
    Teuchos::rcp_dynamic_cast<FSI::DirichletNeumannVel>(fsi, true)->set_binning(binningstrategy);
  }

  const int restart = Global::Problem::instance()->restart();
  if (restart)
  {
    // read the restart information, set vectors and variables
    fsi->read_restart(restart);
  }

  fsi->timeloop(fsi);

  // create result tests for single fields
  Global::Problem::instance()->add_field_test(fsi->mb_fluid_field()->create_field_test());
  Global::Problem::instance()->add_field_test(fsi->structure_field()->create_field_test());

  // do the actual testing
  Global::Problem::instance()->test_all(comm);
  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
}
/*----------------------------------------------------------------------*/
// entry point for FSI using ALE in DRT
/*----------------------------------------------------------------------*/
void fsi_ale_drt()
{
  Global::Problem* problem = Global::Problem::instance();

  Teuchos::RCP<Core::FE::Discretization> structdis = problem->get_dis("structure");
  const Epetra_Comm& comm = structdis->get_comm();

  // make sure the three discretizations are filled in the right order
  // this creates dof numbers with
  //
  //       structure dof < fluid dof < ale dof
  //
  // We rely on this ordering in certain non-intuitive places!
  auto correct_node = [](const Core::Nodes::Node& node) -> decltype(auto)
  {
    const Core::Elements::Element* element = node.elements()[0];
    const auto* beamelement = dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(element);
    if (beamelement != nullptr && !beamelement->is_centerline_node(node))
      return *element->nodes()[0];
    else
      return node;
  };

  auto determine_relevant_points = [correct_node](const Core::FE::Discretization& discret,
                                       const Core::Elements::Element& ele,
                                       Teuchos::RCP<const Core::LinAlg::Vector<double>> disnp)
      -> std::vector<std::array<double, 3>>
  {
    if (dynamic_cast<const Discret::ELEMENTS::Beam3Base*>(&ele))
    {
      return Core::Binstrategy::DefaultRelevantPoints{
          .correct_node = correct_node,
      }(discret, ele, disnp);
    }
    else
      return Core::Binstrategy::DefaultRelevantPoints{}(discret, ele, disnp);
  };

  if (structdis->get_condition("PointCoupling") != nullptr)
  {
    structdis->fill_complete(false, false, false);
    Teuchos::ParameterList binning_params = Global::Problem::instance()->binning_strategy_params();
    Core::Utils::add_enum_class_to_parameter_list<Core::FE::ShapeFunctionType>(
        "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
        binning_params);

    Core::Rebalance::rebalance_discretizations_by_binning(binning_params,
        Global::Problem::instance()->output_control_file(), {structdis}, correct_node,
        determine_relevant_points, true);
  }
  else if (not structdis->filled() || not structdis->have_dofs())
  {
    structdis->fill_complete();
  }

  if (problem->x_fluid_dynamic_params().sublist("GENERAL").get<bool>("XFLUIDFLUID"))
  {
    FLD::XFluid::setup_fluid_discretization();
  }
  else
    problem->get_dis("fluid")->fill_complete();

  problem->get_dis("ale")->fill_complete();

  // get discretizations
  Teuchos::RCP<Core::FE::Discretization> fluiddis = problem->get_dis("fluid");
  Teuchos::RCP<Core::FE::Discretization> aledis = problem->get_dis("ale");

  // create ale elements if the ale discretization is empty
  if (aledis->num_global_nodes() == 0)  // empty ale discretization
  {
    Core::FE::clone_discretization<ALE::Utils::AleCloneStrategy>(
        *fluiddis, *aledis, Global::Problem::instance()->cloning_material_map());
    aledis->fill_complete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->evaluate(params);
  }
  else  // filled ale discretization (i.e. read from input file)
  {
    if (!FSI::Utils::fluid_ale_nodes_disjoint(*fluiddis, *aledis))
      FOUR_C_THROW(
          "Fluid and ALE nodes have the same node numbers. "
          "This it not allowed since it causes problems with Dirichlet BCs. "
          "Use either the ALE cloning functionality or ensure non-overlapping node numbering!");

    if ((not problem->fsi_dynamic_params().get<bool>("MATCHGRID_FLUIDALE")) or
        (not problem->fsi_dynamic_params().get<bool>("MATCHGRID_STRUCTALE")))
    {
      // create vector of discr.
      std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
      dis.push_back(structdis);
      dis.push_back(fluiddis);
      dis.push_back(aledis);

      std::vector<Teuchos::RCP<Epetra_Map>> stdelecolmap;
      std::vector<Teuchos::RCP<Epetra_Map>> stdnodecolmap;

      // redistribute discr. with help of binning strategy
      if (structdis->get_comm().NumProc() > 1)
      {
        // binning strategy is created and parallel redistribution is performed
        Teuchos::ParameterList binning_params =
            Global::Problem::instance()->binning_strategy_params();
        Core::Utils::add_enum_class_to_parameter_list<Core::FE::ShapeFunctionType>(
            "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
            binning_params);
        Core::Binstrategy::BinningStrategy binningstrategy(binning_params,
            Global::Problem::instance()->output_control_file(), comm, comm.MyPID(), correct_node,
            determine_relevant_points, dis);
        binningstrategy
            .do_weighted_partitioning_of_bins_and_extend_ghosting_of_discret_to_one_bin_layer(
                dis, stdelecolmap, stdnodecolmap);
      }
    }
  }

  const Teuchos::ParameterList& fsidyn = problem->fsi_dynamic_params();

  auto coupling = Teuchos::getIntegralValue<FSI_COUPLING>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
    case fsi_iter_mortar_monolithicstructuresplit:
    case fsi_iter_mortar_monolithicfluidsplit:
    case fsi_iter_mortar_monolithicfluidsplit_saddlepoint:
    case fsi_iter_fluidfluid_monolithicfluidsplit:
    case fsi_iter_fluidfluid_monolithicstructuresplit:
    case fsi_iter_sliding_monolithicfluidsplit:
    case fsi_iter_sliding_monolithicstructuresplit:
    {
      // monolithic solver settings
      const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

      Teuchos::RCP<FSI::Monolithic> fsi;

      auto linearsolverstrategy =
          Teuchos::getIntegralValue<Inpar::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

      // call constructor to initialize the base class
      if (coupling == fsi_iter_monolithicfluidsplit)
      {
        fsi = Teuchos::make_rcp<FSI::MonolithicFluidSplit>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_monolithicstructuresplit)
      {
        fsi = Teuchos::make_rcp<FSI::MonolithicStructureSplit>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_mortar_monolithicstructuresplit)
      {
        fsi = Teuchos::make_rcp<FSI::MortarMonolithicStructureSplit>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_mortar_monolithicfluidsplit)
      {
        fsi = Teuchos::make_rcp<FSI::MortarMonolithicFluidSplit>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_mortar_monolithicfluidsplit_saddlepoint)
      {
        fsi = Teuchos::make_rcp<FSI::MortarMonolithicFluidSplitSaddlePoint>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicfluidsplit)
      {
        fsi = Teuchos::make_rcp<FSI::FluidFluidMonolithicFluidSplit>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicstructuresplit)
      {
        fsi = Teuchos::make_rcp<FSI::FluidFluidMonolithicStructureSplit>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_sliding_monolithicfluidsplit)
      {
        fsi = Teuchos::make_rcp<FSI::SlidingMonolithicFluidSplit>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_sliding_monolithicstructuresplit)
      {
        fsi = Teuchos::make_rcp<FSI::SlidingMonolithicStructureSplit>(comm, fsidyn);
      }
      else
      {
        FOUR_C_THROW(
            "Cannot find appropriate monolithic solver for coupling %d and linear strategy %d",
            coupling, linearsolverstrategy);
      }

      // read the restart information, set vectors and variables ---
      // be careful, dofmaps might be changed here in a redistribute() call
      const int restart = Global::Problem::instance()->restart();
      if (restart)
      {
        fsi->read_restart(restart);
      }

      // now do the coupling setup and create the combined dofmap
      fsi->setup_system();

      // here we go...
      fsi->timeloop(fsi);

      // calculate errors in comparison to analytical solution
      fsi->fluid_field()->calculate_error();

      // create result tests for single fields
      Global::Problem::instance()->add_field_test(fsi->ale_field()->create_field_test());
      Global::Problem::instance()->add_field_test(fsi->fluid_field()->create_field_test());
      Global::Problem::instance()->add_field_test(fsi->structure_field()->create_field_test());

      // create fsi specific result test
      Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::make_rcp<FSI::FSIResultTest>(fsi, fsidyn);
      Global::Problem::instance()->add_field_test(fsitest);

      // do the actual testing
      Global::Problem::instance()->test_all(comm);

      break;
    }
    case fsi_iter_fluidfluid_monolithicfluidsplit_nonox:
    case fsi_iter_fluidfluid_monolithicstructuresplit_nonox:
    {
      Teuchos::RCP<FSI::MonolithicNoNOX> fsi;
      if (coupling == fsi_iter_fluidfluid_monolithicfluidsplit_nonox)
      {
        fsi = Teuchos::make_rcp<FSI::FluidFluidMonolithicFluidSplitNoNOX>(comm, fsidyn);
      }
      else if (coupling == fsi_iter_fluidfluid_monolithicstructuresplit_nonox)
      {
        fsi = Teuchos::make_rcp<FSI::FluidFluidMonolithicStructureSplitNoNOX>(comm, fsidyn);
      }
      else
        FOUR_C_THROW("Unsupported monolithic XFFSI scheme");

      // read the restart information, set vectors and variables ---
      // be careful, dofmaps might be changed here in a redistribute() call
      const int restart = Global::Problem::instance()->restart();
      if (restart)
      {
        fsi->read_restart(restart);
      }

      // now do the coupling setup and create the combined dofmap
      fsi->setup_system();

      // here we go...
      fsi->timeloop();

      // calculate errors in comparison to analytical solution
      fsi->fluid_field()->calculate_error();

      // create result tests for single fields
      Global::Problem::instance()->add_field_test(fsi->fluid_field()->create_field_test());
      Global::Problem::instance()->add_field_test(fsi->structure_field()->create_field_test());

      // create fsi specific result test
      Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::make_rcp<FSI::FSIResultTest>(fsi, fsidyn);
      Global::Problem::instance()->add_field_test(fsitest);

      // do the actual testing
      Global::Problem::instance()->test_all(comm);

      break;
    }
    default:
    {
      // Any partitioned algorithm.

      Teuchos::RCP<FSI::Partitioned> fsi;

      auto method = Teuchos::getIntegralValue<Inpar::FSI::PartitionedCouplingMethod>(
          fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED");

      switch (method)
      {
        case Inpar::FSI::DirichletNeumann:
        case Inpar::FSI::DirichletNeumannSlideale:
        case Inpar::FSI::DirichletNeumannVolCoupl:
          fsi = FSI::DirichletNeumannFactory::create_algorithm(comm, fsidyn);
          Teuchos::rcp_dynamic_cast<FSI::DirichletNeumann>(fsi, true)->setup();
          break;
        default:
          FOUR_C_THROW("unsupported partitioned FSI scheme");
          break;
      }
      const int restart = Global::Problem::instance()->restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        fsi->read_restart(restart);
      }

      fsi->timeloop(fsi);

      // create result tests for single fields
      Global::Problem::instance()->add_field_test(fsi->mb_fluid_field()->create_field_test());
      Global::Problem::instance()->add_field_test(fsi->structure_field()->create_field_test());

      // do the actual testing
      Global::Problem::instance()->test_all(comm);

      break;
    }
  }

  Teuchos::TimeMonitor::summarize(std::cout, false, true, false);
}

/*----------------------------------------------------------------------*/
// entry point for FSI using XFEM in DRT (also for ale case)
/*----------------------------------------------------------------------*/
void xfsi_drt()
{
  const Epetra_Comm& comm = Global::Problem::instance()->get_dis("structure")->get_comm();

  if (comm.MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "       @..@    " << std::endl;
    std::cout << "      (----)      " << std::endl;
    std::cout << "     ( >__< )   " << std::endl;
    std::cout << "     ^^ ~~ ^^  " << std::endl;
    std::cout << "     _     _ _______ _______ _____" << std::endl;
    std::cout << "      \\\\__/  |______ |______   |  " << std::endl;
    std::cout << "     _/  \\\\_ |       ______| __|__" << std::endl;
    std::cout << std::endl << std::endl;
  }

  Global::Problem* problem = Global::Problem::instance();
  const Teuchos::ParameterList& fsidyn = problem->fsi_dynamic_params();

  Teuchos::RCP<Core::FE::Discretization> soliddis = problem->get_dis("structure");
  soliddis->fill_complete();

  FLD::XFluid::setup_fluid_discretization();

  Teuchos::RCP<Core::FE::Discretization> fluiddis = Global::Problem::instance()->get_dis(
      "fluid");  // at the moment, 'fluid'-discretization is used for ale!!!

  // CREATE ALE
  const Teuchos::ParameterList& xfdyn = problem->x_fluid_dynamic_params();
  bool ale = xfdyn.sublist("GENERAL").get<bool>("ALE_XFluid");
  Teuchos::RCP<Core::FE::Discretization> aledis;
  if (ale)
  {
    aledis = problem->get_dis("ale");
    if (aledis == Teuchos::null) FOUR_C_THROW("XFSI DYNAMIC: ALE discretization empty!!!");

    aledis->fill_complete(true, true, true);

    // Create ALE elements if the ale discretization is empty
    if (aledis->num_global_nodes() == 0)  // ALE discretization still empty
    {
      Core::FE::clone_discretization<ALE::Utils::AleCloneStrategy>(
          *fluiddis, *aledis, Global::Problem::instance()->cloning_material_map());
      aledis->fill_complete();
      // setup material in every ALE element
      Teuchos::ParameterList params;
      params.set<std::string>("action", "setup_material");
      aledis->evaluate(params);
    }
    else  // ALE discretization already filled
    {
      if (!FSI::Utils::fluid_ale_nodes_disjoint(*fluiddis, *aledis))
        FOUR_C_THROW(
            "Fluid and ALE nodes have the same node numbers. "
            "This it not allowed since it causes problems with Dirichlet BCs. "
            "Use the ALE cloning functionality or ensure non-overlapping node numbering!");
    }
  }

  const auto coupling = Teuchos::getIntegralValue<FsiCoupling>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_xfem_monolithic:
    {
      // monolithic solver settings
      const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");

      auto linearsolverstrategy =
          Teuchos::getIntegralValue<Inpar::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

      if (linearsolverstrategy != Inpar::FSI::PreconditionedKrylov)
        FOUR_C_THROW("Only Newton-Krylov scheme with XFEM fluid");

      // create the MonolithicXFEM object that does the whole work
      FSI::MonolithicXFEM fsi(comm, fsidyn);

      // read the restart information, set vectors and variables ---
      // be careful, dofmaps might be changed here in a redistribute() call
      const int restart = Global::Problem::instance()->restart();
      if (restart)
      {
        fsi.read_restart(restart);
      }

      // setup the system (block-DOF-row maps, systemmatrix etc.) for the monolithic XFEM system
      fsi.setup_system();

      // here we go...
      fsi.timeloop();

      Global::Problem::instance()->add_field_test(fsi.fluid_field()->create_field_test());
      fsi.structure_poro()->test_results(Global::Problem::instance());

      //    // create FSI specific result test
      //    Teuchos::RCP<FSI::FSIResultTest> fsitest = Teuchos::rcp(new
      //    FSI::FSIResultTest(fsi,fsidyn)); Global::Problem::instance()->AddFieldTest(fsitest);

      // do the actual testing
      Global::Problem::instance()->test_all(comm);

      break;
    }
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
      FOUR_C_THROW("Unreasonable choice");
      break;
    default:
    {
      // Any partitioned algorithm. Stable of working horses.

      Teuchos::RCP<FSI::Partitioned> fsi;

      auto method = Teuchos::getIntegralValue<Inpar::FSI::PartitionedCouplingMethod>(
          fsidyn.sublist("PARTITIONED SOLVER"), "PARTITIONED");

      switch (method)
      {
        case Inpar::FSI::DirichletNeumann:
          fsi = FSI::DirichletNeumannFactory::create_algorithm(comm, fsidyn);
          Teuchos::rcp_dynamic_cast<FSI::DirichletNeumann>(fsi, true)->setup();
          break;
        default:
          FOUR_C_THROW("only Dirichlet-Neumann partitioned schemes with XFEM");
          break;
      }

      const int restart = Global::Problem::instance()->restart();
      if (restart)
      {
        // read the restart information, set vectors and variables
        fsi->read_restart(restart);
      }

      fsi->timeloop(fsi);

      Global::Problem::instance()->add_field_test(fsi->mb_fluid_field()->create_field_test());
      Global::Problem::instance()->add_field_test(fsi->structure_field()->create_field_test());
      Global::Problem::instance()->test_all(comm);

      break;
    }
  }

  Teuchos::TimeMonitor::summarize();
}

/*----------------------------------------------------------------------*/
// entry point for FPSI using XFEM in DRT
/*----------------------------------------------------------------------*/
void xfpsi_drt()
{
  const Epetra_Comm& comm = Global::Problem::instance()->get_dis("structure")->get_comm();

  if (comm.MyPID() == 0)
  {
    std::cout << std::endl;
    std::cout << "       @..@    " << std::endl;
    std::cout << "      (----)      " << std::endl;
    std::cout << "     ( >__< )   " << std::endl;
    std::cout << "     ^^ ~~ ^^  " << std::endl;
    std::cout << "     _     _ _______  ______  _______ _____" << std::endl;
    std::cout << "      \\\\__/  |______ ||____|  |______   |  " << std::endl;
    std::cout << "     _/  \\\\_ |       ||       ______| __|__" << std::endl;
    std::cout << std::endl << std::endl;
  }
  Global::Problem* problem = Global::Problem::instance();

  // 1.-Initialization.
  // setup of the discretizations, including clone strategy
  PoroElast::Utils::setup_poro<PoroElast::Utils::PoroelastCloneStrategy>();

  // setup of discretization for xfluid
  FLD::XFluid::setup_fluid_discretization();
  Teuchos::RCP<Core::FE::Discretization> fluiddis = Global::Problem::instance()->get_dis(
      "fluid");  // at the moment, 'fluid'-discretization is used for ale!!!

  Teuchos::RCP<Core::FE::Discretization> aledis;
  const Teuchos::ParameterList& xfdyn = problem->x_fluid_dynamic_params();
  bool ale = xfdyn.sublist("GENERAL").get<bool>("ALE_XFluid");
  if (ale)
  {
    aledis = problem->get_dis("ale");
    if (aledis == Teuchos::null) FOUR_C_THROW("Ale discretization empty!");

    aledis->fill_complete(true, true, true);

    // 3.- Create ALE elements if the ale discretization is empty
    if (aledis->num_global_nodes() == 0)  // ALE discretization still empty
    {
      Core::FE::clone_discretization<ALE::Utils::AleCloneStrategy>(
          *fluiddis, *aledis, Global::Problem::instance()->cloning_material_map());
      aledis->fill_complete();
      // setup material in every ALE element
      Teuchos::ParameterList params;
      params.set<std::string>("action", "setup_material");
      aledis->evaluate(params);
    }
    else  // ALE discretization already filled
    {
      if (!FSI::Utils::fluid_ale_nodes_disjoint(*fluiddis, *aledis))
        FOUR_C_THROW(
            "Fluid and ALE nodes have the same node numbers. "
            "This it not allowed since it causes problems with Dirichlet BCs. "
            "Use the ALE cloning functionality or ensure non-overlapping node numbering!");
    }
  }

  // print all dofsets
  fluiddis->get_dof_set_proxy()->print_all_dofsets(fluiddis->get_comm());

  // 2.- Parameter reading
  const Teuchos::ParameterList& fsidyn = problem->fsi_dynamic_params();
  const auto coupling = Teuchos::getIntegralValue<FsiCoupling>(fsidyn, "COUPALGO");
  switch (coupling)
  {
    case fsi_iter_xfem_monolithic:
    {
      // monolithic solver settings
      const Teuchos::ParameterList& fsimono = fsidyn.sublist("MONOLITHIC SOLVER");
      auto linearsolverstrategy =
          Teuchos::getIntegralValue<Inpar::FSI::LinearBlockSolver>(fsimono, "LINEARBLOCKSOLVER");

      if (linearsolverstrategy != Inpar::FSI::PreconditionedKrylov)
        FOUR_C_THROW("Only Newton-Krylov scheme with XFEM fluid");

      FSI::MonolithicXFEM fsi(comm, fsidyn, Adapter::FieldWrapper::type_PoroField);

      // read the restart information, set vectors and variables ---

      // be careful, dofmaps might be changed here in a redistribute() call
      const int restart =
          Global::Problem::instance()
              ->restart();  // not adapated at the moment .... Todo check it .. ChrAg
      if (restart)
      {
        fsi.read_restart(restart);
      }

      fsi.setup_system();


      // 3.2.- redistribute the FPSI interface
      // Todo .... fsi->redistribute_interface(); // this is required for paralles fpi-condition
      // (not included in this commit)

      // here we go...
      fsi.timeloop();

      Global::Problem::instance()->add_field_test(fsi.fluid_field()->create_field_test());
      fsi.structure_poro()->test_results(Global::Problem::instance());

      // do the actual testing
      Global::Problem::instance()->test_all(comm);
      break;
    }
    case fsi_iter_monolithicfluidsplit:
    case fsi_iter_monolithicstructuresplit:
      FOUR_C_THROW("Unreasonable choice");
      break;
    default:
    {
      FOUR_C_THROW("FPSI_XFEM: No Partitioned Algorithms implemented !!!");
      break;
    }
  }
  Teuchos::TimeMonitor::summarize();
}

FOUR_C_NAMESPACE_CLOSE
