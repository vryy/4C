/*----------------------------------------------------------------------*/
/*! \file
\brief Control routine for structure with ale problems.

\level 2


*/
/*----------------------------------------------------------------------*/

/*----------------------------------------------------------------------*
 |  headers                                                  mgit 04/11 |
 *----------------------------------------------------------------------*/
#include "4C_wear_dyn.hpp"

#include "4C_ale_utils_clonestrategy.hpp"
#include "4C_fem_condition_utils.hpp"
#include "4C_fem_general_utils_createdis.hpp"
#include "4C_fsi_utils.hpp"
#include "4C_global_data.hpp"
#include "4C_inpar_wear.hpp"
#include "4C_rebalance_binning_based.hpp"
#include "4C_so3_base.hpp"
#include "4C_solid_3D_ele.hpp"
#include "4C_wear_partitioned.hpp"

#include <Epetra_MpiComm.h>
#include <Teuchos_TimeMonitor.hpp>

FOUR_C_NAMESPACE_OPEN

/*----------------------------------------------------------------------*
 | entry point for structure ale in DRT                      mgit 04/11 |
 *----------------------------------------------------------------------*/
void wear_dyn_drt(int restart)
{
  // create a communicator
  const Epetra_Comm& comm = Global::Problem::instance()->get_dis("structure")->get_comm();

  // ***********************************************************
  // Setup the problem
  // ***********************************************************
  const Teuchos::ParameterList& sdyn = Global::Problem::instance()->structural_dynamic_params();
  const Teuchos::ParameterList& wearpara = Global::Problem::instance()->wear_params();

  // check if quasistatic analysis
  if (sdyn.get<std::string>("DYNAMICTYP") != "Statics")
  {
    std::cout << "WARNING: wear without dynamic effects!!!" << std::endl;
    // FOUR_C_THROW ("ERROR: Structure with ale only for quasistatic analysis so in new sti so
    // far.");
  }

  // access the structure discretization, make sure it is filled
  Teuchos::RCP<Core::FE::Discretization> structdis = Teuchos::null;
  structdis = Global::Problem::instance()->get_dis("structure");
  // set degrees of freedom in the discretization
  if (!structdis->filled() or !structdis->have_dofs()) structdis->fill_complete();

  // access the ale discretization
  Teuchos::RCP<Core::FE::Discretization> aledis = Teuchos::null;
  aledis = Global::Problem::instance()->get_dis("ale");
  if (!aledis->filled()) aledis->fill_complete();

  // we use the structure discretization as layout for the ale discretization
  if (structdis->num_global_nodes() == 0) FOUR_C_THROW("Structure discretization is empty!");

  // clone ale mesh from structure discretization
  if (aledis->num_global_nodes() == 0)
  {
    Core::FE::CloneDiscretization<ALE::UTILS::AleCloneStrategy>(
        structdis, aledis, Global::Problem::instance()->cloning_material_map());
    aledis->fill_complete();
    // setup material in every ALE element
    Teuchos::ParameterList params;
    params.set<std::string>("action", "setup_material");
    aledis->evaluate(params);
  }
  else  // filled ale discretization
  {
    // if we have non-matching meshes:
    if (!Core::UTILS::IntegralValue<bool>(
            Global::Problem::instance()->wear_params(), "MATCHINGGRID"))
    {
      // create vector of discr.
      std::vector<Teuchos::RCP<Core::FE::Discretization>> dis;
      dis.push_back(structdis);
      dis.push_back(aledis);

      Teuchos::ParameterList binning_params =
          Global::Problem::instance()->binning_strategy_params();
      Core::UTILS::AddEnumClassToParameterList<Core::FE::ShapeFunctionType>(
          "spatial_approximation_type", Global::Problem::instance()->spatial_approximation_type(),
          binning_params);
      auto element_filter = [](const Core::Elements::Element* element)
      { return Core::Binstrategy::Utils::SpecialElement::none; };
      auto rigid_sphere_radius = [](const Core::Elements::Element* element) { return 0.0; };
      auto correct_beam_center_node = [](const Core::Nodes::Node* node) { return node; };
      Core::Rebalance::RebalanceDiscretizationsByBinning(binning_params,
          Global::Problem::instance()->output_control_file(), dis, element_filter,
          rigid_sphere_radius, correct_beam_center_node, false);
    }
  }
  // ***********************************************************

  Teuchos::RCP<Wear::Algorithm> stru_ale = Teuchos::null;

  // structure ale object
  if (Core::UTILS::IntegralValue<Inpar::Wear::WearCoupAlgo>(wearpara, "WEAR_COUPALGO") ==
      Inpar::Wear::wear_stagg)
  {
    stru_ale = Teuchos::rcp(new Wear::Partitioned(comm));
  }
  else if (Core::UTILS::IntegralValue<Inpar::Wear::WearCoupAlgo>(wearpara, "WEAR_COUPALGO") ==
           Inpar::Wear::wear_iterstagg)
  {
    stru_ale = Teuchos::rcp(new Wear::Partitioned(comm));
  }
  else
  {
    FOUR_C_THROW("Chosen algorithm not supported");
  }

  // read restart before joining the time loop
  if (restart != 0) stru_ale->read_restart(restart);

  // solve the whole problem
  stru_ale->time_loop();

  // summarize the performance measurements
  Teuchos::TimeMonitor::summarize();

  // perform the result test
  Global::Problem::instance()->add_field_test(stru_ale->structure_field()->create_field_test());
  Global::Problem::instance()->test_all(comm);

  return;
}  // wear_dyn_drt()

/*----------------------------------------------------------------------*/

FOUR_C_NAMESPACE_CLOSE
