/*-----------------------------------------------------------*/
/*! \file
\brief Structure specific result test class


\level 3

*/
/*-----------------------------------------------------------*/

#include "4C_structure_new_resulttest.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_linalg_fixedsizematrix_voigt_notation.hpp"
#include "4C_structure_new_model_evaluator_data.hpp"
#include "4C_structure_new_timint_base.hpp"
#include "4C_utils_exceptions.hpp"

#include <Epetra_MultiVector.h>
#include <Teuchos_RCPDecl.hpp>

#include <optional>
#include <string>

FOUR_C_NAMESPACE_OPEN

namespace
{
  struct QuantityNameAndComponent
  {
    std::string name;
    int component;
  };

  [[nodiscard]] std::optional<QuantityNameAndComponent> GetGaussPointDataNameAndComponent(
      const std::string& name_with_component,
      const std::unordered_map<std::string, int>& quantities)
  {
    for (const auto& [quantity_name, num_quantities] : quantities)
    {
      if (num_quantities == 1)
      {
        if (name_with_component == quantity_name)
          return std::make_optional<QuantityNameAndComponent>({quantity_name, 0});
        else
          continue;
      }

      for (auto i = 0; i < num_quantities; ++i)
      {
        if (name_with_component == quantity_name + "_" + std::to_string(i + 1))
        {
          return std::make_optional<QuantityNameAndComponent>({quantity_name, i});
        }
      }
    }
    return std::nullopt;
  }

  [[nodiscard]] double GetGaussPointDataValue(const QuantityNameAndComponent& name_and_component,
      int node_id,
      const std::unordered_map<std::string, Teuchos::RCP<Epetra_MultiVector>>& all_data)
  {
    const Epetra_MultiVector& data = *all_data.at(name_and_component.name);

    int local_id = data.Map().LID(node_id);

    if (local_id < 0)
    {
      FOUR_C_THROW("You tried to test %s on a proc that does not own node %i.",
          name_and_component.name.c_str(), node_id);
    }

    return data[name_and_component.component][local_id];
  }

  /*!
   * @brief Returns the stress or strain component requested in label at the given node
   *
   * @param prefix (in) : prefix (either stress or strain)
   * @param label (in) : label of the quantity, e.g. stress_xx
   * @param node_id (in) : Id of the node
   * @param nodal_data (in) : Nodal data
   * @return double
   */
  double GetNodalStressStrainComponent(const std::string& prefix, const std::string& label,
      int node_id, const Epetra_MultiVector& nodal_data)
  {
    int voigt_index = -1;
    if (label == prefix + "_xx")
    {
      voigt_index = Core::LinAlg::Voigt::IndexMappings::symmetric_tensor_to_voigt6_index(0, 0);
    }
    else if (label == prefix + "_yy")
    {
      voigt_index = Core::LinAlg::Voigt::IndexMappings::symmetric_tensor_to_voigt6_index(1, 1);
    }
    else if (label == prefix + "_zz")
    {
      voigt_index = Core::LinAlg::Voigt::IndexMappings::symmetric_tensor_to_voigt6_index(2, 2);
    }
    else if (label == prefix + "_xy")
    {
      voigt_index = Core::LinAlg::Voigt::IndexMappings::symmetric_tensor_to_voigt6_index(0, 1);
    }
    else if (label == prefix + "_xz")
    {
      voigt_index = Core::LinAlg::Voigt::IndexMappings::symmetric_tensor_to_voigt6_index(0, 2);
    }
    else if (label == prefix + "_yz")
    {
      voigt_index = Core::LinAlg::Voigt::IndexMappings::symmetric_tensor_to_voigt6_index(1, 2);
    }

    if (voigt_index < 0)
    {
      FOUR_C_THROW(
          "You try to test an unknown %s component %s. Use one of %s{_xx, _yy, _zz, _xy, _xz, "
          "_yz}.",
          label.c_str(), prefix.c_str(), prefix.c_str());
    }

    int local_id = nodal_data.Map().LID(node_id);

    if (local_id < 0)
    {
      FOUR_C_THROW(
          "You tried to test %s on a proc that does not own node %i.", label.c_str(), node_id);
    }

    return nodal_data[voigt_index][local_id];
  }
}  // namespace

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ResultTest::ResultTest()
    : Core::UTILS::ResultTest("STRUCTURE"),
      isinit_(false),
      issetup_(false),
      strudisc_(Teuchos::null),
      disn_(Teuchos::null),
      dismatn_(Teuchos::null),
      veln_(Teuchos::null),
      accn_(Teuchos::null),
      reactn_(Teuchos::null),
      gstate_(Teuchos::null)
{
  // empty constructor
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ResultTest::init(
    const STR::TimeInt::BaseDataGlobalState& gstate, const STR::MODELEVALUATOR::Data& data)
{
  issetup_ = false;

  disn_ = gstate.get_dis_n();
  veln_ = gstate.get_vel_n();
  accn_ = gstate.get_acc_n();
  reactn_ = gstate.get_freact_n();
  gstate_ = Teuchos::rcpFromRef(gstate);
  data_ = Teuchos::rcpFromRef(data);
  strudisc_ = gstate.get_discret();

  if (Global::Problem::Instance()->GetProblemType() == Core::ProblemType::struct_ale and
      (Global::Problem::Instance()->WearParams()).get<double>("WEARCOEFF") > 0.0)
    FOUR_C_THROW("Material displ. are not yet considered!");
  else
    dismatn_ = Teuchos::null;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ResultTest::setup()
{
  check_init();
  // currently unused
  issetup_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::ResultTest::test_node(Input::LineDefinition& res, int& nerr, int& test_count)
{
  check_init_setup();

  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.extract_string("DIS", dis);
  if (dis != strudisc_->Name()) return;

  int node;
  res.extract_int("NODE", node);
  node -= 1;

  int havenode(strudisc_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  strudisc_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW(
        "Node %d does not belong to discretization %s", node + 1, strudisc_->Name().c_str());
  }

  std::string position;
  res.extract_string("QUANTITY", position);

  if (strudisc_->HaveGlobalNode(node))
  {
    double result;
    int error_code = get_nodal_result(result, node, position);

    if (error_code == 0)
    {
      // compare values
      const int err = compare_values(result, "NODE", res);
      nerr += err;
      test_count++;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
int STR::ResultTest::get_nodal_result(
    double& result, const int node, const std::string& position) const
{
  result = 0.0;

  const Core::Nodes::Node* actnode = strudisc_->gNode(node);

  // Here we are just interested in the nodes that we own (i.e. a row node)!
  if (actnode->Owner() != strudisc_->Comm().MyPID()) return -1;

  bool unknownpos = true;  // make sure the result value std::string can be handled

  // test displacements or pressure
  if (disn_ != Teuchos::null)
  {
    const Epetra_BlockMap& disnpmap = disn_->Map();
    int idx = -1;
    if (position == "dispx")
      idx = 0;
    else if (position == "dispy")
      idx = 1;
    else if (position == "dispz")
      idx = 2;
    else if (position == "press")
      idx = 3;

    if (idx >= 0)
    {
      unknownpos = false;
      int lid = disnpmap.LID(strudisc_->Dof(0, actnode, idx));
      if (lid < 0)
        FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
            actnode->Id());
      result = (*disn_)[lid];
    }
  }

  // test material displacements
  if (!dismatn_.is_null())
  {
    const Epetra_BlockMap& dismpmap = dismatn_->Map();
    int idx = -1;
    if (position == "dispmx")
      idx = 0;
    else if (position == "dispmy")
      idx = 1;
    else if (position == "dispmz")
      idx = 2;

    if (idx >= 0)
    {
      unknownpos = false;
      int lid = dismpmap.LID(strudisc_->Dof(0, actnode, idx));
      if (lid < 0)
        FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
            actnode->Id());
      result = (*dismatn_)[lid];
    }
  }

  // test velocities
  if (veln_ != Teuchos::null)
  {
    const Epetra_BlockMap& velnpmap = veln_->Map();
    int idx = -1;
    if (position == "velx")
      idx = 0;
    else if (position == "vely")
      idx = 1;
    else if (position == "velz")
      idx = 2;

    if (idx >= 0)
    {
      unknownpos = false;
      int lid = velnpmap.LID(strudisc_->Dof(0, actnode, idx));
      if (lid < 0)
        FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
            actnode->Id());
      result = (*veln_)[lid];
    }
  }

  // test accelerations
  if (accn_ != Teuchos::null)
  {
    const Epetra_BlockMap& accnpmap = accn_->Map();
    int idx = -1;
    if (position == "accx")
      idx = 0;
    else if (position == "accy")
      idx = 1;
    else if (position == "accz")
      idx = 2;

    if (idx >= 0)
    {
      unknownpos = false;
      int lid = accnpmap.LID(strudisc_->Dof(0, actnode, idx));
      if (lid < 0)
        FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
            actnode->Id());
      result = (*accn_)[lid];
    }
  }

  // test nodal stresses
  if (position.rfind("stress", 0) == 0)
  {
    if (data_->get_stress_data_node_postprocessed() == Teuchos::null)
    {
      FOUR_C_THROW(
          "It looks like you don't write stresses. You have to specify the stress type in "
          "IO->STRUCT_STRESS");
    }
    result = GetNodalStressStrainComponent(
        "stress", position, node, *data_->get_stress_data_node_postprocessed());
    unknownpos = false;
  }

  // test nodal strain
  if (position.rfind("strain", 0) == 0)
  {
    if (data_->get_stress_data_node_postprocessed() == Teuchos::null)
    {
      FOUR_C_THROW(
          "It looks like you don't write strains. You have to specify the strain type in "
          "IO->STRUCT_STRAIN");
    }
    result = GetNodalStressStrainComponent(
        "strain", position, node, *data_->get_strain_data_node_postprocessed());
    unknownpos = false;
  }

  // test for any postprocessed gauss point data
  if (data_->get_gauss_point_data_output_manager_ptr() != Teuchos::null)
  {
    std::optional<QuantityNameAndComponent> name_and_component = GetGaussPointDataNameAndComponent(
        position, data_->get_gauss_point_data_output_manager_ptr()->get_quantities());
    if (name_and_component.has_value())
    {
      result = GetGaussPointDataValue(*name_and_component, node,
          data_->get_gauss_point_data_output_manager_ptr()->get_nodal_data());
      unknownpos = false;
    }
  }

  // test reaction
  if (reactn_ != Teuchos::null)
  {
    const Epetra_BlockMap& reactmap = reactn_->Map();
    int idx = -1;
    if (position == "reactx")
      idx = 0;
    else if (position == "reacty")
      idx = 1;
    else if (position == "reactz")
      idx = 2;

    if (idx >= 0)
    {
      unknownpos = false;
      int lid = reactmap.LID(strudisc_->Dof(0, actnode, idx));
      if (lid < 0)
        FOUR_C_THROW("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
            actnode->Id());
      result = (*reactn_)[lid];
    }
  }

  // catch position std::strings, which are not handled by structure result test
  if (unknownpos)
    FOUR_C_THROW("Quantity '%s' not supported in structure testing", position.c_str());

  return 0;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::ResultTest::test_node_on_geometry(Input::LineDefinition& res, int& nerr, int& test_count,
    const std::vector<std::vector<std::vector<int>>>& nodeset)
{
  check_init_setup();

  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.extract_string("DIS", dis);
  if (dis != strudisc_->Name()) return;

  auto get_geometry_info = [res]() -> std::tuple<int, int>
  {
    int geometry_type = 0, geometry_id = 0;
    if (res.has_named("NODE"))
    {
      res.extract_int("NODE", geometry_id);
      geometry_id -= 1;
      geometry_type = 0;
    }
    else if (res.has_named("LINE"))
    {
      res.extract_int("LINE", geometry_id);
      geometry_id -= 1;
      geometry_type = 1;
    }
    else if (res.has_named("SURFACE"))
    {
      res.extract_int("SURFACE", geometry_id);
      geometry_id -= 1;
      geometry_type = 2;
    }
    else if (res.has_named("VOLUME"))
    {
      res.extract_int("VOLUME", geometry_id);
      geometry_id -= 1;
      geometry_type = 3;
    }
    else
      FOUR_C_THROW("Invalid line definition found");

    return {geometry_type, geometry_id};
  };
  const auto [geometry_type, geometry_id] = get_geometry_info();

  FOUR_C_ASSERT(geometry_type >= 0 && geometry_type <= 3, "The geometry type is invalid");
  if (geometry_id < 0 || geometry_id >= static_cast<int>(nodeset[geometry_type].size()))
    FOUR_C_THROW("Invalid geometry id %d", geometry_id);

  const std::vector<int>& nodes = nodeset[geometry_type][geometry_id];

  std::string op_str;
  TestOp op = TestOp::unknown;
  res.extract_string("OP", op_str);
  if (op_str == "sum")
    op = TestOp::sum;
  else if (op_str == "max")
    op = TestOp::max;
  else if (op_str == "min")
    op = TestOp::min;
  else
    FOUR_C_THROW("Invalid operation %s", op_str.c_str());

  std::string position;
  res.extract_string("QUANTITY", position);

  // collect the local result
  double tmp_result;
  switch (op)
  {
    case TestOp::max:
      tmp_result = -1e99;
      break;
    case TestOp::min:
      tmp_result = 1e99;
      break;
    default:
      tmp_result = 0.0;
      break;
  }
  for (int node : nodes)
  {
    double tmp = 0.0;
    if (strudisc_->HaveGlobalNode(node)) get_nodal_result(tmp, node, position);

    switch (op)
    {
      case TestOp::sum:
        tmp_result += tmp;
        break;
      case TestOp::max:
        if (tmp > tmp_result) tmp_result = tmp;
        break;
      case TestOp::min:
        if (tmp < tmp_result) tmp_result = tmp;
        break;
      default:
        break;
    }
  }

  // gather the result across processes
  auto gather_result = [op](const Core::FE::Discretization& disc,
                           const double local_result) -> double
  {
    double tmp_result = local_result, result = 0.0;
    switch (op)
    {
      case TestOp::sum:
        disc.Comm().SumAll(&tmp_result, &result, 1);
        break;
      case TestOp::max:
        disc.Comm().MaxAll(&tmp_result, &result, 1);
        break;
      case TestOp::min:
        disc.Comm().MinAll(&tmp_result, &result, 1);
        break;
      default:
        break;
    }
    return result;
  };
  const double result = gather_result(*strudisc_, tmp_result);

  // test the value; we only do at rank 0
  if (strudisc_->Comm().MyPID() == 0)
  {
    int err = 0;
    switch (geometry_type)
    {
      case 0:
        err = compare_values(result, "NODE", res);
        break;
      case 1:
        err = compare_values(result, "LINE", res);
        break;
      case 2:
        err = compare_values(result, "SURFACE", res);
        break;
      case 3:
        err = compare_values(result, "VOLUME", res);
        break;
      default:
        break;
    }
    nerr += err;
    test_count++;
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ResultTest::TestSpecial(
    Input::LineDefinition& res, int& nerr, int& test_count, int& uneval_test_count)
{
  check_init_setup();

  std::string quantity;
  res.extract_string("QUANTITY", quantity);

  Status special_status = Status::unevaluated;
  const std::optional<double> result = get_special_result(quantity, special_status);

  switch (special_status)
  {
    case Status::evaluated:
    {
      if (result.has_value())
        nerr += compare_values(*result, "SPECIAL", res);
      else
        FOUR_C_THROW(
            "STR::ResultTest::TestSpecial: Special result has no defined value assigned to it!");
      ++test_count;
      break;
    }
    case Status::unevaluated:
    {
      ++uneval_test_count;
      break;
    }
    default:
    {
      FOUR_C_THROW(
          "STR::ResultTest::TestSpecial: Undefined status type (enum=%d)!", special_status);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::optional<double> STR::ResultTest::get_special_result(
    const std::string& quantity, Status& special_status) const
{
  if (quantity.find("num_iter_step_") != quantity.npos)
  {
    return get_nln_iteration_number(quantity, special_status);
  }
  else if (quantity.find("lin_iter_step_") != quantity.npos)
  {
    return get_last_lin_iteration_number(quantity, special_status);
  }
  else if (quantity.find("nodes_proc") != quantity.npos)
  {
    return get_nodes_per_proc_number(quantity, special_status);
  }
  else if (quantity == "internal_energy" or quantity == "kinetic_energy" or
           quantity == "total_energy" or quantity == "beam_contact_penalty_potential" or
           quantity == "beam_interaction_potential" or
           quantity == "beam_to_beam_link_internal_energy" or
           quantity == "beam_to_beam_link_kinetic_energy" or
           quantity == "beam_to_sphere_link_internal_energy" or
           quantity == "beam_to_sphere_link_kinetic_energy")
  {
    return get_energy(quantity, special_status);
  }
  else
    FOUR_C_THROW(
        "Quantity '%s' not supported by special result testing functionality "
        "for structure field!",
        quantity.c_str());

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::optional<int> STR::ResultTest::get_last_lin_iteration_number(
    const std::string& quantity, Status& special_status) const
{
  std::optional<int> result = std::nullopt;

  if (strudisc_->Comm().MyPID() == 0)
  {
    const int stepn = GetIntegerNumberAtLastPositionOfName(quantity);

    const int restart = Global::Problem::Instance()->Restart();
    if (stepn <= restart) return -1;

    special_status = Status::evaluated;
    result = gstate_->get_last_lin_iteration_number(stepn);
  }

  return result;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::optional<int> STR::ResultTest::get_nln_iteration_number(
    const std::string& quantity, Status& special_status) const
{
  std::optional<int> result = std::nullopt;

  if (strudisc_->Comm().MyPID() == 0)
  {
    const int stepn = GetIntegerNumberAtLastPositionOfName(quantity);

    const int restart = Global::Problem::Instance()->Restart();
    if (stepn <= restart) return -1;

    special_status = Status::evaluated;
    result = gstate_->get_nln_iteration_number(stepn);
  }

  return result;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::optional<int> STR::ResultTest::get_nodes_per_proc_number(
    const std::string& quantity, Status& special_status) const
{
  std::optional<int> result = std::nullopt;

  std::string proc_string = quantity.substr(quantity.find("nodes_proc") + 10);
  const int proc_num = std::stoi(proc_string);

  // extract processor ID
  if (proc_num >= strudisc_->Comm().NumProc())
    FOUR_C_THROW("STR::ResultTest::get_nodes_per_proc_number: Invalid processor ID!");

  if (strudisc_->Comm().MyPID() == proc_num)
  {
    // extract number of nodes owned by specified processor
    special_status = Status::evaluated;
    result = strudisc_->NodeRowMap()->NumMyElements();
  }

  return result;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
std::optional<double> STR::ResultTest::get_energy(
    const std::string& quantity, Status& special_status) const
{
  std::optional<double> result = std::nullopt;

  if (strudisc_->Comm().MyPID() == 0)
  {
    special_status = Status::evaluated;
    result = data_->get_energy_data(quantity);
  }

  return result;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::GetIntegerNumberAtLastPositionOfName(const std::string& quantity)
{
  std::stringstream ss(quantity);
  std::string s;

  std::vector<std::string> split_strings;
  while (std::getline(ss, s, '_')) split_strings.push_back(s);

  try
  {
    return std::stoi(split_strings.back());
  }
  catch (const std::invalid_argument& e)
  {
    FOUR_C_THROW(
        "You provided the wrong format. The integer number must be "
        "at the very last position of the name, separated by an underscore. "
        "The correct format is:\n"
        "\"<prefix_name>_<number>\"");
  }
  exit(EXIT_FAILURE);
}

FOUR_C_NAMESPACE_CLOSE
