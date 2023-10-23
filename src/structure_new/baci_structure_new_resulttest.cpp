/*-----------------------------------------------------------*/
/*! \file
\brief Structure specific result test class


\level 3

*/
/*-----------------------------------------------------------*/


#include "baci_structure_new_resulttest.H"

#include "baci_io_linedefinition.H"
#include "baci_lib_discret.H"
#include "baci_lib_globalproblem.H"
#include "baci_linalg_fixedsizematrix_voigt_notation.H"
#include "baci_structure_new_model_evaluator_data.H"
#include "baci_structure_new_timint_base.H"

#include <Epetra_MultiVector.h>
#include <Teuchos_RCPDecl.hpp>

#include <optional>
#include <string>


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
      dserror("You tried to test %s on a proc that does not own node %i.",
          name_and_component.name.c_str(), node_id);
    }

    return data[name_and_component.component][local_id];
  }
}  // namespace

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
STR::ResultTest::ResultTest()
    : DRT::ResultTest("STRUCTURE"),
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
void STR::ResultTest::Init(
    const STR::TIMINT::BaseDataGlobalState& gstate, const STR::MODELEVALUATOR::Data& data)
{
  issetup_ = false;

  disn_ = gstate.GetDisN();
  veln_ = gstate.GetVelN();
  accn_ = gstate.GetAccN();
  reactn_ = gstate.GetFreactN();
  gstate_ = Teuchos::rcpFromRef(gstate);
  data_ = Teuchos::rcpFromRef(data);
  strudisc_ = gstate.GetDiscret();

  if (DRT::Problem::Instance()->GetProblemType() == ProblemType::struct_ale and
      (DRT::Problem::Instance()->WearParams()).get<double>("WEARCOEFF") > 0.0)
    dserror("Material displ. are not yet considered!");
  else
    dismatn_ = Teuchos::null;

  isinit_ = true;
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ResultTest::Setup()
{
  CheckInit();
  // currently unused
  issetup_ = true;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::ResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  CheckInitSetup();

  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != strudisc_->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(strudisc_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  strudisc_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", node + 1, strudisc_->Name().c_str());
  }
  else
  {
    if (strudisc_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = strudisc_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != strudisc_->Comm().MyPID()) return;

      std::string position;
      res.ExtractString("QUANTITY", position);
      bool unknownpos = true;  // make sure the result value std::string can be handled
      double result = 0.0;     // will hold the actual result of run

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
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
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
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
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
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
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
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
                actnode->Id());
          result = (*accn_)[lid];
        }
      }

      // test nodal stresses
      if (position.rfind("stress", 0) == 0)
      {
        result = GetNodalStressComponent(position, node);
        unknownpos = false;
      }

      // test for any postprocessed gauss point data
      if (data_->GetGaussPointDataOutputManagerPtr() != Teuchos::null)
      {
        std::optional<QuantityNameAndComponent> name_and_component =
            GetGaussPointDataNameAndComponent(
                position, data_->GetGaussPointDataOutputManagerPtr()->GetQuantities());
        if (name_and_component.has_value())
        {
          result = GetGaussPointDataValue(*name_and_component, node,
              data_->GetGaussPointDataOutputManagerPtr()->GetNodalData());
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
            dserror("You tried to test %s on nonexistent dof %d on node %d", position.c_str(), idx,
                actnode->Id());
          result = (*reactn_)[lid];
        }
      }

      // catch position std::strings, which are not handled by structure result test
      if (unknownpos) dserror("Quantity '%s' not supported in structure testing", position.c_str());

      // compare values
      const int err = CompareValues(result, "NODE", res);
      nerr += err;
      test_count++;
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
void STR::ResultTest::TestSpecial(
    DRT::INPUT::LineDefinition& res, int& nerr, int& test_count, int& uneval_test_count)
{
  CheckInitSetup();

  if (strudisc_->Comm().MyPID() != 0) return;

  std::string quantity;
  res.ExtractString("QUANTITY", quantity);

  Status special_status = Status::unevaluated;
  const double result = GetSpecialResult(quantity, special_status);
  switch (special_status)
  {
    case Status::evaluated:
    {
      nerr += CompareValues(result, "SPECIAL", res);
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
      dserror("What shall be done for this Status type? (enum=%d)", special_status);
      exit(EXIT_FAILURE);
    }
  }
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::ResultTest::GetSpecialResult(const std::string& quantity, Status& special_status) const
{
  if (quantity.find("num_iter_step_") != quantity.npos)
  {
    return GetNlnIterationNumber(quantity, special_status);
  }
  else if (quantity.find("lin_iter_step_") != quantity.npos)
  {
    return GetLastLinIterationNumber(quantity, special_status);
  }
  else if (quantity == "internal_energy" or quantity == "kinetic_energy" or
           quantity == "total_energy" or quantity == "beam_contact_penalty_potential" or
           quantity == "beam_interaction_potential" or
           quantity == "beam_to_beam_link_internal_energy" or
           quantity == "beam_to_beam_link_kinetic_energy" or
           quantity == "beam_to_sphere_link_internal_energy" or
           quantity == "beam_to_sphere_link_kinetic_energy")
  {
    return GetEnergy(quantity, special_status);
  }
  else
    dserror(
        "Quantity '%s' not supported by special result testing functionality "
        "for structure field!",
        quantity.c_str());

  exit(EXIT_FAILURE);
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
int STR::ResultTest::GetLastLinIterationNumber(
    const std::string& quantity, Status& special_status) const
{
  const int stepn = GetIntegerNumberAtLastPositionOfName(quantity);

  const int restart = DRT::Problem::Instance()->Restart();
  if (stepn <= restart) return -1;

  special_status = Status::evaluated;
  return static_cast<double>(gstate_->GetLastLinIterationNumber(stepn));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::ResultTest::GetNlnIterationNumber(
    const std::string& quantity, Status& special_status) const
{
  const int stepn = GetIntegerNumberAtLastPositionOfName(quantity);

  const int restart = DRT::Problem::Instance()->Restart();
  if (stepn <= restart) return -1.0;

  special_status = Status::evaluated;
  return static_cast<double>(gstate_->GetNlnIterationNumber(stepn));
}

/*----------------------------------------------------------------------------*
 *----------------------------------------------------------------------------*/
double STR::ResultTest::GetEnergy(const std::string& quantity, Status& special_status) const
{
  special_status = Status::evaluated;
  return data_->GetEnergyData(quantity);
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
    dserror(
        "You provided the wrong format. The integer number must be "
        "at the very last position of the name, separated by an underscore. "
        "The correct format is:\n"
        "\"<prefix_name>_<number>\"");
  }
  exit(EXIT_FAILURE);
}

double STR::ResultTest::GetNodalStressComponent(const std::string& label, int node_id) const
{
  int stress_voigt_index = -1;
  if (label == "stress_xx")
  {
    stress_voigt_index = CORE::LINALG::VOIGT::IndexMappings::SymToVoigt6(0, 0);
  }
  else if (label == "stress_yy")
  {
    stress_voigt_index = CORE::LINALG::VOIGT::IndexMappings::SymToVoigt6(1, 1);
  }
  else if (label == "stress_zz")
  {
    stress_voigt_index = CORE::LINALG::VOIGT::IndexMappings::SymToVoigt6(2, 2);
  }
  else if (label == "stress_xy")
  {
    stress_voigt_index = CORE::LINALG::VOIGT::IndexMappings::SymToVoigt6(0, 1);
  }
  else if (label == "stress_xz")
  {
    stress_voigt_index = CORE::LINALG::VOIGT::IndexMappings::SymToVoigt6(0, 2);
  }
  else if (label == "stress_yz")
  {
    stress_voigt_index = CORE::LINALG::VOIGT::IndexMappings::SymToVoigt6(1, 2);
  }

  if (stress_voigt_index < 0)
  {
    dserror(
        "You try to test an unknown stress component %s. Use one of [stress_xx, stress_yy, "
        "stress_zz, stress_xy, stress_xz, stress_yz]",
        label.c_str());
  }

  Teuchos::RCP<Epetra_MultiVector> nodalStressData = data_->GetStressDataNodePostprocessed();

  if (Teuchos::is_null(nodalStressData))
  {
    dserror(
        "It looks like you don't write stresses. You have to specify the stress type in "
        "IO->STRUCT_STRESS");
  }

  int local_id = nodalStressData->Map().LID(node_id);

  if (local_id < 0)
  {
    dserror("You tried to test %s on a proc that does not own node %i.", label.c_str(), node_id);
  }

  return (*nodalStressData)[stress_voigt_index][local_id];
}
