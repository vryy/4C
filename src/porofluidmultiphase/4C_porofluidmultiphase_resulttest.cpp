/*----------------------------------------------------------------------*/
/*! \file
 \brief result test for multiphase porous flow

   \level 3

 *----------------------------------------------------------------------*/


#include "4C_porofluidmultiphase_resulttest.hpp"

#include "4C_fem_discretization.hpp"
#include "4C_global_data.hpp"
#include "4C_io_linedefinition.hpp"
#include "4C_porofluidmultiphase_meshtying_strategy_base.hpp"
#include "4C_porofluidmultiphase_timint_implicit.hpp"

FOUR_C_NAMESPACE_OPEN


/*----------------------------------------------------------------------*
 | ctor                                                     vuong 08/16 |
 *----------------------------------------------------------------------*/
POROFLUIDMULTIPHASE::ResultTest::ResultTest(TimIntImpl& porotimint)
    : Core::UTILS::ResultTest("POROFLUIDMULTIPHASE"), porotimint_(porotimint)
{
  return;
}


/*----------------------------------------------------------------------*
 | test node                                                vuong 08/16 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::ResultTest::test_node(
    Input::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != porotimint_.discretization()->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(porotimint_.discretization()->HaveGlobalNode(node));
  int isnodeofanybody(0);
  porotimint_.discretization()->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    FOUR_C_THROW("Node %d does not belong to discretization %s", node + 1,
        porotimint_.discretization()->Name().c_str());
  }
  else
  {
    if (porotimint_.discretization()->HaveGlobalNode(node))
    {
      Core::Nodes::Node* actnode = porotimint_.discretization()->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != porotimint_.discretization()->Comm().MyPID()) return;

      // extract name of quantity to be tested
      std::string quantity;
      res.ExtractString("QUANTITY", quantity);

      // get result to be tested
      const double result = result_node(quantity, actnode);

      nerr += compare_values(result, "NODE", res);
      test_count++;
    }
  }

  return;
}

/*----------------------------------------------------------------------*
 | test element                                        kremheller 10/19 |
 *----------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::ResultTest::TestElement(
    Input::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);

  if (dis != porotimint_.discretization()->Name()) return;

  int element;
  res.ExtractInt("ELEMENT", element);
  element -= 1;

  int haveelement(porotimint_.discretization()->HaveGlobalElement(element));
  int iselementofanybody(0);
  porotimint_.discretization()->Comm().SumAll(&haveelement, &iselementofanybody, 1);

  if (iselementofanybody == 0)
  {
    FOUR_C_THROW("Element %d does not belong to discretization %s", element + 1,
        porotimint_.discretization()->Name().c_str());
  }
  else
  {
    if (porotimint_.discretization()->HaveGlobalElement(element))
    {
      const Core::Elements::Element* actelement = porotimint_.discretization()->gElement(element);

      // Here we are just interested in the elements that we own (i.e. a row element)!
      if (actelement->Owner() != porotimint_.discretization()->Comm().MyPID()) return;

      // extract name of quantity to be tested
      std::string quantity;
      res.ExtractString("QUANTITY", quantity);

      // get result to be tested
      const double result = result_element(quantity, actelement);

      nerr += compare_values(result, "ELEMENT", res);
      test_count++;
    }
  }

  return;
}


/*----------------------------------------------------------------------*
 | get nodal result to be tested                            vuong 08/16 |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::ResultTest::result_node(
    const std::string quantity, Core::Nodes::Node* node) const
{
  // initialize variable for result
  double result(0.);

  // extract row map from solution vector
  const Epetra_BlockMap& phinpmap = porotimint_.Phinp()->Map();

  // test result value of phi field
  if (quantity == "phi")
    result = (*porotimint_.Phinp())[phinpmap.LID(porotimint_.discretization()->Dof(0, node, 0))];

  // test result value for a system of scalars
  else if (!quantity.compare(0, 3, "phi"))
  {
    // read species ID
    std::string k_string = quantity.substr(3);
    char* locator(nullptr);
    int k = strtol(k_string.c_str(), &locator, 10) - 1;
    if (locator == k_string.c_str()) FOUR_C_THROW("Couldn't read species ID!");

    if (porotimint_.discretization()->NumDof(0, node) <= k)
      FOUR_C_THROW("Species ID is larger than number of DOFs of node!");

    // extract result
    result = (*porotimint_.Phinp())[phinpmap.LID(porotimint_.discretization()->Dof(0, node, k))];
  }

  // test result value of phi field
  else if (quantity == "pressure")
    result = (*porotimint_.Pressure())[phinpmap.LID(porotimint_.discretization()->Dof(0, node, 0))];

  // test result value for a system of scalars
  else if (!quantity.compare(0, 8, "pressure"))
  {
    // read species ID
    std::string k_string = quantity.substr(8);
    char* locator(nullptr);
    int k = strtol(k_string.c_str(), &locator, 10) - 1;
    if (locator == k_string.c_str()) FOUR_C_THROW("Couldn't read pressure ID!");

    if (porotimint_.discretization()->NumDof(0, node) <= k)
      FOUR_C_THROW("Pressure ID is larger than number of DOFs of node!");

    // extract result
    result = (*porotimint_.Pressure())[phinpmap.LID(porotimint_.discretization()->Dof(0, node, k))];
  }

  // test result value of phi field
  else if (quantity == "saturation")
    result =
        (*porotimint_.Saturation())[phinpmap.LID(porotimint_.discretization()->Dof(0, node, 0))];

  // test result value for a system of scalars
  else if (!quantity.compare(0, 10, "saturation"))
  {
    // read species ID
    std::string k_string = quantity.substr(10);
    char* locator(nullptr);
    int k = strtol(k_string.c_str(), &locator, 10) - 1;
    if (locator == k_string.c_str()) FOUR_C_THROW("Couldn't read saturation ID!");

    if (porotimint_.discretization()->NumDof(0, node) <= k)
      FOUR_C_THROW("Saturation ID is larger than number of DOFs of node!");

    // extract result
    result =
        (*porotimint_.Saturation())[phinpmap.LID(porotimint_.discretization()->Dof(0, node, k))];
  }

  // catch unknown quantity strings
  else
    FOUR_C_THROW("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}  // POROFLUIDMULTIPHASE::ResultTest::ResultNode

/*----------------------------------------------------------------------*
 | get element result to be tested                     kremheller 10/19 |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::ResultTest::result_element(
    const std::string quantity, const Core::Elements::Element* element) const
{
  // initialize variable for result
  double result(0.);

  if (quantity == "bloodvesselvolfrac")
  {
    result = (*porotimint_.MeshTyingStrategy()->blood_vessel_volume_fraction())
        [porotimint_.discretization()->ElementRowMap()->LID(element->Id())];
  }
  else if (!quantity.compare(0, 13, "phasevelocity"))
  {
    const int num_dim = Global::Problem::Instance()->NDim();
    // get phase ID
    // example: "phasevelocity3x" -> k = 2 (phase IDs start at index 0)
    std::string k_string = quantity.substr(13);
    char* locator(nullptr);
    auto idx_poro_dof = int(strtol(k_string.c_str(), &locator, 13) - 1);
    if (locator == k_string.c_str()) FOUR_C_THROW("Could not read phase ID in result test!");

    // get spatial dimension
    int idx_dim(-1);
    if (!quantity.compare(14, 15, "x"))
      idx_dim = 0;
    else if (!quantity.compare(14, 15, "y"))
      idx_dim = 1;
    else if (!quantity.compare(14, 15, "z"))
      idx_dim = 2;

    result = ((*porotimint_.PhaseVelocity())[idx_poro_dof * num_dim + idx_dim])
        [porotimint_.discretization()->ElementRowMap()->LID(element->Id())];
  }
  // catch unknown quantity strings
  else
    FOUR_C_THROW("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}

/*-------------------------------------------------------------------------------------*
 | test special quantity not associated with a particular element or node  vuong 08/16 |
 *-------------------------------------------------------------------------------------*/
void POROFLUIDMULTIPHASE::ResultTest::TestSpecial(
    Input::LineDefinition& res, int& nerr, int& test_count)
{
  // make sure that quantity is tested only once
  if (porotimint_.discretization()->Comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY", quantity);

    // get result to be tested
    const double result = result_special(quantity);

    // compare values
    const int err = compare_values(result, "SPECIAL", res);
    nerr += err;
    test_count++;
  }

  return;
}


/*----------------------------------------------------------------------*
 | get special result to be tested                          vuong 08/16 |
 *----------------------------------------------------------------------*/
double POROFLUIDMULTIPHASE::ResultTest::result_special(
    const std::string quantity  //! name of quantity to be tested
) const
{
  // initialize variable for result
  double result(0.);

  if (quantity == "numiterlastnewton") result = (double)porotimint_.IterNum();
  // result test of domain integrals
  else if (!quantity.compare(0, 22, "domain_integral_value_"))
  {
    // get the index of the value which should be checked
    std::string suffix = quantity.substr(22);
    int idx = -1;
    try
    {
      idx = std::stoi(suffix);
    }
    catch (const std::invalid_argument& e)
    {
      FOUR_C_THROW(
          "You provided the wrong format for output of domain_integral_values. The integer number "
          "must be at the very last position of the name, separated by an underscore.\n"
          "The correct format is: domain_integral_value_<number>");
    }

    // index should be in range [0, number_functions - 1]
    if (idx < 0 || idx >= porotimint_.num_domain_int_functions())
      FOUR_C_THROW("detected wrong index %i, index should be in range [0,%i]", idx,
          porotimint_.num_domain_int_functions() - 1);

    // return the result
    result = (*porotimint_.DomainIntValues())[idx];
  }
  // catch unknown quantity strings
  else
    FOUR_C_THROW("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}  // POROFLUIDMULTIPHASE::ResultTest::result_special

FOUR_C_NAMESPACE_CLOSE
