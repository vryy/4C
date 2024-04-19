/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of general result test framework

\level 1


*/
/*----------------------------------------------------------------------*/

#include "baci_lib_resulttest.hpp"

#include "baci_io_inputreader.hpp"
#include "baci_io_linedefinition.hpp"
#include "baci_io_pstream.hpp"
#include "baci_utils_exceptions.hpp"

#include <utility>

FOUR_C_NAMESPACE_OPEN


DRT::ResultTest::ResultTest(std::string name) : myname_(std::move(name)) {}

void DRT::ResultTest::TestElement(INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  FOUR_C_THROW("no element test available");
}

void DRT::ResultTest::TestNode(INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  FOUR_C_THROW("no node test available");
}

void DRT::ResultTest::TestSpecial(INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  FOUR_C_THROW("no special case test available");
}

void DRT::ResultTest::TestSpecial(
    INPUT::LineDefinition& res, int& nerr, int& test_count, int& unevaluated_test_count)
{
  TestSpecial(res, nerr, test_count);
}


/*----------------------------------------------------------------------*/
/*!
 \brief Compare \a actresult with \a givenresult and return 0 if they are
 considered to be equal.

 Compare \a actresult with \a givenresult and return 0 if they are
 considered to be equal.

 \param err        the file where to document both values
 \param res        the describtion of the expected result including name and tolerance

 \author uk
 \date 06/04
 */
/*----------------------------------------------------------------------*/
int DRT::ResultTest::CompareValues(double actresult, std::string type, INPUT::LineDefinition& res)
{
  int gid;

  if (type != "SPECIAL")
  {
    res.ExtractInt(type, gid);
  }
  std::string quantity;
  res.ExtractString("QUANTITY", quantity);
  double givenresult;
  res.ExtractDouble("VALUE", givenresult);
  double tolerance;
  res.ExtractDouble("TOLERANCE", tolerance);
  // safety check
  if (tolerance <= 0.) FOUR_C_THROW("Tolerance for result test must be strictly positive!");
  // name is an optional input argument!
  std::string name = "";
  if (res.HaveNamed("NAME")) res.ExtractString("NAME", name);

  // return value (0 if results are correct, 1 if results are not correct)
  int ret = 0;

  // prepare std::string stream 'msghead' containing general information on the current test
  std::stringstream msghead;
  msghead << std::left << std::setw(9) << myname_ << ": " << std::left << std::setw(8)
          << quantity.c_str();

  if (name != "") msghead << "(" << name << ")";

  if (type != "SPECIAL")
  {
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    msghead << " at " << type << " " << std::right << std::setw(3) << gid;
  }
  else
  {
    msghead << "\t";
  }

  // write something to screen depending if the result check was ok or not
  if (std::isnan(actresult))
  {
    // Result is 'not a number'
    std::cout << msghead.str() << "\t is NAN!\n";
    ret = 1;
  }
  else if (std::abs(actresult - givenresult) > tolerance)
  {
    // Result is wrong
    std::cout << msghead.str() << "\t is WRONG --> actresult=" << std::setw(24)
              << std::setprecision(17) << std::scientific << actresult
              << ", givenresult=" << std::setw(24) << givenresult << ", abs(diff)=" << std::setw(24)
              << std::abs(actresult - givenresult) << " >" << std::setw(24) << tolerance << "\n";

    ret = 1;
  }
  else
  {
    // Result is correct
    std::cout << msghead.str() << "\t is CORRECT"
              << ", abs(diff)=" << std::setw(24) << std::setprecision(17) << std::scientific
              << std::abs(actresult - givenresult) << " <" << std::setw(24) << tolerance << "\n";
  }

  return ret;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ResultTest::Match(INPUT::LineDefinition& res) { return res.HaveNamed(myname_); }


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ResultTestManager::AddFieldTest(Teuchos::RCP<ResultTest> test)
{
  fieldtest_.push_back(test);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ResultTestManager::TestAll(const Epetra_Comm& comm)
{
  int nerr = 0;                      // number of tests with errors
  int test_count = 0;                // number of tests performed
  int uneval_test_count = 0;         // number of unevaluated tests
  const int size = results_.size();  // total number of tests

  if (comm.MyPID() == 0) IO::cout << "\nChecking results of " << size << " tests:\n";

  for (auto& result : results_)
  {
    for (const auto& fieldtest : fieldtest_)
    {
      if (fieldtest->Match(result))
      {
        if (result.HaveNamed("ELEMENT"))
          fieldtest->TestElement(result, nerr, test_count);
        else if (result.HaveNamed("NODE"))
          fieldtest->TestNode(result, nerr, test_count);
        else
          fieldtest->TestSpecial(result, nerr, test_count, uneval_test_count);
      }
    }
  }

  // print number of unevaluated tests to screen
  int guneval_test_count = 0;
  comm.SumAll(&uneval_test_count, &guneval_test_count, 1);
  if (guneval_test_count > 0 and comm.MyPID() == 0)
    IO::cout << guneval_test_count << " tests stay unevaluated" << IO::endl;

  // determine the total number of errors
  int numerr;
  comm.SumAll(&nerr, &numerr, 1);

  if (numerr > 0)
  {
    FOUR_C_THROW("Result check failed with %d errors out of %d tests", numerr, size);
  }

  /* test_count == -1 means we had a special test routine. It's thus
   * illegal to use both a special routine and single tests. But who
   * wants that? */
  int count;
  if (test_count > -1)
  {
    int lcount = test_count + uneval_test_count;
    comm.SumAll(&lcount, &count, 1);

    /* It's indeed possible to count more tests than expected if you
     * happen to test values of a boundary element. We don't want this
     * FOUR_C_THROW to go off in that case. */
    if (count < size)
    {
      FOUR_C_THROW("expected %d tests but performed %d", size, count);
    }
  }

  if (comm.MyPID() == 0) IO::cout << "\nOK (" << count << ")\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
INPUT::Lines DRT::ResultTestManager::ValidResultLines()
{
  INPUT::LineDefinition structure = INPUT::LineDefinition::Builder()
                                        .AddTag("STRUCTURE")
                                        .AddNamedString("DIS")
                                        .AddNamedInt("NODE")
                                        .AddNamedString("QUANTITY")
                                        .AddNamedDouble("VALUE")
                                        .AddNamedDouble("TOLERANCE")
                                        .AddOptionalNamedString("NAME")
                                        .Build();

  INPUT::LineDefinition structure_special = INPUT::LineDefinition::Builder()
                                                .AddTag("STRUCTURE")
                                                .AddTag("SPECIAL")
                                                .AddNamedString("QUANTITY")
                                                .AddNamedDouble("VALUE")
                                                .AddNamedDouble("TOLERANCE")
                                                .AddOptionalNamedString("NAME")
                                                .Build();

  INPUT::LineDefinition fluid_node = INPUT::LineDefinition::Builder()
                                         .AddTag("FLUID")
                                         .AddNamedString("DIS")
                                         .AddNamedInt("NODE")
                                         .AddNamedString("QUANTITY")
                                         .AddNamedDouble("VALUE")
                                         .AddNamedDouble("TOLERANCE")
                                         .AddOptionalNamedString("NAME")
                                         .Build();

  INPUT::LineDefinition fluid_ele = INPUT::LineDefinition::Builder()
                                        .AddTag("FLUID")
                                        .AddNamedString("DIS")
                                        .AddNamedInt("ELEMENT")
                                        .AddNamedString("QUANTITY")
                                        .AddNamedDouble("VALUE")
                                        .AddNamedDouble("TOLERANCE")
                                        .AddOptionalNamedString("NAME")
                                        .Build();

  INPUT::LineDefinition xfluid_node = INPUT::LineDefinition::Builder()
                                          .AddTag("XFLUID")
                                          .AddNamedString("DIS")
                                          .AddNamedInt("NODE")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .AddOptionalNamedString("NAME")
                                          .Build();

  INPUT::LineDefinition ale = INPUT::LineDefinition::Builder()
                                  .AddTag("ALE")
                                  .AddNamedString("DIS")
                                  .AddNamedInt("NODE")
                                  .AddNamedString("QUANTITY")
                                  .AddNamedDouble("VALUE")
                                  .AddNamedDouble("TOLERANCE")
                                  .AddOptionalNamedString("NAME")
                                  .Build();

  INPUT::LineDefinition thermal = INPUT::LineDefinition::Builder()
                                      .AddTag("THERMAL")
                                      .AddNamedString("DIS")
                                      .AddNamedInt("NODE")
                                      .AddNamedString("QUANTITY")
                                      .AddNamedDouble("VALUE")
                                      .AddNamedDouble("TOLERANCE")
                                      .AddOptionalNamedString("NAME")
                                      .Build();

  INPUT::LineDefinition lubrication = INPUT::LineDefinition::Builder()
                                          .AddTag("LUBRICATION")
                                          .AddNamedString("DIS")
                                          .AddNamedInt("NODE")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .AddOptionalNamedString("NAME")
                                          .Build();

  INPUT::LineDefinition porofluidmultiphase_node = INPUT::LineDefinition::Builder()
                                                       .AddTag("POROFLUIDMULTIPHASE")
                                                       .AddNamedString("DIS")
                                                       .AddNamedInt("NODE")
                                                       .AddNamedString("QUANTITY")
                                                       .AddNamedDouble("VALUE")
                                                       .AddNamedDouble("TOLERANCE")
                                                       .AddOptionalNamedString("NAME")
                                                       .Build();

  INPUT::LineDefinition porofluidmultiphase_ele = INPUT::LineDefinition::Builder()
                                                      .AddTag("POROFLUIDMULTIPHASE")
                                                      .AddNamedString("DIS")
                                                      .AddNamedInt("ELEMENT")
                                                      .AddNamedString("QUANTITY")
                                                      .AddNamedDouble("VALUE")
                                                      .AddNamedDouble("TOLERANCE")
                                                      .AddOptionalNamedString("NAME")
                                                      .Build();

  INPUT::LineDefinition porofluidmultiphase_special = INPUT::LineDefinition::Builder()
                                                          .AddTag("POROFLUIDMULTIPHASE")
                                                          .AddNamedString("DIS")
                                                          .AddTag("SPECIAL")
                                                          .AddNamedString("QUANTITY")
                                                          .AddNamedDouble("VALUE")
                                                          .AddNamedDouble("TOLERANCE")
                                                          .AddOptionalNamedString("NAME")
                                                          .Build();

  INPUT::LineDefinition scatra = INPUT::LineDefinition::Builder()
                                     .AddTag("SCATRA")
                                     .AddNamedString("DIS")
                                     .AddNamedInt("NODE")
                                     .AddNamedString("QUANTITY")
                                     .AddNamedDouble("VALUE")
                                     .AddNamedDouble("TOLERANCE")
                                     .AddOptionalNamedString("NAME")
                                     .Build();

  INPUT::LineDefinition scatra_special = INPUT::LineDefinition::Builder()
                                             .AddTag("SCATRA")
                                             .AddNamedString("DIS")
                                             .AddTag("SPECIAL")
                                             .AddNamedString("QUANTITY")
                                             .AddNamedDouble("VALUE")
                                             .AddNamedDouble("TOLERANCE")
                                             .AddOptionalNamedString("NAME")
                                             .Build();

  INPUT::LineDefinition ssi = INPUT::LineDefinition::Builder()
                                  .AddTag("SSI")
                                  .AddNamedString("DIS")
                                  .AddNamedInt("NODE")
                                  .AddNamedString("QUANTITY")
                                  .AddNamedDouble("VALUE")
                                  .AddNamedDouble("TOLERANCE")
                                  .AddOptionalNamedString("NAME")
                                  .Build();

  INPUT::LineDefinition ssi_special = INPUT::LineDefinition::Builder()
                                          .AddTag("SSI")
                                          .AddTag("SPECIAL")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .Build();

  INPUT::LineDefinition ssti_special = INPUT::LineDefinition::Builder()
                                           .AddTag("SSTI")
                                           .AddTag("SPECIAL")
                                           .AddNamedString("QUANTITY")
                                           .AddNamedDouble("VALUE")
                                           .AddNamedDouble("TOLERANCE")
                                           .Build();

  INPUT::LineDefinition sti_special = INPUT::LineDefinition::Builder()
                                          .AddTag("STI")
                                          .AddTag("SPECIAL")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .Build();

  INPUT::LineDefinition red_airway = INPUT::LineDefinition::Builder()
                                         .AddTag("RED_AIRWAY")
                                         .AddNamedString("DIS")
                                         .AddNamedInt("NODE")
                                         .AddNamedString("QUANTITY")
                                         .AddNamedDouble("VALUE")
                                         .AddNamedDouble("TOLERANCE")
                                         .AddOptionalNamedString("NAME")
                                         .Build();

  INPUT::LineDefinition red_airway_ele = INPUT::LineDefinition::Builder()
                                             .AddTag("RED_AIRWAY")
                                             .AddNamedString("DIS")
                                             .AddNamedInt("ELEMENT")
                                             .AddNamedString("QUANTITY")
                                             .AddNamedDouble("VALUE")
                                             .AddNamedDouble("TOLERANCE")
                                             .AddOptionalNamedString("NAME")
                                             .Build();

  INPUT::LineDefinition art_net_node = INPUT::LineDefinition::Builder()
                                           .AddTag("ARTNET")
                                           .AddNamedString("DIS")
                                           .AddNamedInt("NODE")
                                           .AddNamedString("QUANTITY")
                                           .AddNamedDouble("VALUE")
                                           .AddNamedDouble("TOLERANCE")
                                           .AddOptionalNamedString("NAME")
                                           .Build();

  INPUT::LineDefinition art_net_ele = INPUT::LineDefinition::Builder()
                                          .AddTag("ARTNET")
                                          .AddNamedString("DIS")
                                          .AddNamedInt("ELEMENT")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .AddOptionalNamedString("NAME")
                                          .Build();

  INPUT::LineDefinition fld_adj = INPUT::LineDefinition::Builder()
                                      .AddTag("ADJOINT")
                                      .AddNamedString("DIS")
                                      .AddNamedInt("NODE")
                                      .AddNamedString("QUANTITY")
                                      .AddNamedDouble("VALUE")
                                      .AddNamedDouble("TOLERANCE")
                                      .AddOptionalNamedString("NAME")
                                      .Build();

  INPUT::LineDefinition opti_node = INPUT::LineDefinition::Builder()
                                        .AddTag("OPTI")
                                        .AddNamedString("DIS")
                                        .AddNamedInt("NODE")
                                        .AddNamedString("QUANTITY")
                                        .AddNamedDouble("VALUE")
                                        .AddNamedDouble("TOLERANCE")
                                        .AddOptionalNamedString("NAME")
                                        .Build();

  INPUT::LineDefinition opti_ele = INPUT::LineDefinition::Builder()
                                       .AddTag("OPTI")
                                       .AddNamedString("DIS")
                                       .AddNamedInt("ELEMENT")
                                       .AddNamedString("QUANTITY")
                                       .AddNamedDouble("VALUE")
                                       .AddNamedDouble("TOLERANCE")
                                       .AddOptionalNamedString("NAME")
                                       .Build();

  INPUT::LineDefinition fsi_node = INPUT::LineDefinition::Builder()
                                       .AddTag("FSI")
                                       .AddNamedInt("NODE")
                                       .AddNamedString("QUANTITY")
                                       .AddNamedDouble("VALUE")
                                       .AddNamedDouble("TOLERANCE")
                                       .AddOptionalNamedString("NAME")
                                       .Build();

  INPUT::LineDefinition fsi_special = INPUT::LineDefinition::Builder()
                                          .AddTag("FSI")
                                          .AddTag("SPECIAL")
                                          .AddNamedString("QUANTITY")
                                          .AddNamedDouble("VALUE")
                                          .AddNamedDouble("TOLERANCE")
                                          .AddOptionalNamedString("NAME")
                                          .Build();

  INPUT::LineDefinition particle = INPUT::LineDefinition::Builder()
                                       .AddTag("PARTICLE")
                                       .AddNamedInt("ID")
                                       .AddNamedString("QUANTITY")
                                       .AddNamedDouble("VALUE")
                                       .AddNamedDouble("TOLERANCE")
                                       .Build();

  INPUT::LineDefinition particlewall_node = INPUT::LineDefinition::Builder()
                                                .AddTag("PARTICLEWALL")
                                                .AddNamedString("DIS")
                                                .AddNamedInt("NODE")
                                                .AddNamedString("QUANTITY")
                                                .AddNamedDouble("VALUE")
                                                .AddNamedDouble("TOLERANCE")
                                                .Build();

  INPUT::LineDefinition particlewall_special = INPUT::LineDefinition::Builder()
                                                   .AddTag("PARTICLEWALL")
                                                   .AddNamedString("DIS")
                                                   .AddTag("SPECIAL")
                                                   .AddNamedString("QUANTITY")
                                                   .AddNamedDouble("VALUE")
                                                   .AddNamedDouble("TOLERANCE")
                                                   .Build();

  INPUT::LineDefinition rigidbody = INPUT::LineDefinition::Builder()
                                        .AddTag("RIGIDBODY")
                                        .AddNamedInt("ID")
                                        .AddNamedString("QUANTITY")
                                        .AddNamedDouble("VALUE")
                                        .AddNamedDouble("TOLERANCE")
                                        .Build();

  INPUT::LineDefinition elemag = INPUT::LineDefinition::Builder()
                                     .AddTag("ELECTROMAGNETIC")
                                     .AddNamedString("DIS")
                                     .AddNamedInt("NODE")
                                     .AddNamedString("QUANTITY")
                                     .AddNamedDouble("VALUE")
                                     .AddNamedDouble("TOLERANCE")
                                     .AddOptionalNamedString("NAME")
                                     .Build();

  INPUT::LineDefinition cardiovascular0d = INPUT::LineDefinition::Builder()
                                               .AddTag("CARDIOVASCULAR0D")
                                               .AddNamedString("DIS")
                                               .AddTag("SPECIAL")
                                               .AddNamedString("QUANTITY")
                                               .AddNamedDouble("VALUE")
                                               .AddNamedDouble("TOLERANCE")
                                               .AddOptionalNamedString("NAME")
                                               .Build();

  INPUT::Lines lines("RESULT DESCRIPTION",
      "The result of the simulation with respect to specific quantities at concrete points "
      "can be tested against particular values with a given tolerance.");
  lines.Add(structure);
  lines.Add(structure_special);
  lines.Add(fluid_node);
  lines.Add(fluid_ele);
  lines.Add(xfluid_node);
  lines.Add(ale);
  lines.Add(thermal);
  lines.Add(lubrication);
  lines.Add(porofluidmultiphase_node);
  lines.Add(porofluidmultiphase_ele);
  lines.Add(porofluidmultiphase_special);
  lines.Add(scatra);
  lines.Add(scatra_special);
  lines.Add(ssi);
  lines.Add(ssi_special);
  lines.Add(ssti_special);
  lines.Add(sti_special);
  lines.Add(red_airway);
  lines.Add(red_airway_ele);
  lines.Add(art_net_node);
  lines.Add(art_net_ele);
  lines.Add(fld_adj);
  lines.Add(opti_node);
  lines.Add(opti_ele);
  lines.Add(fsi_node);
  lines.Add(fsi_special);
  lines.Add(particle);
  lines.Add(particlewall_node);
  lines.Add(particlewall_special);
  lines.Add(rigidbody);
  lines.Add(elemag);
  lines.Add(cardiovascular0d);

  return lines;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ResultTestManager::ReadInput(INPUT::DatFileReader& reader)
{
  INPUT::Lines lines = ValidResultLines();
  results_ = lines.Read(reader);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintResultDescrDatHeader()
{
  DRT::ResultTestManager resulttestmanager;
  INPUT::Lines lines = resulttestmanager.ValidResultLines();

  lines.Print(std::cout);
}

FOUR_C_NAMESPACE_CLOSE
