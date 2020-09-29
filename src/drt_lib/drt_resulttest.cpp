/*----------------------------------------------------------------------*/
/*! \file

\brief Implementation of general result test framework

\level 1


*/
/*----------------------------------------------------------------------*/

#include "drt_resulttest.H"
#include "drt_dserror.H"
#include "../drt_lib/drt_colors.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_io/io_control.H"
#include "../drt_io/io_pstream.H"
#include "drt_linedefinition.H"
#include "drt_inputreader.H"


DRT::ResultTest::ResultTest(const std::string name) : myname_(name) {}

DRT::ResultTest::~ResultTest() {}

void DRT::ResultTest::TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  dserror("no element test available");
}

void DRT::ResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  dserror("no node test available");
}

void DRT::ResultTest::TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  dserror("no special case test available");
}

void DRT::ResultTest::TestSpecial(
    DRT::INPUT::LineDefinition& res, int& nerr, int& test_count, int& unevaluated_test_count)
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
int DRT::ResultTest::CompareValues(
    double actresult, std::string type, DRT::INPUT::LineDefinition& res)
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
  if (tolerance <= 0.) dserror("Tolerance for result test must be strictly positive!");
  // name is an optional input argument!
  std::string name = "";
  if (res.HaveNamed("NAME")) res.ExtractString("NAME", name);

  // return value (0 if results are correct, 1 if results are not correct)
  int ret = 0;

  // write to error file
  FILE* err = DRT::Problem::Instance()->ErrorFile()->Handle();
  if (err != NULL)
  {
    fprintf(err, "actual = %.17e, given = %.17e, diff = %.17e\n", actresult, givenresult,
        actresult - givenresult);
  }

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
  if (!(fabs(fabs(actresult - givenresult) - fabs(actresult - givenresult)) < tolerance))
  {
    // Result is 'not a number'
    std::cout << msghead.str() << "\t is NAN!\n";
    ret = 1;
  }
  else if (fabs(actresult - givenresult) > tolerance)
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
bool DRT::ResultTest::Match(DRT::INPUT::LineDefinition& res) { return res.HaveNamed(myname_); }


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

  for (unsigned i = 0; i < results_.size(); ++i)
  {
    DRT::INPUT::LineDefinition& res = *results_[i];

    for (unsigned j = 0; j < fieldtest_.size(); ++j)
    {
      if (fieldtest_[j]->Match(res))
      {
        if (res.HaveNamed("ELEMENT"))
          fieldtest_[j]->TestElement(res, nerr, test_count);
        else if (res.HaveNamed("NODE"))
          fieldtest_[j]->TestNode(res, nerr, test_count);
        else
          fieldtest_[j]->TestSpecial(res, nerr, test_count, uneval_test_count);
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
    dserror("Result check failed with %d errors out of %d tests", numerr, size);
  }
  else
  {
    FILE* err = DRT::Problem::Instance()->ErrorFile()->Handle();
    if (err != NULL) fprintf(err, "===========================================\n");
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
     * dserror to go off in that case. */
    if (count < size)
    {
      dserror("expected %d tests but performed %d", size, count);
    }
  }

  if (comm.MyPID() == 0) IO::cout << "\n" GREEN_LIGHT "OK (" << count << ")" END_COLOR "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::INPUT::Lines> DRT::ResultTestManager::ValidResultLines()
{
  DRT::INPUT::LineDefinition structure;
  structure.AddTag("STRUCTURE")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition structure_special;
  structure_special.AddTag("STRUCTURE")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition fluid_node;
  fluid_node.AddTag("FLUID")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition fluid_ele;
  fluid_ele.AddTag("FLUID")
      .AddNamedString("DIS")
      .AddNamedInt("ELEMENT")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition xfluid_node;
  xfluid_node.AddTag("XFLUID")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition ale;
  ale.AddTag("ALE")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition thermal;
  thermal.AddTag("THERMAL")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition lubrication;
  lubrication.AddTag("LUBRICATION")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition porofluidmultiphase_node;
  porofluidmultiphase_node.AddTag("POROFLUIDMULTIPHASE")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition porofluidmultiphase_ele;
  porofluidmultiphase_ele.AddTag("POROFLUIDMULTIPHASE")
      .AddNamedString("DIS")
      .AddNamedInt("ELEMENT")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition porofluidmultiphase_special;
  porofluidmultiphase_special.AddTag("POROFLUIDMULTIPHASE")
      .AddNamedString("DIS")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition scatra;
  scatra.AddTag("SCATRA")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition scatra_special;
  scatra_special.AddTag("SCATRA")
      .AddNamedString("DIS")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition ssi;
  ssi.AddTag("SSI")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition ssi_special;
  ssi_special.AddTag("SSI")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE");

  DRT::INPUT::LineDefinition sti_special;
  sti_special.AddTag("STI")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE");

  DRT::INPUT::LineDefinition red_airway;
  red_airway.AddTag("RED_AIRWAY")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition red_airway_ele;
  red_airway_ele.AddTag("RED_AIRWAY")
      .AddNamedString("DIS")
      .AddNamedInt("ELEMENT")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition art_net_node;
  art_net_node.AddTag("ARTNET")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition art_net_ele;
  art_net_ele.AddTag("ARTNET")
      .AddNamedString("DIS")
      .AddNamedInt("ELEMENT")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition fld_adj;
  fld_adj.AddTag("ADJOINT")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition opti_node;
  opti_node.AddTag("OPTI")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition opti_ele;
  opti_ele.AddTag("OPTI")
      .AddNamedString("DIS")
      .AddNamedInt("ELEMENT")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition fsi_node;
  fsi_node.AddTag("FSI")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition fsi_special;
  fsi_special.AddTag("FSI")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition invana;
  invana.AddTag("INVANA")
      .AddNamedString("DIS")
      .AddNamedInt("ELEMENT")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition invana_special;
  invana_special.AddTag("INVANA")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition particle;
  particle.AddTag("PARTICLE")
      .AddNamedInt("ID")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE");

  DRT::INPUT::LineDefinition particlewall_node;
  particlewall_node.AddTag("PARTICLEWALL")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE");

  DRT::INPUT::LineDefinition particlewall_special;
  particlewall_special.AddTag("PARTICLEWALL")
      .AddNamedString("DIS")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE");

  DRT::INPUT::LineDefinition rigidbody;
  rigidbody.AddTag("RIGIDBODY")
      .AddNamedInt("ID")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE");

  DRT::INPUT::LineDefinition elemag;
  elemag.AddTag("ELECTROMAGNETIC")
      .AddNamedString("DIS")
      .AddNamedInt("NODE")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  DRT::INPUT::LineDefinition cardiovascular0d;
  cardiovascular0d.AddTag("CARDIOVASCULAR0D")
      .AddNamedString("DIS")
      .AddTag("SPECIAL")
      .AddNamedString("QUANTITY")
      .AddNamedDouble("VALUE")
      .AddNamedDouble("TOLERANCE")
      .AddOptionalNamedString("NAME");

  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("RESULT DESCRIPTION"));
  lines->Add(structure);
  lines->Add(structure_special);
  lines->Add(fluid_node);
  lines->Add(fluid_ele);
  lines->Add(xfluid_node);
  lines->Add(ale);
  lines->Add(thermal);
  lines->Add(lubrication);
  lines->Add(porofluidmultiphase_node);
  lines->Add(porofluidmultiphase_ele);
  lines->Add(porofluidmultiphase_special);
  lines->Add(scatra);
  lines->Add(scatra_special);
  lines->Add(ssi);
  lines->Add(ssi_special);
  lines->Add(sti_special);
  lines->Add(red_airway);
  lines->Add(red_airway_ele);
  lines->Add(art_net_node);
  lines->Add(art_net_ele);
  lines->Add(fld_adj);
  lines->Add(opti_node);
  lines->Add(opti_ele);
  lines->Add(fsi_node);
  lines->Add(fsi_special);
  lines->Add(invana);
  lines->Add(invana_special);
  lines->Add(particle);
  lines->Add(particlewall_node);
  lines->Add(particlewall_special);
  lines->Add(rigidbody);
  lines->Add(elemag);
  lines->Add(cardiovascular0d);

  return lines;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void DRT::ResultTestManager::ReadInput(DRT::INPUT::DatFileReader& reader)
{
  Teuchos::RCP<DRT::INPUT::Lines> lines = ValidResultLines();
  results_ = lines->Read(reader);
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void PrintResultDescrDatHeader()
{
  DRT::ResultTestManager resulttestmanager;
  Teuchos::RCP<DRT::INPUT::Lines> lines = resulttestmanager.ValidResultLines();

  lines->Print(std::cout);
}
