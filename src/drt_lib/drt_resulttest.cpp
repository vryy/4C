/*----------------------------------------------------------------------*/
/*!
\file drt_resulttest.cpp

\brief general result test framework

<pre>
Maintainer: Ulrich Kuettler
            kuettler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/Members/kuettler
            089 - 289-15238
</pre>
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


DRT::ResultTest::ResultTest(const std::string name)
  : myname_(name)
{
}

DRT::ResultTest::~ResultTest()
{
}

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
int DRT::ResultTest::CompareValues(double actresult, DRT::INPUT::LineDefinition& res)
{
  int node;
  res.ExtractInt("NODE",node);
  std::string quantity;
  res.ExtractString("QUANTITY",quantity);
  double givenresult;
  res.ExtractDouble("VALUE",givenresult);
  double tolerance;
  res.ExtractDouble("TOLERANCE",tolerance);
  // name is an optional input argument!
  std::string name = "";
  if (res.HaveNamed("NAME"))
    res.ExtractString("NAME",name);

  // return value (0 if results are correct, 1 if results are not correct)
  int ret = 0;

  // write to error file
  FILE *err = DRT::Problem::Instance()->ErrorFile()->Handle();
  if (err != NULL)
  {
    fprintf(err,"actual = %.17e, given = %.17e, diff = %.17e\n",
            actresult, givenresult, actresult-givenresult);
  }

  // prepare string stream 'msghead' containing general information on the current test
  std::stringstream msghead;
  msghead << std::left << std::setw(9) << myname_
          << ": "
          << std::left << std::setw(8) << quantity.c_str();

  if (name != "")
    msghead << "(" << name << ")";

  msghead << " at node "
          << std::right << std::setw(3)<< node;

  // write something to screen depending if the result check was ok or not
  if (!(fabs(fabs(actresult-givenresult)-fabs(actresult-givenresult)) < tolerance) )
  {
    // Result is 'not a number'
    IO::cout  << msghead.str()
              << " is NAN!\n";
    ret = 1;
  }
  else if (fabs(actresult-givenresult) > tolerance)
  {
    // Result is wrong
    IO::cout  << msghead.str()
              << " is WRONG --> actresult="
              << std::setw(24) << std::setprecision(17) << std::scientific << actresult
              << ", givenresult="
              << std::setw(24) << std::setprecision(17) << std::scientific << givenresult
              << "\n";

    ret = 1;
  }
  else
  {
    // Result is correct
    IO::cout  << msghead.str()
              << " is CORRECT\n";
  }

  return ret;
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
bool DRT::ResultTest::Match(DRT::INPUT::LineDefinition& res)
{
  return res.HaveNamed(myname_);
}


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
  int nerr = 0;                     // number of tests with errors
  int test_count = 0;               // number of tests performed
  const int size = results_.size(); // total number of tests

  if (comm.MyPID()==0)
    IO::cout << "\nChecking results of "<< size <<" tests:\n";

  for (unsigned i=0; i<results_.size(); ++i)
  {
    DRT::INPUT::LineDefinition& res = *results_[i];

    for (unsigned j=0; j<fieldtest_.size(); ++j)
    {
      if (fieldtest_[j]->Match(res))
      {
        if (res.HaveNamed("ELEMENT"))
          fieldtest_[j]->TestElement(res,nerr,test_count);
        else if (res.HaveNamed("NODE"))
          fieldtest_[j]->TestNode(res,nerr,test_count);
        else
          fieldtest_[j]->TestSpecial(res,nerr,test_count);
      }
    }
  }

  // determine the total number of errors
  int numerr;
  comm.SumAll(&nerr,&numerr,1);

  if (numerr > 0)
  {
    dserror("Result check failed with %d errors out of %d tests", numerr, size);
  }
  else
  {
    FILE *err = DRT::Problem::Instance()->ErrorFile()->Handle();
    if (err != NULL)
      fprintf(err,"===========================================\n");
  }

  /* test_count == -1 means we had a special test routine. It's thus
   * illegal to use both a special routine and single tests. But who
   * wants that? */
  int count;
  if (test_count > -1)
  {
    comm.SumAll(&test_count,&count,1);

    /* It's indeed possible to count more tests than expected if you
     * happen to test values of a boundary element. We don't want this
     * dserror to go off in that case. */
    if (count < size)
    {
      dserror("expected %d tests but performed %d", size, count);
    }
  }

  if (comm.MyPID()==0)
    IO::cout << "\n" GREEN_LIGHT "OK (" << count << ")" END_COLOR "\n";
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
Teuchos::RCP<DRT::INPUT::Lines> DRT::ResultTestManager::ValidResultLines()
{
  DRT::INPUT::LineDefinition structure;
  structure
    .AddTag("STRUCTURE")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition fluid;
  fluid
    .AddTag("FLUID")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition ale;
  ale
    .AddTag("ALE")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition thermal;
  thermal
    .AddTag("THERMAL")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition scatra;
  scatra
    .AddTag("SCATRA")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition red_airway;
  red_airway
    .AddTag("RED_AIRWAY")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition art_net;
  art_net
    .AddTag("ARTNET")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition fld_adj;
  fld_adj
    .AddTag("ADJOINT")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition opti;
  opti
    .AddTag("OPTI")
    .AddNamedString("DIS")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  DRT::INPUT::LineDefinition fsi;
  fsi
    .AddTag("FSI")
    .AddNamedInt("NODE")
    .AddNamedString("QUANTITY")
    .AddNamedDouble("VALUE")
    .AddNamedDouble("TOLERANCE")
    .AddOptionalNamedString("NAME")
    ;

  Teuchos::RCP<DRT::INPUT::Lines> lines = Teuchos::rcp(new DRT::INPUT::Lines("RESULT DESCRIPTION"));
  lines->Add(structure);
  lines->Add(fluid);
  lines->Add(ale);
  lines->Add(thermal);
  lines->Add(scatra);
  lines->Add(red_airway);
  lines->Add(art_net);
  lines->Add(fld_adj);
  lines->Add(opti);
  lines->Add(fsi);

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
