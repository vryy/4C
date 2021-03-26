/*----------------------------------------------------------------------*/
/*! \file
\brief result testing functionality for scalar-structure interaction problems

\level 2


*/
/*----------------------------------------------------------------------*/
#include "ssi_resulttest.H"

#include "ssi_monolithic.H"
#include "ssi_str_model_evaluator_monolithic.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_linedefinition.H"

#include "../drt_scatra/scatra_timint_implicit.H"

#include "../linalg/linalg_solver.H"

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
SSI::SSIResultTest::SSIResultTest(const Teuchos::RCP<const SSI::SSIBase> ssi_base)
    : DRT::ResultTest("SSI"), ssi_base_(ssi_base)
{
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double SSI::SSIResultTest::ResultNode(const std::string& quantity, DRT::Node* node) const
{
  // initialize variable for result
  double result(0.0);

  // extract result
  if (!quantity.compare(0, 6, "stress"))
  {
    // extract nodal stresses
    const Epetra_MultiVector& stresses = dynamic_cast<const STR::MODELEVALUATOR::MonolithicSSI&>(
        ssi_base_->StructureField()->ModelEvaluator(INPAR::STR::model_monolithic_coupling))
                                             .Stresses();

    // extract local node ID
    const int lid = stresses.Map().LID(node->Id());
    if (lid < 0) dserror("Invalid node ID!");

    // extract desired value
    if (quantity == "stress_xx")
      result = (*stresses(0))[lid];
    else if (quantity == "stress_yy")
      result = (*stresses(1))[lid];
    else if (quantity == "stress_zz")
      result = (*stresses(2))[lid];
    else if (quantity == "stress_xy")
      result = (*stresses(3))[lid];
    else if (quantity == "stress_yz")
      result = (*stresses(4))[lid];
    else if (quantity == "stress_xz")
      result = (*stresses(5))[lid];
    else
      dserror("Invalid stress specification!");
  }

  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
}

/*----------------------------------------------------------------------*
 *----------------------------------------------------------------------*/
double SSI::SSIResultTest::ResultSpecial(const std::string& quantity) const
{
  // initialize variable for result
  double result(0.0);

  // number of outer coupling iterations (partitioned SSI) or Newton-Raphson iterations (monolithic
  // SSI) in last time step
  if (quantity == "numiterlastnonlinearsolve")
  {
    result = static_cast<double>(ssi_base_->IterationCount());
  }

  // number of iterations performed by linear solver during last Newton-Raphson iteration
  // (monolithic SSI only)
  else if (quantity == "numiterlastlinearsolve")
  {
    // safety check
    if (SSIMono().Solver().Params().get("solver", "none") != "aztec")
    {
      dserror(
          "Must have Aztec solver for result test involving number of solver iterations during "
          "last Newton-Raphson iteration!");
    }
    result = static_cast<double>(SSIMono().Solver().getNumIters());
  }

  // test total number of time steps
  else if (!quantity.compare(0, 7, "numstep"))
  {
    result = static_cast<double>(ssi_base_->Step());
  }
  // catch unknown quantity strings
  else
  {
    dserror(
        "Quantity '%s' not supported by result testing functionality for scalar-structure "
        "interaction!",
        quantity.c_str());
  }

  return result;
}

/*---------------------------------------------------------------------------------*
 *---------------------------------------------------------------------------------*/
const SSI::SSIMono& SSI::SSIResultTest::SSIMono() const
{
  const auto* const ssi_mono = dynamic_cast<const SSI::SSIMono* const>(ssi_base_.get());
  if (ssi_mono == nullptr)
    dserror("Couldn't access time integrator for monolithic scalar-structure interaction!");
  return *ssi_mono;
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
void SSI::SSIResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // determine discretization
  std::string dis;
  res.ExtractString("DIS", dis);
  const DRT::Discretization* discretization(nullptr);
  if (dis == ssi_base_->ScaTraField()->Discretization()->Name())
    discretization = ssi_base_->ScaTraField()->Discretization().get();
  else if (dis == ssi_base_->StructureField()->Discretization()->Name())
    discretization = ssi_base_->StructureField()->Discretization().get();
  else
    dserror("Invalid discretization name!");

  // extract global node ID
  int node(-1);
  res.ExtractInt("NODE", node);
  --node;

  // check global node ID
  int isowner(discretization->HaveGlobalNode(node) and
              discretization->gNode(node)->Owner() == discretization->Comm().MyPID());
  int hasowner(-1);
  discretization->Comm().SumAll(&isowner, &hasowner, 1);
  if (hasowner != 1) dserror("Node must have exactly one owner!");

  // result test performed by node owner
  if (isowner)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY", quantity);

    // perform result test
    nerr += CompareValues(ResultNode(quantity, discretization->gNode(node)), "NODE", res);
    ++test_count;
  }
}

/*-------------------------------------------------------------------------------------*
 *-------------------------------------------------------------------------------------*/
void SSI::SSIResultTest::TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // make sure that quantity is tested only by one processor
  if (ssi_base_->Comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY", quantity);

    // get result to be tested
    const double result = ResultSpecial(quantity);

    // compare values
    const int err = CompareValues(result, "SPECIAL", res);
    nerr += err;
    ++test_count;
  }
}
