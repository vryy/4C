/*----------------------------------------------------------------------*/
/*!
\file ssi_monolithic_resulttest.cpp

\brief result testing functionality for monolithic scalar-structure interaction problems

\level 2

<pre>
\maintainer Rui Fang
            fang@lnm.mw.tum.de
            http://www.lnm.mw.tum.de/
            089-289-15251
</pre>
*/
/*----------------------------------------------------------------------*/
#include "ssi_monolithic_resulttest.H"

#include "ssi_monolithic.H"
#include "ssi_str_model_evaluator_monolithic.H"

#include "../drt_adapter/ad_str_ssiwrapper.H"
#include "../drt_adapter/adapter_scatra_base_algorithm.H"

#include "../drt_lib/drt_linedefinition.H"

#include "../drt_scatra/scatra_timint_implicit.H"

/*----------------------------------------------------------------------*
 | constructor                                               fang 11/17 |
 *----------------------------------------------------------------------*/
SSI::SSI_Mono_ResultTest::SSI_Mono_ResultTest(
    const Teuchos::RCP<const SSI::SSI_Mono>   ssi_mono   //!< time integrator for monolithic scalar-structure interaction
    )
    // call base class constructor
  : DRT::ResultTest("SSI"),

    // store pointer to time integrator for monolithic scalar-structure interaction
    ssi_mono_(ssi_mono)
{
  return;
}


/*----------------------------------------------------------------------*
 | get nodal result to be tested                             fang 12/17 |
 *----------------------------------------------------------------------*/
double SSI::SSI_Mono_ResultTest::ResultNode(
    const std::string   quantity,   //!< name of quantity to be tested
    DRT::Node*          node        //!< node carrying the result to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // extract result
  if(!quantity.compare(0,6,"stress"))
  {
    // extract nodal stresses
    const Epetra_MultiVector& stresses = dynamic_cast<const STR::MODELEVALUATOR::MonolithicSSI&>(ssi_mono_->StructureField()->ModelEvaluator(INPAR::STR::model_monolithic_coupling)).Stresses();

    // extract local node ID
    const int lid = stresses.Map().LID(node->Id());
    if(lid < 0)
      dserror("Invalid node ID!");

    // extract desired value
    if(quantity == "stress_xx")
      result = (*stresses(0))[lid];
    else if(quantity == "stress_yy")
      result = (*stresses(1))[lid];
    else if(quantity == "stress_zz")
      result = (*stresses(2))[lid];
    else if(quantity == "stress_xy")
      result = (*stresses(3))[lid];
    else if(quantity == "stress_yz")
      result = (*stresses(4))[lid];
    else if(quantity == "stress_xz")
      result = (*stresses(5))[lid];
    else
      dserror("Invalid stress specification!");
  }

  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported in result test!", quantity.c_str());

  return result;
} // SSI::SSI_Mono_ResultTest::ResultNode


/*----------------------------------------------------------------------*
 | get special result to be tested                           fang 11/17 |
 *----------------------------------------------------------------------*/
double SSI::SSI_Mono_ResultTest::ResultSpecial(
    const std::string&   quantity   //!< name of quantity to be tested
    ) const
{
  // initialize variable for result
  double result(0.);

  // number of Newton-Raphson iterations in last time step
  if(quantity == "numiterlastnonlinearsolve")
    result = (double) ssi_mono_->Iter();

  // catch unknown quantity strings
  else
    dserror("Quantity '%s' not supported by result testing functionality for monolithic scalar-structure interaction!",quantity.c_str());

  return result;
} // SSI::SSI_Mono_ResultTest::ResultSpecial


/*-------------------------------------------------------------------------------------*
 | test quantity associated with a particular node                          fang 12/17 |
 *-------------------------------------------------------------------------------------*/
void SSI::SSI_Mono_ResultTest::TestNode(
    DRT::INPUT::LineDefinition&   res,         //!< input file line containing result test specification
    int&                          nerr,        //!< number of failed result tests
    int&                          test_count   //!< number of result tests
    )
{
  // determine discretization
  std::string dis;
  res.ExtractString("DIS",dis);
  const DRT::Discretization* discretization(NULL);
  if(dis == ssi_mono_->ScaTraField()->ScaTraField()->Discretization()->Name())
    discretization = ssi_mono_->ScaTraField()->ScaTraField()->Discretization().get();
  else if(dis == ssi_mono_->StructureField()->Discretization()->Name())
    discretization = ssi_mono_->StructureField()->Discretization().get();
  else
    dserror("Invalid discretization name!");

  // extract global node ID
  int node(-1);
  res.ExtractInt("NODE",node);
  --node;

  // check global node ID
  int isowner(discretization->HaveGlobalNode(node) and discretization->gNode(node)->Owner() == discretization->Comm().MyPID());
  int hasowner(-1);
  discretization->Comm().SumAll(&isowner,&hasowner,1);
  if(hasowner != 1)
    dserror("Node must have exactly one owner!");

  // result test performed by node owner
  if(isowner)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY",quantity);

    // perform result test
    nerr += CompareValues(ResultNode(quantity,discretization->gNode(node)),"NODE",res);
    ++test_count;
  }

  return;
} // SSI::SSI_Mono_ResultTest::TestNode


/*-------------------------------------------------------------------------------------*
 | test special quantity not associated with a particular element or node   fang 11/17 |
 *-------------------------------------------------------------------------------------*/
void SSI::SSI_Mono_ResultTest::TestSpecial(
    DRT::INPUT::LineDefinition&   res,         //!< input file line containing result test specification
    int&                          nerr,        //!< number of failed result tests
    int&                          test_count   //!< number of result tests
    )
{
  // make sure that quantity is tested only by one processor
  if(ssi_mono_->Comm().MyPID() == 0)
  {
    // extract name of quantity to be tested
    std::string quantity;
    res.ExtractString("QUANTITY",quantity);

    // get result to be tested
    const double result = ResultSpecial(quantity);

    // compare values
    const int err = CompareValues(result,"SPECIAL",res);
    nerr += err;
    ++test_count;
  }

  return;
} // SSI::SSI_Mono_ResultTest::TestSpecial
