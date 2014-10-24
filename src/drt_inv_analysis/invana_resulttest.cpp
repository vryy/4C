
/*----------------------------------------------------------------------*/
/*!
\file invana_resulttest.cpp

\brief tesing of inverse analysis results
<pre>
Maintainer: Sebastian Kehl
            kehl@mhpc.mw.tum.de
            089 - 289-10361
</pre>

*/
/*----------------------------------------------------------------------*/


#include <string>

#include "invana_resulttest.H"
#include "../drt_lib/drt_linedefinition.H"
#include "invana_base.H"
#include "matpar_manager.H"
#include "optimizer_base.H"
#include "invana_utils.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
INVANA::InvanaResultTest::InvanaResultTest(InvanaBase& ia)
  : DRT::ResultTest("INVANA"),
    ia_(ia)
{
    discret_= ia_.Discret();
    matman_ = ia_.Matman();
    mysol_ = matman_->GetMatParams();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void INVANA::InvanaResultTest::TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != discret_->Name())
    return;

  int element;
  res.ExtractInt("ELEMENT",element);
  element -= 1;

  int haveelement(discret_->HaveGlobalElement(element));
  int iselementofanybody(0);
  discret_->Comm().SumAll(&haveelement,&iselementofanybody,1);

  if (iselementofanybody==0)
  {
    dserror("Element %d does not belong to discretization %s",element+1,discret_->Name().c_str());
  }
  else
  {
    if (discret_->HaveGlobalElement(element))
    {
      const DRT::Element* actelement = discret_->gElement(element);

      // Here we are just interested in elements we own
      if (actelement->Owner() != discret_->Comm().MyPID())
        return;

      // find out which position the quantity to test has in the material parameter vector
      std::string position;
      res.ExtractString("QUANTITY",position);

      // get the value
      int location = matman_->GetParameterLocation(element,position);
      double result = (*(*mysol_)(location))[discret_->ElementColMap()->LID(element)];

      nerr += CompareValues(result, "ELEMENT", res);
      test_count++;
    }
  }
}

void INVANA::InvanaResultTest::TestSpecial(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // get the quantity to test
  std::string quantity;
  res.ExtractString("QUANTITY",quantity);

  double result=0.0;

  if (quantity == "gradient")
    INVANA::MVNorm(ia_.Optimizer()->GetGradientView(),ia_.Optimizer()->SolLayoutMapUnique(),2,&result);
  else if (quantity == "error")
    result=ia_.Optimizer()->GetObjFunctValView();
  else
    dserror("given quantity to test not to be found yet!");

  nerr += CompareValues(result, "SPECIAL", res);
  test_count++;

}
