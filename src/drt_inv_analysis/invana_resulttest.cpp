
/*----------------------------------------------------------------------*/
/*!
\file invana_resulttest.cpp

\brief tesing of inverse analysis results

<pre>
Maintainer: Sebastian Kehl
kehl@lnm.mw.tum.de
089 - 289-15249
</pre>
*/
/*----------------------------------------------------------------------*/


#include <string>

#include "invana_resulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_inv_analysis/stat_inv_analysis.H"
#include "../drt_inv_analysis/matpar_manager.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
STR::INVANA::InvAnaResultTest::InvAnaResultTest(StatInvAnalysis& ia)
  : DRT::ResultTest("INVANA")
{
    discret_= ia.Discretization();
    matman_ = ia.MatParManager();
    mysol_ = matman_->GetParams();
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void STR::INVANA::InvAnaResultTest::TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
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

      // Here we are just interested in the elements that we own (i.e. a row element)!
      if (actelement->Owner() != discret_->Comm().MyPID())
        return;

      std::string position;
      res.ExtractString("QUANTITY",position);

      int location = matman_->GetParameterLocation(element,position);


      double result = (*(*mysol_)(location))[discret_->ElementColMap()->LID(element)];

      nerr += CompareValues(result, "ELEMENT", res);
      test_count++;
    }
  }
}



