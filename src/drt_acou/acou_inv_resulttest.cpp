/*----------------------------------------------------------------------*/
/*!
\file acou_inv_resulttest.cpp

\brief testing of inverse photoacoustic reconstruction results

<pre>
        Maintainer: Svenja Schoeder
        schoeder@lnm.mw.tum.de
        http://www.lnm.mw.tum.de/staff/svenja-schoeder/
        089 - 289-15271
</pre>
*/
/*----------------------------------------------------------------------*/

#include "acou_inv_resulttest.H"
#include "acou_inv_analysis.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "../linalg/linalg_utils.H"

/*----------------------------------------------------------------------*
 |                                                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
ACOU::AcouInvResultTest::AcouInvResultTest(InvAnalysis&  invalgo)
  : DRT::ResultTest("ACOUSTIC_INVANA")
{
  dis_ = invalgo.ScatraDiscretization();
  mysol_ =  invalgo.ElementMatVec();
}

/*----------------------------------------------------------------------*
 |                                                       schoeder 01/14 |
 *----------------------------------------------------------------------*/
void ACOU::AcouInvResultTest::TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);

  if (dis != dis_->Name())
    return;

  int element;
  res.ExtractInt("ELEMENT",element);
  element -= 1;

  int haveelement(dis_->HaveGlobalElement(element));
  int iselementofanybody(0);
  dis_->Comm().SumAll(&haveelement,&iselementofanybody,1);

  if (iselementofanybody==0)
  {
    dserror("Element %d does not belong to discretization %s",element+1,dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalElement(element))
    {
      const DRT::Element* actelement = dis_->gElement(element);

      // Here we are just interested in the elements that we own (i.e. a row element)!
      if (actelement->Owner() != dis_->Comm().MyPID())
        return;

      double result = 0.;
      std::string position;
      res.ExtractString("QUANTITY",position);

      if (position == "absorptioncoeff")
      {
        result = mysol_->operator ()(0)->operator [](dis_->ElementColMap()->LID(element));
      }
      else
      {
        dserror("Quantity '%s' not supported in result-test of photoacoustic reconstruction problems", position.c_str());
      }

      nerr += CompareValues(result, "ELEMENT", res);
      test_count++;
    }
  }
}
