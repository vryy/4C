/*----------------------------------------------------------------------*/
/*!
\file red_airway_resulttest.cpp

\brief testing of Red_Airway calculation results

\maintainer Lena Yoshihara

\level 3
*/
/*----------------------------------------------------------------------*/


#include "airwayimplicitintegration.H"
#include "red_airway_resulttest.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
AIRWAY::RedAirwayResultTest::RedAirwayResultTest(RedAirwayImplicitTimeInt& airways)
  : DRT::ResultTest("RED_AIRWAY")
{
  dis_    = airways.Discretization();
  mysol_  = airways.Pnp();
  myelemsol_ = airways.AciniVolume();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void AIRWAY::RedAirwayResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS",dis);
  if (dis != dis_->Name())
    return;

  int node;
  res.ExtractInt("NODE",node);
  node -= 1;

  int havenode(dis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  dis_->Comm().SumAll(&havenode,&isnodeofanybody,1);

  if (isnodeofanybody==0)
  {
    dserror("Node %d does not belong to discretisation %s",node+1,dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalNode(node))
    {
      DRT::Node* actnode = dis_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != dis_->Comm().MyPID())
        return;

      double result = 0.;
      const Epetra_BlockMap& pnpmap = mysol_->Map();
      std::string position;
      res.ExtractString("QUANTITY",position);

      // test result value of single scalar field
      if (position=="pressure")
        result = (*mysol_)[pnpmap.LID(dis_->Dof(actnode,0))];
      // test result values for a system of scalars
      else
      {
        dserror("Quantity '%s' not supported in result-test of red_airway problems", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}


/*----------------------------------------------------------------------*
 * Element based result test for red-airway problems. Tests the results *
 * of acini_volume.                                                     *
 *                                                         roth 11/2014 *
 *----------------------------------------------------------------------*/
void AIRWAY::RedAirwayResultTest::TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
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
    dserror("Node %d does not belong to discretisation %s",element+1,dis_->Name().c_str());
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
      const Epetra_BlockMap& acvolnp_map = myelemsol_->Map();
      std::string position;
      res.ExtractString("QUANTITY",position);
      if (position=="acini_volume")
      {
        result = (*myelemsol_)[acvolnp_map.LID(actelement->Id())];
      }
      else
      {
        dserror("Quantity '%s' not supported in result-test of red_airway problems.", position.c_str());
      }

      nerr += CompareValues(result, "ELEMENT", res);
      test_count++;
    }
  }
}
