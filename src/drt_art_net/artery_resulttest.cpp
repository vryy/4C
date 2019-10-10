/*----------------------------------------------------------------------*/
/*! \file

\brief testing of artnet calculation results

\level 3

\maintainer Johannes Kremheller
*/
/*----------------------------------------------------------------------*/

#include "artnetexplicitintegration.H"
#include "artery_resulttest.H"
#include "../drt_lib/drt_linedefinition.H"
#include "../drt_lib/drt_discret.H"
#include "art_net_impl_stationary.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ART::ArteryResultTest::ArteryResultTest(ArtNetExplicitTimeInt& art_net) : DRT::ResultTest("ARTNET")
{
  dis_ = art_net.Discretization();
  mysol_ = art_net.QAnp();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
ART::ArteryResultTest::ArteryResultTest(ArtNetImplStationary& art_net) : DRT::ResultTest("ARTNET")
{
  dis_ = art_net.Discretization();
  mysol_ = art_net.Pressurenp();
  myelevolflow_ = art_net.EleVolflow();
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ART::ArteryResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != dis_->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(dis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  dis_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", node + 1, dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalNode(node))
    {
      DRT::Node* actnode = dis_->gNode(node);

      // Strange! It seems we might actually have a global node around
      // even if it does not belong to us. But here we are just
      // interested in our nodes!
      if (actnode->Owner() != dis_->Comm().MyPID()) return;

      double result = 0.;
      const Epetra_BlockMap& pnpmap = mysol_->Map();
      std::string position;
      res.ExtractString("QUANTITY", position);

      // test result value of single scalar field
      if (position == "area")
        result = (*mysol_)[pnpmap.LID(dis_->Dof(actnode, 0))];
      else if (position == "pressure")
        result = (*mysol_)[pnpmap.LID(dis_->Dof(0, actnode, 0))];
      else if (position == "flowrate")
        result = (*mysol_)[pnpmap.LID(dis_->Dof(actnode, 1))];
      else
      {
        dserror("Quantity '%s' not supported in result-test of artery transport problems",
            position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void ART::ArteryResultTest::TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);

  if (dis != dis_->Name()) return;

  int element;
  res.ExtractInt("ELEMENT", element);
  element -= 1;

  int haveelement(dis_->HaveGlobalElement(element));
  int iselementofanybody(0);
  dis_->Comm().SumAll(&haveelement, &iselementofanybody, 1);

  if (iselementofanybody == 0)
  {
    dserror("Element %d does not belong to discretization %s", element + 1, dis_->Name().c_str());
  }
  else
  {
    if (dis_->HaveGlobalElement(element))
    {
      const DRT::Element* actelement = dis_->gElement(element);

      // Here we are just interested in the elements that we own (i.e. a row element)!
      if (actelement->Owner() != dis_->Comm().MyPID()) return;

      // extract name of quantity to be tested
      std::string quantity;
      res.ExtractString("QUANTITY", quantity);

      double result = 0.;
      // test result value of single scalar field
      if (quantity == "volflow")
      {
        if (myelevolflow_ == Teuchos::null) dserror("Element volume flow not available");
        result = (*myelevolflow_)[dis_->ElementRowMap()->LID(actelement->Id())];
      }
      else
      {
        dserror("Quantity '%s' not supported in result-test of artery transport problems",
            quantity.c_str());
      }

      nerr += CompareValues(result, "ELEMENT", res);
      test_count++;
    }
  }
}
