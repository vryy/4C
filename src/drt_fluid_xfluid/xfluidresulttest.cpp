/*----------------------------------------------------------------------*/
/*! \file

\brief xfem based fluid result tests

\level 0

 */
/*----------------------------------------------------------------------*/

#include <string>

#include "xfluid.H"
#include "xfluidfluid.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"

#include "xfluidresulttest.H"


FLD::XFluidResultTest::XFluidResultTest(const FLD::XFluid& xfluid)
    : DRT::ResultTest("XFLUID"),
      discret_(xfluid.discret_),
      velnp_(xfluid.state_->velnp_),
      node_from_zero_(true)
{
}

FLD::XFluidResultTest::XFluidResultTest(const FLD::XFluidFluid& xfluid)
    : DRT::ResultTest("XFLUID"),
      discret_(xfluid.discret_),
      velnp_(xfluid.state_->velnp_),
      coupl_discret_(xfluid.embedded_fluid_->Discretization()),
      coupl_velnp_(xfluid.embedded_fluid_->Velnp()),
      node_from_zero_(false)
{
  // Todo: remove the "node_from_zero" flag for fluidfluid:
  // adapt the test cases!
}

void FLD::XFluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);

  int node;
  res.ExtractInt("NODE", node);

  // Todo: remove!
  if (node_from_zero_) node -= 1;

  if (dis == discret_->Name())
  {
    TestNode(res, nerr, test_count, node, discret_, velnp_);
  }
  else if (dis == coupl_discret_->Name())
  {
    TestNode(res, nerr, test_count, node, coupl_discret_, coupl_velnp_);
  }
  else
    return;
}

void FLD::XFluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count,
    int node, const Teuchos::RCP<const DRT::Discretization>& discret,
    const Teuchos::RCP<const Epetra_Vector>& velnp)
{
  int havenode(discret->HaveGlobalNode(node));
  int isnodeofanybody(0);
  discret->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", node + 1, discret->Name().c_str());
  }
  else
  {
    if (discret->HaveGlobalNode(node))
    {
      DRT::Node* actnode = discret->gNode(node);

      if (actnode->Owner() != discret->Comm().MyPID()) return;

      double result = 0.;

      const Epetra_BlockMap& velnpmap = velnp->Map();

      std::string position;
      res.ExtractString("QUANTITY", position);
      if (position == "velx")
      {
        result = (*velnp)[velnpmap.LID(discret->Dof(0, actnode, 0))];
      }
      else if (position == "vely")
      {
        result = (*velnp)[velnpmap.LID(discret->Dof(0, actnode, 1))];
      }
      else if (position == "velz")
      {
        result = (*velnp)[velnpmap.LID(discret->Dof(0, actnode, 2))];
      }
      else if (position == "pressure")
      {
        result = (*velnp)[velnpmap.LID(discret->Dof(0, actnode, 3))];
      }
      else
      {
        dserror("Quantity '%s' not supported in ale testing", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}
