/*---------------------------------------------------------------------*/
/*! \file

\brief Result test for the adjoint fluid equations

\maintainer Martin Kronbichler

\level 3

*/
/*---------------------------------------------------------------------*/

#include <string>

#include "topopt_fluidAdjointResulttest.H"
#include "topopt_fluidAdjointImplTimeIntegration.H"
#include "../drt_lib/drt_discret.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
TOPOPT::ADJOINT::FluidAdjointResultTest::FluidAdjointResultTest(const ImplicitTimeInt& adjointfluid)
    : DRT::ResultTest("ADJOINT"),
      fluiddis_(adjointfluid.Discretization()),
      mysol_(adjointfluid.Velnp())
{
  return;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void TOPOPT::ADJOINT::FluidAdjointResultTest::TestNode(
    DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != fluiddis_->Name()) return;

  int nodeGid;
  res.ExtractInt("NODE", nodeGid);
  nodeGid -= 1;

  int havenode(fluiddis_->HaveGlobalNode(nodeGid));
  int isnodeofanybody(0);
  fluiddis_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", nodeGid + 1, fluiddis_->Name().c_str());
  }
  else
  {
    if (fluiddis_->HaveGlobalNode(nodeGid))
    {
      const DRT::Node* node = fluiddis_->gNode(nodeGid);

      // Test only, if actnode is a row node
      if (node->Owner() != fluiddis_->Comm().MyPID()) return;

      double result = 0.;

      const Epetra_BlockMap& velnpmap = mysol_->Map();

      const int numdim = DRT::Problem::Instance()->NDim();

      std::string position;
      res.ExtractString("QUANTITY", position);
      if (position == "velx")
      {
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, node, 0))];
      }
      else if (position == "vely")
      {
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, node, 1))];
      }
      else if (position == "velz")
      {
        if (numdim == 2) dserror("Cannot test result for velz in 2D case.");
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, node, 2))];
      }
      else if (position == "pressure")
      {
        if (numdim == 2)
        {
          if (fluiddis_->NumDof(0, node) < 3)
            dserror("too few dofs at node %d for pressure testing", node->Id());
          result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, node, 2))];
        }
        else
        {
          if (fluiddis_->NumDof(0, node) < 4)
            dserror("too few dofs at node %d for pressure testing", node->Id());
          result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, node, 3))];
        }
      }
      else
      {
        dserror("Quantity '%s' not supported in fluid testing", position.c_str());
      }

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}
