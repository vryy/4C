/*----------------------------------------------------------------------*/
/*!
\file fluidresulttest.cpp

\brief testing of fluid calculation results

\level 1

\maintainer Martin Kronbichler
            kronbichler@lnm.mw.tum.de
            http://www.lnm.mw.tum.de
            089 - 289-15235
*/
/*----------------------------------------------------------------------*/


#include <string>

#include "fluidresulttest.H"
#include "../drt_lib/drt_globalproblem.H"
#include "../drt_lib/drt_linedefinition.H"
#include "fluidimplicitintegration.H"
#include "fluid_utils.H"

/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
FLD::FluidResultTest::FluidResultTest(FluidImplicitTimeInt& fluid) : DRT::ResultTest("FLUID")
{
  fluiddis_ = fluid.discret_;
  mysol_ = fluid.velnp_;
  mytraction_ = fluid.stressmanager_->GetPreCalcStresses(fluid.trueresidual_);
  mywss_ = fluid.stressmanager_->GetPreCalcWallShearStresses(fluid.trueresidual_);
  myerror_ = fluid.EvaluateErrorComparedToAnalyticalSol();
  mydivu_ = fluid.EvaluateDivU();
  mydensity_scaling_ = fluid.density_scaling_;
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidResultTest::TestNode(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != fluiddis_->Name()) return;

  int node;
  res.ExtractInt("NODE", node);
  node -= 1;

  int havenode(fluiddis_->HaveGlobalNode(node));
  int isnodeofanybody(0);
  fluiddis_->Comm().SumAll(&havenode, &isnodeofanybody, 1);

  if (isnodeofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", node + 1, fluiddis_->Name().c_str());
  }
  else
  {
    if (fluiddis_->HaveGlobalNode(node))
    {
      const DRT::Node* actnode = fluiddis_->gNode(node);

      // Here we are just interested in the nodes that we own (i.e. a row node)!
      if (actnode->Owner() != fluiddis_->Comm().MyPID()) return;

      double result = 0.;

      const Epetra_BlockMap& velnpmap = mysol_->Map();

      const int numdim = DRT::Problem::Instance()->NDim();

      std::string position;
      res.ExtractString("QUANTITY", position);
      if (position == "velx")
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, actnode, 0))];
      else if (position == "vely")
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, actnode, 1))];
      else if (position == "velz")
      {
        if (numdim == 2) dserror("Cannot test result for velz in 2D case.");
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, actnode, 2))];
      }
      else if (position == "pressure")
      {
        if (fluiddis_->NumDof(0, actnode) < (numdim + 1))
          dserror("too few dofs at node %d for pressure testing", actnode->Id());
        result = (*mysol_)[velnpmap.LID(fluiddis_->Dof(0, actnode, numdim))];
      }
      else if (position == "tractionx")
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(0, actnode, 0))];
      else if (position == "tractiony")
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(0, actnode, 1))];
      else if (position == "tractionz")
      {
        if (numdim == 2) dserror("Cannot test result for tractionz in 2D case.");
        result = (*mytraction_)[(mytraction_->Map()).LID(fluiddis_->Dof(0, actnode, 2))];
      }
      else if (position == "wssx")
        result = (*mywss_)[(mywss_->Map()).LID(fluiddis_->Dof(0, actnode, 0))];
      else if (position == "wssy")
        result = (*mywss_)[(mywss_->Map()).LID(fluiddis_->Dof(0, actnode, 1))];
      else if (position == "wssz")
      {
        if (numdim == 2) dserror("Cannot test result for wssz in 2D case.");
        result = (*mywss_)[(mywss_->Map()).LID(fluiddis_->Dof(0, actnode, 2))];
      }
      else if (position == "L2errvel")
        result = (*myerror_)[0];
      else if (position == "L2errpre")
        result = (*myerror_)[1];
      else if (position == "H1errvel")
        result = (*myerror_)[2];
      else if (position == "divu")
        result = (*mydivu_);
      else
        dserror("Quantity '%s' not supported in fluid testing", position.c_str());

      nerr += CompareValues(result, "NODE", res);
      test_count++;
    }
  }
}


/*----------------------------------------------------------------------*/
/*----------------------------------------------------------------------*/
void FLD::FluidResultTest::TestElement(DRT::INPUT::LineDefinition& res, int& nerr, int& test_count)
{
  // care for the case of multiple discretizations of the same field type
  std::string dis;
  res.ExtractString("DIS", dis);
  if (dis != fluiddis_->Name()) return;

  int element;
  res.ExtractInt("ELEMENT", element);
  element -= 1;

  int haveelement(fluiddis_->HaveGlobalElement(element));
  int iselementofanybody(0);
  fluiddis_->Comm().SumAll(&haveelement, &iselementofanybody, 1);

  if (iselementofanybody == 0)
  {
    dserror("Node %d does not belong to discretization %s", element + 1, fluiddis_->Name().c_str());
  }
  else
  {
    if (fluiddis_->HaveGlobalElement(element))
    {
      const DRT::Element* actelement = fluiddis_->gElement(element);

      // Here we are just interested in the elements that we own (i.e. a row element)!
      if (actelement->Owner() != fluiddis_->Comm().MyPID()) return;

      double result = 0.;

      const Epetra_BlockMap& elerowmap = mydensity_scaling_->Map();

      std::string position;
      res.ExtractString("QUANTITY", position);
      if (position == "fluidfraction")
      {
        result = (*mydensity_scaling_)[elerowmap.LID(actelement->Id())];
      }
      else
      {
        dserror("Quantity '%s' not supported in fluid testing", position.c_str());
      }

      nerr += CompareValues(result, "ELEMENT", res);
      test_count++;
    }
  }
}
